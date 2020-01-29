package slsqp.optimize;

import slsqp.optimize.constraints.ConstraintType;
import slsqp.optimize.constraints.ScalarConstraint;
import slsqp.optimize.constraints.VectorConstraint;
import slsqp.optimize.functions.CallBackFunc;
import slsqp.optimize.functions.Vector2ScalarFunc;
import slsqp.optimize.functions.Vector2VectorFunc;
import slsqp.optimize.functions.WrappedVector2ScalarFunction;
import slsqp.optimize.util.NativeUtils;

import java.util.Arrays;

public final class Slsqp
{
    public static class SlsqpBuilder
    {
        private Vector2ScalarFunc objectiveFunc;
        private Vector2VectorFunc objectiveFuncJacobian;
        private double[][] bounds;
        private ScalarConstraint[] scalarConstraints;
        private double tolerance;
        private int maxIterations;
        private CallBackFunc callBackFunc;
        private double[] objectiveFunctionArgs;
        private VectorConstraint[] vectorConstraints;

        public SlsqpBuilder withObjectiveFunction(Vector2ScalarFunc objectiveFunc, double... objectiveFunctionArgs)
        {
            this.objectiveFunc = objectiveFunc;
            this.objectiveFunctionArgs = objectiveFunctionArgs;
            return this;
        }

        public SlsqpBuilder withJacobian(Vector2VectorFunc objectiveFuncJacobian)
        {
            this.objectiveFuncJacobian = objectiveFuncJacobian;
            return this;
        }

        public SlsqpBuilder withBounds(double[][] bounds)
        {
            this.bounds = bounds;
            return this;
        }

        public SlsqpBuilder addScalarConstraint(ScalarConstraint[] scalarConstraints)
        {
            this.scalarConstraints = scalarConstraints;
            return this;
        }

        public SlsqpBuilder addVectorConstraint(VectorConstraint[] vectorConstraints)
        {
            this.vectorConstraints = vectorConstraints;
            return this;
        }

        public SlsqpBuilder withTolerance(double tolerance)
        {
            this.tolerance = tolerance;
            return this;
        }

        public SlsqpBuilder withMaxIterations(int maxIterations)
        {
            this.maxIterations = maxIterations;
            return this;
        }

        public SlsqpBuilder withCallBackFunction(CallBackFunc callBackFunc)
        {
            this.callBackFunc = callBackFunc;
            return this;
        }

        public Slsqp build()
        {
            final Slsqp slsqp = new Slsqp();
            slsqp.objectiveFunc = this.objectiveFunc;
            slsqp.objectiveFuncJacobian = this.objectiveFuncJacobian;
            slsqp.bounds = this.bounds;
            slsqp.scalarConstraints = this.scalarConstraints;
            slsqp.tolerance = this.tolerance;
            slsqp.maxIterations = this.maxIterations;
            slsqp.callBackFunc = this.callBackFunc;
            slsqp.objectiveFunctionArgs = this.objectiveFunctionArgs;
            slsqp.vectorConstraints = this.vectorConstraints;
            return slsqp;
        }
    }

    private Slsqp()
    {
    }

    private Vector2ScalarFunc objectiveFunc;
    private Vector2VectorFunc objectiveFuncJacobian;
    private double[][] bounds;
    private ScalarConstraint[] scalarConstraints;
    private double tolerance;
    private int maxIterations;
    private CallBackFunc callBackFunc;
    private double[] objectiveFunctionArgs;
    private VectorConstraint[] vectorConstraints;

    private final double[] alpha = new double[]{0};
    private final double[] f0 = new double[]{0};
    private final double[] gs = new double[]{0};
    private final double[] h1 = new double[]{0};
    private final double[] h2 = new double[]{0};
    private final double[] h3 = new double[]{0};
    private final double[] h4 = new double[]{0};
    private final double[] t = new double[]{0};
    private final double[] t0 = new double[]{0};
    private final double[] tol = new double[]{0};
    private final int[] iexact = new int[]{0};
    private final int[] incons = new int[]{0};
    private final int[] ireset = new int[]{0};
    private final int[] itermx = new int[]{0};
    private final int[] line = new int[]{0};
    private final int[] n1 = new int[]{0};
    private final int[] n2 = new int[]{0};
    private final int[] n3 = new int[]{0};
    private final int[] mode = new int[]{0};

    public OptimizeResult optimize(double[] x)
    {
        if (this.vectorConstraints == null)
        {
            return optimizeWithScalarConstraints(x);
        }
        else
        {
            return optimizeWithVectorConstraints(x);
        }
    }

    private OptimizeResult optimizeWithVectorConstraints(double[] x)
    {
        final WrappedVector2ScalarFunction wrappedObjectiveFunction =
            new WrappedVector2ScalarFunction(objectiveFunc, objectiveFuncJacobian, objectiveFunctionArgs);

        final int n = x.length;

        final int nPlus1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        int m = 0;
        int meq = 0;

        // get the number of constraints
        for (final VectorConstraint constraint : vectorConstraints)
        {
            m += constraint.apply(x).length;
            if (constraint.getConstraintType() == ConstraintType.EQ)
            {
                meq++;
            }
        }
        final int la = Math.max(1, m);

        final int mineq = m - meq + nPlus1 + nPlus1;

        final int lenW = (3 * nPlus1 + m) * (nPlus1 + 1) +
            (nPlus1 - meq + 1) * (mineq + 2) + 2 * mineq +
            (nPlus1 + mineq) * (nPlus1 - meq) + 2 * meq +
            nPlus1 * n / 2 + 2 * m + 3 * n + 4 * nPlus1 + 1;

        final double[] w = new double[lenW];
        final int[] jw = new int[mineq];

        final int[] majIter = new int[] {maxIterations};
        int majIterPrev = 0;
        final double[] acc = new double[] {tolerance};

        double fx = 0;

        final double[] g = new double[nPlus1];

        final double[] c = new double[la];

        // Note that Fortran expects arrays to be laid out in column-major order.
        final double[][] a = new double[nPlus1][la];

        double[] fprime;
        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                // calculate the value of the objective function
                fx = wrappedObjectiveFunction.apply(x);

                int constraintIndex = 0;
                // first get the equality constraints
                for (final VectorConstraint constraint : vectorConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[] constraintVec = constraint.apply(x);
                        System.arraycopy(constraintVec, 0, c, constraintIndex, constraintVec.length);
                        constraintIndex += constraintVec.length;
                    }
                }
                // then get the inequality constraints
                for (final VectorConstraint constraint : vectorConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[] constraintVec = constraint.apply(x);
                        System.arraycopy(constraintVec, 0, c, constraintIndex, constraintVec.length);
                        constraintIndex += constraintVec.length;
                    }
                }
            }

            if (mode[0] == 0 || mode[0] == -1)
            {
                fprime = wrappedObjectiveFunction.getJacobian(x);
                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                int i = 0;
                for (final VectorConstraint constraint : vectorConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[][] constraintJac = constraint.getJacobian(x);

                        // copy the constraint jacobian matrix into the array a
                        for (int l = 0; l < constraintJac[0].length; l++)
                        {
                            for (int j = 0; j < constraintJac.length; j++)
                            {
                                if (constraintJac[j].length > 1)
                                {
                                    a[j][i] = constraintJac[j][l]; // a is in column-major order
                                }
                                else
                                {
                                    a[j][i] = constraintJac[j][0];
                                }
                            }
                            i++;
                        }
                    }
                }
                for (final VectorConstraint constraint : vectorConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[][] constraintJac = constraint.getJacobian(x);

                        // copy the constraint jacobian matrix into the array a
                        for (int l = 0; l < constraintJac[0].length; l++)
                        {
                            for (int j = 0; j < constraintJac.length; j++)
                            {
                                if (constraintJac[j].length > 1)
                                {
                                    a[j][i] = constraintJac[j][l]; // a is in column-major order
                                }
                                else
                                {
                                    a[j][i] = constraintJac[j][0];
                                }
                            }
                            i++;
                        }
                    }
                }
            }

            NativeUtils.slsqp(m, meq, la, x, xl, xu, new double[] {fx}, c, g, a, acc, majIter, mode, w, jw,
                alpha, f0, gs, h1, h2, h3, h4, t, t0, tol,
                iexact, incons, ireset, itermx, line,
                n1, n2, n3);

            if (callBackFunc != null && majIter[0] > majIterPrev)
            {
                callBackFunc.callback(Arrays.copyOf(x, x.length));
            }

            if (Math.abs(mode[0]) != 1)
            {
                break;
            }

            majIterPrev = majIter[0];
        }
        return new OptimizeResult(x, fx, g, majIter[0], mode[0], mode[0] == 0, a);
    }

    private OptimizeResult optimizeWithScalarConstraints(double[] x)
    {
        final WrappedVector2ScalarFunction wrappedObjectiveFunction =
            new WrappedVector2ScalarFunction(objectiveFunc, objectiveFuncJacobian, objectiveFunctionArgs);

        final int n = x.length;

        final int nPlus1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        final int m = scalarConstraints.length;
        final int meq = (int)Arrays.stream(scalarConstraints).filter(
            p -> p.getConstraintType() == ConstraintType.EQ
        ).count();

        final int mineq = m - meq + nPlus1 + nPlus1;
        final int lenW = (3 * nPlus1 + m) * (nPlus1 + 1) +
            (nPlus1 - meq + 1) * (mineq + 2) + 2 * mineq +
            (nPlus1 + mineq) * (nPlus1 - meq) + 2 * meq +
            nPlus1 * n / 2 + 2 * m + 3 * n + 4 * nPlus1 + 1;

        final double[] w = new double[lenW];
        final int[] jw = new int[mineq];

        final int[] majIter = new int[] {maxIterations};
        int majIterPrev = 0;
        final double[] acc = new double[] {tolerance};

        double fx = 0;
        final int la = Math.max(1, m);
        final double[] c = new double[la];

        final double[] g = new double[nPlus1];

        // Note that Fortran expects arrays to be laid out in column-major order.
        final double[][] a = new double[nPlus1][la];
        double[] fprime;

        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                // calculate the value of the objective function
                fx = wrappedObjectiveFunction.apply(x);

                int constraintIndex = 0;
                // first get the equality constraints
                for (final ScalarConstraint constraint : scalarConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double constraintVal = constraint.apply(x);
                        c[constraintIndex] = constraintVal;
                        constraintIndex++;
                    }
                }
                // then get the inequality constraints
                for (final ScalarConstraint constraint : scalarConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double constraintVal = constraint.apply(x);
                        c[constraintIndex] = constraintVal;
                        constraintIndex++;
                    }
                }
            }

            if (mode[0] == 0 || mode[0] == -1)
            {
                fprime = wrappedObjectiveFunction.getJacobian(x);

                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                int i = 0;
                for (final ScalarConstraint constraint : scalarConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[] constraintJac = constraint.getJacobian(x);

                        // copy the constraint jacobian matrix into the array a
                        for (int j = 0; j < constraintJac.length; j++)
                        {
                            a[j][i] = constraintJac[j];
                        }
                        i++;
                    }
                }

                for (final ScalarConstraint constraint : scalarConstraints)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[] constraintJac = constraint.getJacobian(x);

                        // copy the constraint jacobian matrix into the array a
                        for (int j = 0; j < constraintJac.length; j++)
                        {
                            a[j][i] = constraintJac[j];
                        }
                        i++;
                    }
                }
            }

            NativeUtils.slsqp(m, meq, la, x, xl, xu, new double[] {fx}, c, g, a, acc, majIter, mode, w, jw,
                alpha, f0, gs, h1, h2, h3, h4, t, t0, tol,
                iexact, incons, ireset, itermx, line,
                n1, n2, n3);

            if (callBackFunc != null && majIter[0] > majIterPrev)
            {
                callBackFunc.callback(Arrays.copyOf(x, x.length));
            }

            if (Math.abs(mode[0]) != 1)
            {
                break;
            }

            majIterPrev = majIter[0];
        }
        return new OptimizeResult(x, fx, g, majIter[0], mode[0], mode[0] == 0, a);
    }

    private void computeBounds(int n, double[][] bounds, double[] xl, double[] xu)
    {
        if (bounds == null)
        {
            Arrays.fill(xl, Double.NaN);
            Arrays.fill(xu, Double.NaN);
        }
        else
        {
            System.arraycopy(bounds[0], 0, xl, 0, n);
            System.arraycopy(bounds[1], 0, xu, 0, n);
        }
        for (int i = 0; i < n; i++)
        {
            if (xl[i] == Double.NEGATIVE_INFINITY || xl[i] == Double.POSITIVE_INFINITY)
            {
                xl[i] = Double.NaN;
            }
            if (xu[i] == Double.NEGATIVE_INFINITY || xu[i] == Double.POSITIVE_INFINITY)
            {
                xu[i] = Double.NaN;
            }
        }

    }

    private void clip(double[] x, double[] xl, double[] xu)
    {
        final int n = x.length;
        for (int i = 0; i < n; i++)
        {
            if (x[i] < xl[i])
            {
                x[i] = xl[i];
            }
            else if (x[i] > xu[i])
            {
                x[i] = xu[i];
            }
            else if (x[i] >= Double.POSITIVE_INFINITY)
            {
                x[i] = Double.POSITIVE_INFINITY;
            }
            else if (x[i] <= Double.NEGATIVE_INFINITY)
            {
                x[i] = Double.NEGATIVE_INFINITY;
            }
        }
    }
}
