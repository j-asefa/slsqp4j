package sci4j.optimize.slsqp;

import sci4j.optimize.slsqp.constraints.ConstraintType;
import sci4j.optimize.slsqp.constraints.ScalarConstraint;
import sci4j.optimize.slsqp.constraints.VectorConstraint;
import sci4j.optimize.slsqp.functions.CallBackFunc;
import sci4j.optimize.slsqp.functions.Vector2ScalarFunc;
import sci4j.optimize.slsqp.functions.Vector2VectorFunc;
import sci4j.optimize.slsqp.functions.WrappedScalarFunction;
import sci4j.optimize.slsqp.util.NativeUtils;

import java.util.Arrays;
import java.util.List;

public class Slsqp implements SlsqpSolver
{
    private Vector2ScalarFunc objectiveFunc;
    private Vector2VectorFunc objectiveFuncJacobian;
    private double[] x;

    private double[] alpha = new double[]{0};
    private double[] f0 = new double[]{0};
    private double[] gs = new double[]{0};
    private double[] h1 = new double[]{0};
    private double[] h2 = new double[]{0};
    private double[] h3 = new double[]{0};
    private double[] h4 = new double[]{0};
    private double[] t = new double[]{0};
    private double[] t0 = new double[]{0};
    private double[] tol = new double[]{0};
    private int[] iexact = new int[]{0};
    private int[] incons = new int[]{0};
    private int[] ireset = new int[]{0};
    private int[] itermx = new int[]{0};
    private int[] line = new int[]{0};
    private int[] n1 = new int[]{0};
    private int[] n2 = new int[]{0};
    private int[] n3 = new int[]{0};
    private double[][] a;

    public Slsqp(Vector2ScalarFunc objectiveFunc, Vector2VectorFunc objectiveFuncJacobian, double[] x)
    {
        this.objectiveFunc = objectiveFunc;
        this.objectiveFuncJacobian = objectiveFuncJacobian;
        this.x = x;
    }

    public OptimizeResult minimize_slsqp_with_scalar_constraints(
        double[][] bounds, // bounds is 2xn matrix. bounds[0] is array of lower bounds, bounds[1] is array of upper bounds;
        List<ScalarConstraint> scalarConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc,
        double... objectiveFuncionArgs)
    {
        final WrappedScalarFunction wrappedObjectiveFunction = new WrappedScalarFunction(objectiveFunc, objectiveFuncionArgs);
        final int n = x.length;

        int n_1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        int[] mode = new int[]{0};

        int m = scalarConstraintList.size();
        int meq = (int)scalarConstraintList.stream().filter(p -> p.getConstraintType() == ConstraintType.EQ).count();

        int mineq = m - meq + n_1 + n_1;
        int l_w = (3 * n_1 + m) * (n_1 + 1) +
            (n_1 - meq + 1) * (mineq + 2) + 2 * mineq +
            (n_1 + mineq) * (n_1 - meq) + 2 * meq +
            n_1 * n / 2 + 2 * m + 3 * n + 4 * n_1 + 1;

        int l_jw = mineq;

        double[] w = new double[l_w];
        Arrays.fill(w, 0);
        int[] jw = new int[l_jw];
        Arrays.fill(jw, 0);

        final int iter = maxIterations;

        int[] majiter = new int[] {iter};
        int majiter_prev = 0;
        final double[] acc = new double[] {tolerance};

        double fx = 0;
        int la = Math.max(1, m);
        double[] c = new double[la];

        double[] g = new double[n_1];

        // Note that Fortran expects arrays to be laid out in column-major order.
        double[][] a = new double[n_1][la];
        double[] fprime;

        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                // calculate the value of the objective function
                fx = wrappedObjectiveFunction.apply(x);

                int constraintIndex = 0;
                // first get the equality constraints
                for (ScalarConstraint constraint : scalarConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double constraintVal = constraint.apply(x);
                        c[constraintIndex] = constraintVal;
                        constraintIndex++;
                    }
                }
                // then get the inequality constraints
                for (ScalarConstraint constraint : scalarConstraintList)
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
                if (objectiveFuncJacobian == null)
                {
                    fprime = wrappedObjectiveFunction.approx_jacobian(x);
                }
                else
                {
                    fprime = objectiveFuncJacobian.apply(x, objectiveFuncionArgs);
                }

                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                int i = 0;
                for (ScalarConstraint constraint : scalarConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[] constraintJac = constraint.getJacobian(x);

                        assert (constraintJac.length <= n);
                        // copy the constraint jacobian matrix into the array a
                        // jacDim2 is the number of columns in constraint jacobian
                        for (int j = 0; j < constraintJac.length; j++)
                        {
                            a[j][i] = constraintJac[j];
                        }
                        i++;
                    }
                }

                for (ScalarConstraint constraint : scalarConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[] constraintJac = constraint.getJacobian(x);
                        assert (constraintJac.length <= n);

                        // copy the constraint jacobian matrix into the array a
                        // jacDim2 is the number of columns in constraint jacobian
                        for (int j = 0; j < constraintJac.length; j++)
                        {
                            a[j][i] = constraintJac[j];
                        }
                        i++;
                    }
                }
            }

            NativeUtils.slsqp(m, meq, la, x, xl, xu, new double[] {fx}, c, g, a, acc, majiter, mode, w, jw,
                alpha, f0, gs, h1, h2, h3, h4, t, t0, tol,
                iexact, incons, ireset, itermx, line,
                n1, n2, n3);

            if (callBackFunc != null && majiter[0] > majiter_prev)
            {
                callBackFunc.callback(Arrays.copyOf(x, x.length));
            }

            if (Math.abs(mode[0]) != 1)
            {
                break;
            }

            majiter_prev = majiter[0];
        }
        return new OptimizeResult(x, fx, g, mode[0], majiter[0], mode[0], mode[0] == 0, a);
    }

    public OptimizeResult minimize_slsqp_with_vector_constraints(
        double[][] bounds, // bounds is 2xn matrix. bounds[0] is array of lower bounds, bounds[1] is array of upper bounds;
        List<VectorConstraint> vectorConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc,
        double... objectiveFuncionArgs
        )
    {
        final WrappedScalarFunction wrappedObjectiveFunction = new WrappedScalarFunction(objectiveFunc, objectiveFuncionArgs);
        final int n = x.length;

        int n_1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        int[] mode = new int[] {0};

        int m = 0;
        int meq = 0;

        // get the number of constraints
        for (VectorConstraint constraint : vectorConstraintList)
        {
            m += constraint.apply(x).length;
            if (constraint.getConstraintType() == ConstraintType.EQ)
            {
                meq++;
            }
        }
        int la = Math.max(1, m);

        int mineq = m - meq + n_1 + n_1;

        int l_w = (3 * n_1 + m) * (n_1 + 1) +
            (n_1 - meq + 1) * (mineq + 2) + 2 * mineq +
            (n_1 + mineq) * (n_1 - meq) + 2 * meq +
            n_1 * n / 2 + 2 * m + 3 * n + 4 * n_1 + 1;

        int l_jw = mineq;

        double[] w = new double[l_w];
        int[] jw = new int[l_jw];

        final int iter = maxIterations;

        int[] majiter = new int[] {iter};
        int majiter_prev = 0;
        final double[] acc = new double[] {tolerance};

        double fx = 0;

        double[] g = new double[n_1];

        double[] c = new double[la];
        // Note that Fortran expects arrays to be laid out in column-major order.
        double[][] a = new double[n_1][la];

        double[] fprime = new double[n];
        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                // calculate the value of the objective function
                fx = wrappedObjectiveFunction.apply(x);

                int constraintIndex = 0;
                // first get the equality constraints
                for (VectorConstraint constraint : vectorConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[] constraintVec = constraint.apply(x);
                        System.arraycopy(constraintVec, 0, c, constraintIndex, constraintVec.length);
                        constraintIndex += constraintVec.length;
                    }
                }
                // then get the inequality constraints
                for (VectorConstraint constraint : vectorConstraintList)
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
                if (objectiveFuncJacobian == null)
                {
                    fprime = wrappedObjectiveFunction.approx_jacobian(x);
                }
                else
                {
                    fprime = objectiveFuncJacobian.apply(x, objectiveFuncionArgs);
                }
                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                int i = 0;
                for (VectorConstraint constraint : vectorConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[][] constraintJac = constraint.getJacobian(x);
                        assert (constraintJac.length <= n);

                        // copy the constraint jacobian matrix into the array a
                        // jacDim2 is the number of columns in constraint jacobian
                        final int jacDim2 = constraintJac[0].length;
                        for (int l = 0; l < jacDim2; l++)
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
                        }
                        i += jacDim2;
                    }
                }

                for (VectorConstraint constraint : vectorConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[][] constraintJac = constraint.getJacobian(x);
                        assert (constraintJac.length <= n);

                        // copy the constraint jacobian matrix into the array a
                        // jacDim2 is the number of columns in constraint jacobian
                        final int jacDim2 = constraintJac[0].length;
                        for (int l = 0; l < jacDim2; l++)
                        {
                            for (int j = 0; j < constraintJac.length; j++)
                            {
                                if (constraintJac[j].length > 1)
                                {
                                    a[j][i] = constraintJac[j][l]; // a is in column-major order
                                } else
                                {
                                    a[j][i] = constraintJac[j][0];
                                }
                            }
                        }
                        i += jacDim2;
                    }
                }
            }

            this.a = a;

            NativeUtils.slsqp(m, meq, la, x, xl, xu, new double[] {fx}, c, g, a, acc, majiter, mode, w, jw,
                alpha, f0, gs, h1, h2, h3, h4, t, t0, tol,
                iexact, incons, ireset, itermx, line,
                n1, n2, n3);

            if (callBackFunc != null && majiter[0] > majiter_prev)
            {
                callBackFunc.callback(Arrays.copyOf(x, x.length));
            }

            if (Math.abs(mode[0]) != 1)
            {
                break;
            }

            majiter_prev = majiter[0];
        }
        return new OptimizeResult(x, fx, g, mode[0], majiter[0], mode[0], mode[0] == 0, a);
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
