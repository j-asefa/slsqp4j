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

public class Slsqp
{
    /**
     *
     * @param objectiveFunc objective function of the optimization problem
     * @param objectiveFuncJacobian optional function representing the Jacobian of the objective function. If null,
     *                              numerical approximation is used.
     * @param x input to the optimization problem.
     * @param bounds 2xn matrix. bounds[0] is array of lower bounds, bounds[1] is array of upper bounds
     * @param scalarConstraintList list of scalar-valued constraints
     * @param tolerance accuracy tolerance
     * @param maxIterations max number of iterations to run the solver for
     * @param callBackFunc optional callback function on major iterations. If null, no callback is called.
     * @param objectiveFuncionArgs optional arguments to the objective function
     * @return OptimizeResult result of optimization
     */
    public static OptimizeResult minimizeSlsqpWithScalarConstraints(
        Vector2ScalarFunc objectiveFunc,
        Vector2VectorFunc objectiveFuncJacobian,
        double[] x,
        double[][] bounds,
        List<ScalarConstraint> scalarConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc,
        double... objectiveFuncionArgs)
    {
        final double[] alpha = new double[]{0};
        final double[] f0 = new double[]{0};
        final double[] gs = new double[]{0};
        final double[] h1 = new double[]{0};
        final double[] h2 = new double[]{0};
        final double[] h3 = new double[]{0};
        final double[] h4 = new double[]{0};
        final double[] t = new double[]{0};
        final double[] t0 = new double[]{0};
        final double[] tol = new double[]{0};
        final int[] iexact = new int[]{0};
        final int[] incons = new int[]{0};
        final int[] ireset = new int[]{0};
        final int[] itermx = new int[]{0};
        final int[] line = new int[]{0};
        final int[] n1 = new int[]{0};
        final int[] n2 = new int[]{0};
        final int[] n3 = new int[]{0};

        final WrappedScalarFunction wrappedObjectiveFunction =
            new WrappedScalarFunction(objectiveFunc, objectiveFuncionArgs);

        final int n = x.length;

        final int nPlus1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        final int[] mode = new int[]{0};

        final int m = scalarConstraintList.size();
        final int meq = (int)scalarConstraintList.stream().filter(
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
                for (final ScalarConstraint constraint : scalarConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double constraintVal = constraint.apply(x);
                        c[constraintIndex] = constraintVal;
                        constraintIndex++;
                    }
                }
                // then get the inequality constraints
                for (final ScalarConstraint constraint : scalarConstraintList)
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
                    fprime = wrappedObjectiveFunction.approxJacobian(x);
                }
                else
                {
                    fprime = objectiveFuncJacobian.apply(x, objectiveFuncionArgs);
                }

                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                int i = 0;
                for (final ScalarConstraint constraint : scalarConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[] constraintJac = constraint.getJacobian(x);

                        // copy the constraint jacobian matrix into the array a
                        // jacDim2 is the number of columns in constraint jacobian
                        for (int j = 0; j < constraintJac.length; j++)
                        {
                            a[j][i] = constraintJac[j];
                        }
                        i++;
                    }
                }

                for (final ScalarConstraint constraint : scalarConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[] constraintJac = constraint.getJacobian(x);

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
        return new OptimizeResult(x, fx, g, mode[0], majIter[0], mode[0], mode[0] == 0, a);
    }

    /**
     *
     * @param objectiveFunc objective function of the optimization problem
     * @param objectiveFuncJacobian optional function representing the Jacobian of the objective function. If null,
     *                              numerical approximation is used.
     * @param x input to the optimization problem.
     * @param bounds 2xn matrix. bounds[0] is array of lower bounds, bounds[1] is array of upper bounds
     * @param vectorConstraintList list of vector-valued constraints
     * @param tolerance accuracy tolerance
     * @param maxIterations max number of iterations to run the solver for
     * @param callBackFunc optional callback function on major iterations. If null, no callback is called.
     * @param objectiveFuncionArgs optional arguments to the objective function
     * @return OptimizeResult result of optimization
     */
    public static OptimizeResult minimizeSlsqpWithVectorConstraints(
        Vector2ScalarFunc objectiveFunc,
        Vector2VectorFunc objectiveFuncJacobian,
        double[] x,
        double[][] bounds,
        List<VectorConstraint> vectorConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc,
        double... objectiveFuncionArgs)
    {
        final double[] alpha = new double[]{0};
        final double[] f0 = new double[]{0};
        final double[] gs = new double[]{0};
        final double[] h1 = new double[]{0};
        final double[] h2 = new double[]{0};
        final double[] h3 = new double[]{0};
        final double[] h4 = new double[]{0};
        final double[] t = new double[]{0};
        final double[] t0 = new double[]{0};
        final double[] tol = new double[]{0};
        final int[] iexact = new int[]{0};
        final int[] incons = new int[]{0};
        final int[] ireset = new int[]{0};
        final int[] itermx = new int[]{0};
        final int[] line = new int[]{0};
        final int[] n1 = new int[]{0};
        final int[] n2 = new int[]{0};
        final int[] n3 = new int[]{0};

        final WrappedScalarFunction wrappedObjectiveFunction =
            new WrappedScalarFunction(objectiveFunc, objectiveFuncionArgs);

        final int n = x.length;

        final int nPlus1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        final int[] mode = new int[] {0};

        int m = 0;
        int meq = 0;

        // get the number of constraints
        for (final VectorConstraint constraint : vectorConstraintList)
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
                for (final VectorConstraint constraint : vectorConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[] constraintVec = constraint.apply(x);
                        System.arraycopy(constraintVec, 0, c, constraintIndex, constraintVec.length);
                        constraintIndex += constraintVec.length;
                    }
                }
                // then get the inequality constraints
                for (final VectorConstraint constraint : vectorConstraintList)
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
                    fprime = wrappedObjectiveFunction.approxJacobian(x);
                }
                else
                {
                    fprime = objectiveFuncJacobian.apply(x, objectiveFuncionArgs);
                }
                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                int i = 0;
                for (final VectorConstraint constraint : vectorConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.EQ)
                    {
                        final double[][] constraintJac = constraint.getJacobian(x);

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
                            i++;
                        }
                    }
                }
                for (final VectorConstraint constraint : vectorConstraintList)
                {
                    if (constraint.getConstraintType() == ConstraintType.INEQ)
                    {
                        final double[][] constraintJac = constraint.getJacobian(x);

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
        return new OptimizeResult(x, fx, g, mode[0], majIter[0], mode[0], mode[0] == 0, a);
    }

    private static void computeBounds(int n, double[][] bounds, double[] xl, double[] xu)
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

    private static void clip(double[] x, double[] xl, double[] xu)
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
