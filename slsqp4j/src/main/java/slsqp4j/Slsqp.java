/*
Copyright (c) 1988 Dieter Kraft

Copyright (c) 1994 Association for Computing Machinery

Copyright (c) 2001, 2002 Enthought, Inc.
All rights reserved.

Copyright (c) 2003-2019 SciPy Developers.
All rights reserved.

Copyright (c) 2020, Skew Ltd.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package slsqp4j;

import slsqp4j.constraints.ConstraintType;
import slsqp4j.constraints.ScalarConstraint;
import slsqp4j.constraints.VectorConstraint;
import slsqp4j.functions.CallBackFunc;
import slsqp4j.functions.Vector2ScalarFunc;
import slsqp4j.functions.Vector2VectorFunc;
import slsqp4j.functions.WrappedVector2ScalarFunction;
import slsqp4j.util.NativeUtils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * This module implements the Sequential Least Squares Programming optimization
 * algorithm (SLSQP), originally developed by Dieter Kraft.
 * See http://www.netlib.org/toms/733
 *
 * An instance of {@link Slsqp} corresponds to an instance of the Slsqp solver. A reference to an Slsqp instance
 * is constructed using the {@link SlsqpBuilder} in order to specify parameters, some of which may be optional.
 *
 * Once an instance of {@link Slsqp} is constructed, calls to {@link #minimize(double[])} can occur repeatedly. It is
 * expected that users will call {@link #minimize(double[])} iteratively within their own code until some error level is
 * satisfied or a value of true is returned by {@link OptimizeResult#success()}.
 *
 * Instances of {@link Slsqp} are not thread safe.
 */
public final class Slsqp
{
    public static class SlsqpBuilder
    {
        private static final double DEFAULT_ACCURACY = 1.0E-6;
        private static final int DEFAULT_MAX_ITER = 100;

        private Vector2ScalarFunc objectiveFunc;
        private Vector2VectorFunc objectiveFuncJacobian;
        private double[][] bounds;
        private Set<ScalarConstraint> scalarEqualityConstraints = new HashSet<>();
        private Set<ScalarConstraint> scalarInequalityConstraints = new HashSet<>();
        private Set<VectorConstraint> vectorEqualityConstraints = new HashSet<>();
        private Set<VectorConstraint> vectorInequalityConstraints = new HashSet<>();
        private double accuracy = DEFAULT_ACCURACY;
        private int maxIterations = DEFAULT_MAX_ITER;
        private CallBackFunc callBackFunc;
        private double[] objectiveFunctionArgs;

        /**
         * Set the objective function for this problem to minimize.
         *
         * @param objectiveFunc objective function.
         * @param objectiveFunctionArgs arguments, if any, to the objective function.
         * @return this builder.
         */
        public SlsqpBuilder withObjectiveFunction(Vector2ScalarFunc objectiveFunc, double... objectiveFunctionArgs)
        {
            this.objectiveFunc = objectiveFunc;
            this.objectiveFunctionArgs = objectiveFunctionArgs;
            return this;
        }

        /**
         * Set the analytical Jacobian of this objective function, if any. If no analytical Jacobian is specified,
         * numerical approximation is performed during optimization.
         *
         * @param objectiveFuncJacobian analytical jacobian of the objective function specified in
         * {@link #withObjectiveFunction(Vector2ScalarFunc, double...)}
         * @return this builder.
         */
        public SlsqpBuilder withJacobian(Vector2VectorFunc objectiveFuncJacobian)
        {
            this.objectiveFuncJacobian = objectiveFuncJacobian;
            return this;
        }

        /**
         *
         *
         * @param bounds
         * @return
         */
        public SlsqpBuilder withBounds(double[][] bounds)
        {
            this.bounds = bounds;
            return this;
        }

        /**
         * Add a {@link ScalarConstraint} to this builder.
         *
         * @param constraint scalar constraint to include in this optimization problem.
         * @return this builder.
         */
        public SlsqpBuilder addScalarConstraint(ScalarConstraint constraint)
        {
            if (constraint.getConstraintType() == ConstraintType.EQ)
            {
                this.scalarEqualityConstraints.add(constraint);
            }
            else
            {
                this.scalarInequalityConstraints.add(constraint);
            }
            return this;
        }

        /**
         * Add a {@link VectorConstraint} to this builder.
         *
         * @param constraint vector constraint to include in this optimization problem.
         * @return this builder.
         */
        public SlsqpBuilder addVectorConstraint(VectorConstraint constraint)
        {
            if (constraint.getConstraintType() == ConstraintType.EQ)
            {
                this.vectorEqualityConstraints.add(constraint);
            }
            else
            {
                this.vectorInequalityConstraints.add(constraint);
            }
            return this;
        }

        /**
         * Set the desired accuracy for this problem. Defaults to {@link #DEFAULT_ACCURACY}.
         *
         * @param accuracy desired accuracy for the minimization problem.
         * @return this builder.
         */
        public SlsqpBuilder withAccuracy(double accuracy)
        {
            this.accuracy = accuracy;
            return this;
        }

        /**
         * Set the maximum number of iterations to perform.
         *
         * @param maxIterations maximum number of iterations to perform.
         * @return this builder.
         */
        public SlsqpBuilder withMaxIterations(int maxIterations)
        {
            this.maxIterations = maxIterations;
            return this;
        }

        /**
         * Set a callback function to be called on every major iteration of the solver.
         *
         * @param callBackFunc callback function to call on major iterations.
         * @return this builder.
         */
        public SlsqpBuilder withCallBackFunction(CallBackFunc callBackFunc)
        {
            this.callBackFunc = callBackFunc;
            return this;
        }

        /**
         * Build an instance of a {@link Slsqp}.
         *
         * @return a new {@link Slsqp} with the properties specified in this builder.
         */
        public Slsqp build()
        {
            if (this.objectiveFunc == null)
            {
                throw new IllegalStateException("must specify an objective function");
            }
            if (!this.scalarEqualityConstraints.isEmpty() &&
                (!this.vectorEqualityConstraints.isEmpty() || !this.vectorInequalityConstraints.isEmpty()))
            {
                throw new IllegalStateException("cannot specify both vector and scalar constraints");
            }
            else if (!this.scalarInequalityConstraints.isEmpty() &&
                (!this.vectorEqualityConstraints.isEmpty() || !this.vectorInequalityConstraints.isEmpty()))
            {
                throw new IllegalStateException("cannot specify both vector and scalar constraints");
            }
            final Slsqp slsqp = new Slsqp();
            slsqp.bounds = this.bounds;
            slsqp.scalarEqualityConstraints = this.scalarEqualityConstraints;
            slsqp.scalarInequalityConstraints = this.scalarInequalityConstraints;
            slsqp.vectorEqualityConstraints = this.vectorEqualityConstraints;
            slsqp.vectorInequalityConstraints = this.vectorInequalityConstraints;
            slsqp.tolerance = this.accuracy;
            slsqp.maxIterations = this.maxIterations;
            slsqp.callBackFunc = this.callBackFunc;
            slsqp.wrappedObjectiveFunction =
                new WrappedVector2ScalarFunction(objectiveFunc, objectiveFuncJacobian, objectiveFunctionArgs);

            slsqp.majIter = new int[] {maxIterations};
            slsqp.acc = new double[] {accuracy};
            return slsqp;
        }

        public void reset()
        {
            this.objectiveFunc = null;
            this.objectiveFuncJacobian = null;
            this.bounds = null;
            this.scalarEqualityConstraints.clear();
            this.scalarInequalityConstraints.clear();
            this.vectorEqualityConstraints.clear();
            this.vectorInequalityConstraints.clear();
            this.accuracy = DEFAULT_ACCURACY;
            this.maxIterations = DEFAULT_MAX_ITER;
            this.callBackFunc = null;
            this.objectiveFunctionArgs = null;
        }
    }

    private Slsqp()
    {
    }

    private double[][] bounds;
    private Set<ScalarConstraint> scalarEqualityConstraints;
    private Set<ScalarConstraint> scalarInequalityConstraints;
    private Set<VectorConstraint> vectorEqualityConstraints;
    private Set<VectorConstraint> vectorInequalityConstraints;
    private double tolerance;
    private int maxIterations;
    private CallBackFunc callBackFunc;
    private WrappedVector2ScalarFunction wrappedObjectiveFunction;
    private int[] majIter = new int[] {maxIterations};
    private double[] acc = new double[] {tolerance};

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

    /**
     * Minimize the objective function specified in this instance of {@link Slsqp} subject to the specified constraints.
     *
     * @param x initial guess for the independent variables of the problem.
     * @return an instance of {@link OptimizeResult} with the current state of the problem.
     */
    public OptimizeResult minimize(double[] x)
    {
        if (this.vectorEqualityConstraints.isEmpty() && this.vectorInequalityConstraints.isEmpty())
        {
            return minimizeWithScalarConstraints(x);
        }
        else
        {
            return minimizeWithVectorConstraints(x);
        }
    }

    private OptimizeResult minimizeWithVectorConstraints(double[] x)
    {
        final int n = x.length;

        final int nPlus1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        final int meq = vectorEqualityConstraints.size();

        int m = 0;
        // get the number of constraints
        for (final VectorConstraint constraint : vectorEqualityConstraints)
        {
            m += constraint.apply(x).length;
        }
        for (final VectorConstraint constraint : vectorInequalityConstraints)
        {
            m += constraint.apply(x).length;
        }
        final int la = Math.max(1, m);

        final int mineq = m - meq + nPlus1 + nPlus1;

        final int lenW = (3 * nPlus1 + m) * (nPlus1 + 1) +
            (nPlus1 - meq + 1) * (mineq + 2) + 2 * mineq +
            (nPlus1 + mineq) * (nPlus1 - meq) + 2 * meq +
            nPlus1 * n / 2 + 2 * m + 3 * n + 4 * nPlus1 + 1;

        final double[] w = new double[lenW];
        final int[] jw = new int[mineq];
        final double[] g = new double[nPlus1];
        final double[] c = new double[la];

        // Note that Fortran expects arrays to be laid out in column-major order.
        final double[][] a = new double[nPlus1][la];

        double fx = 0;
        int majIterPrev = 0;
        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                // calculate the value of the objective function
                fx = wrappedObjectiveFunction.apply(x);
                final int inequalityStartIndex = copyVectorConstraints(x, c, vectorEqualityConstraints, 0);
                copyVectorConstraints(x, c, vectorInequalityConstraints, inequalityStartIndex);
            }

            if (mode[0] == 0 || mode[0] == -1)
            {
                final double[] fprime = wrappedObjectiveFunction.getJacobian(x);
                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                final int inequalityStartIndex = copyVectorConstraintJacobians(x, a, vectorEqualityConstraints, 0);
                copyVectorConstraintJacobians(x, a, vectorInequalityConstraints, inequalityStartIndex);
            }

            NativeUtils.slsqp(m, meq, la, x, xl, xu, fx, c, g, a, acc, majIter, mode, w, jw,
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
        return new OptimizeResult(x, majIter[0], mode[0], mode[0] == 0);
    }

    private OptimizeResult minimizeWithScalarConstraints(double[] x)
    {
        final int n = x.length;

        final int nPlus1 = n + 1;

        final double[] xl = new double[n]; // lower bounds
        final double[] xu = new double[n]; // upper bounds

        computeBounds(n, bounds, xl, xu);

        clip(x, xl, xu);

        final int m = scalarEqualityConstraints.size() + scalarInequalityConstraints.size();
        final int meq = scalarEqualityConstraints.size();

        final int mineq = m - meq + nPlus1 + nPlus1;
        final int lenW = (3 * nPlus1 + m) * (nPlus1 + 1) +
            (nPlus1 - meq + 1) * (mineq + 2) + 2 * mineq +
            (nPlus1 + mineq) * (nPlus1 - meq) + 2 * meq +
            nPlus1 * n / 2 + 2 * m + 3 * n + 4 * nPlus1 + 1;

        final double[] w = new double[lenW];
        final int[] jw = new int[mineq];

        final int la = Math.max(1, m);
        final double[] c = new double[la];
        final double[] g = new double[nPlus1];

        // Note that Fortran expects arrays to be laid out in column-major order.
        final double[][] a = new double[nPlus1][la];

        double fx = 0;
        int majIterPrev = 0;
        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                // calculate the value of the objective function
                fx = wrappedObjectiveFunction.apply(x);
                final int inequalityStartIndex = copyScalarConstraints(x, c, scalarEqualityConstraints, 0);
                copyScalarConstraints(x, c, scalarInequalityConstraints, inequalityStartIndex);
            }

            if (mode[0] == 0 || mode[0] == -1)
            {
                final double[] fprime = wrappedObjectiveFunction.getJacobian(x);
                System.arraycopy(fprime, 0, g, 0, n);
                g[n] = 0;
                final int inequalityStartIndex = copyScalarConstraintJacobians(x, a, scalarEqualityConstraints, 0);
                copyScalarConstraintJacobians(x, a, scalarInequalityConstraints, inequalityStartIndex);
            }

            NativeUtils.slsqp(m, meq, la, x, xl, xu, fx, c, g, a, acc, majIter, mode, w, jw,
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
        return new OptimizeResult(x, majIter[0], mode[0], mode[0] == 0);
    }

    private int copyVectorConstraintJacobians(double[] x, double[][] a, Set<VectorConstraint> constraints, int index)
    {
        for (final VectorConstraint constraint : constraints)
        {
            index += copyVectorConstraintJacobian(constraint, x, a, index);
        }
        return index;
    }

    private int copyVectorConstraintJacobian(VectorConstraint constraint, double[] x, double[][] a, int i)
    {
        final double[][] constraintJac = constraint.getJacobian(x);
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
        return constraintJac[0].length;
    }

    private int copyVectorConstraints(double[] x, double[] c, Set<VectorConstraint> constraints, int index)
    {
        for (final VectorConstraint constraint : constraints)
        {
            final double[] constraintVec = constraint.apply(x);
            System.arraycopy(constraintVec, 0, c, index, constraintVec.length);
            index += constraintVec.length;
        }
        return index;
    }

    private int copyScalarConstraintJacobians(double[] x, double[][] a, Set<ScalarConstraint> constraints, int index)
    {
        for (final ScalarConstraint constraint : constraints)
        {
            final double[] constraintJac = constraint.getJacobian(x);
            for (int j = 0; j < constraintJac.length; j++)
            {
                a[j][index] = constraintJac[j];
            }
            index++;
        }
        return index;
    }

    private int copyScalarConstraints(double[] x, double[] c, Set<ScalarConstraint> constraints, int index)
    {
        for (final ScalarConstraint constraint : constraints)
        {
            final double constraintVal = constraint.apply(x);
            c[index++] = constraintVal;
        }
        return index;
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
