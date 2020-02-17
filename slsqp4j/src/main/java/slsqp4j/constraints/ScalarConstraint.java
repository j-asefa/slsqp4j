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
package slsqp4j.constraints;

import slsqp4j.Jacobian;
import slsqp4j.functions.Vector2ScalarFunc;
import slsqp4j.functions.Vector2VectorFunc;

/**
 * A scalar-valued constraint function.
 */
public final class ScalarConstraint
{
    private ConstraintType constraintType;
    private Vector2ScalarFunc constraintFunc;
    private Vector2VectorFunc jacobian;
    private double[] args;

    public static class ScalarConstraintBuilder
    {
        private ConstraintType constraintType;
        private Vector2ScalarFunc constraintFunc;
        private Vector2VectorFunc jacobian;
        private double[] args;

        /**
         * Set the constraint function and any arguments it takes.
         *
         * @param constraintFunc constraint function.
         * @param args arguments, if any, to the constraint function.
         * @return this builder.
         */
        public ScalarConstraintBuilder withConstraintFunction(Vector2ScalarFunc constraintFunc, double... args)
        {
            this.constraintFunc = constraintFunc;
            this.args = args;
            return this;
        }

        /**
         * Set the type of this constraint - equality or inequality. See {@link ConstraintType}
         *
         * @param constraintType type of this constraint.
         * @return this builder.
         */
        public ScalarConstraintBuilder withConstraintType(ConstraintType constraintType)
        {
            this.constraintType = constraintType;
            return this;
        }

        /**
         * Set the analytical Jacobian of this constraint function, if any. If no analytical Jacobian is specified,
         * numerical approximation is performed on calls to {@link ScalarConstraint#getJacobian(double[])}.
         *
         * @param jacobian analytical jacobian of the constraint function specified in
         * {@link #withConstraintFunction(Vector2ScalarFunc, double...)}
         * @return this builder.
         */
        public ScalarConstraintBuilder withJacobian(Vector2VectorFunc jacobian)
        {
            this.jacobian = jacobian;
            return this;
        }

        /**
         * Build an instance of a {@link ScalarConstraint}.
         *
         * @return a new {@link ScalarConstraint} with the properties specified in this builder.
         */
        public ScalarConstraint build()
        {
            if (this.constraintType == null)
            {
                throw new IllegalStateException("must specify a constraint type");
            }
            if (this.constraintFunc == null)
            {
                throw new IllegalStateException("must specify a constraint function");
            }
            final ScalarConstraint scalarConstraint = new ScalarConstraint();
            scalarConstraint.constraintType = this.constraintType;
            scalarConstraint.constraintFunc = this.constraintFunc;
            scalarConstraint.jacobian = this.jacobian;
            scalarConstraint.args = this.args;
            return scalarConstraint;
        }
    }

    private ScalarConstraint()
    {
    }

    /**
     * Get the type of this constraint - equality or inequality. See {@link ConstraintType}
     *
     * @return the type of this constraint.
     */
    public ConstraintType getConstraintType()
    {
        return constraintType;
    }

    /**
     * Returns the Jacobian of this constraint function, evaluated at x. If an analytical Jacobian was included
     * in the builder, returns the analytical Jacobian, otherwise it returns the approximate Jacobian via
     * {@link Jacobian#approxJacobian(double[], Vector2ScalarFunc, double...)}
     *
     * @param x point at which to evaluate the jacobian of this function.
     * @return n x 1 vector consisting of the Jacobian of this constraint function, evaluated at x.
     */
    public double[] getJacobian(double[] x)
    {
        if (jacobian == null)
        {
            return Jacobian.approxJacobian(x, constraintFunc, args);
        }
        else
        {
            return jacobian.apply(x, args);
        }
    }

    /**
     * Apply this constraint function at point x.
     *
     * @param x vector-valued point at which to apply this constraint function.
     * @return the value of this constraint function applied at point x.
     */
    public double apply(double[] x)
    {
        return constraintFunc.apply(x, args);
    }
}
