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
package com.skew.slsqp4j;

import com.skew.slsqp4j.functions.Vector2MatrixFunc;
import com.skew.slsqp4j.functions.Vector2ScalarFunc;
import com.skew.slsqp4j.functions.Vector2VectorFunc;

public class TestUtil
{
    public static final double ERROR = 1.0E-5;

    // this test case is taken from the example at
    // https://stackoverflow.com/questions/26882087/python-scipy-optimization-minimize-using-slsqp-showing-maximized-results
    public static class TestConstraintFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return x[0] + x[1] - 5;
        }
    }

    public static class TestInputFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return x[0] * x[1];
        }
    }

    public static final class Fun implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] arr, double... arg)
        {
            if (arg.length > 1)
            {
                throw new IllegalArgumentException("optional argument must be one of {1, -1}");
            }

            final double x = arr[0];
            final double y = arr[1];
            double sign = 1;
            if (arg.length > 0)
            {
                sign = arg[0];
            }

            return sign * (2 * x * y + 2 * x - Math.pow(x, 2) - 2 * Math.pow(y, 2));
        }
    }

    public static final class Jac implements Vector2VectorFunc
    {
        // returns the analytical jacobian of Fun above.
        @Override
        public double[] apply(double[] arr, double... arg)
        {
            if (arg.length > 1)
            {
                throw new IllegalArgumentException("optional argument must be one of {1, -1}");
            }
            final double x = arr[0];
            final double y = arr[1];

            double sign = 1;
            if (arg.length > 0)
            {
                sign = arg[0];
            }
            final double dfdx = sign * (-2 * x + 2 * y + 2);
            final double dfdy = sign * (2 * x - 4 * y);
            return new double[] {dfdx, dfdy};
        }
    }

    public static final class Fecon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {x[0] - x[1]};
        }
    }

    public static final class FprimeEcon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {1, -1};
        }
    }

    public static final class FprimeEconScalar implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new FprimeEcon().apply(x, arg);
        }
    }

    public static final class FeconScalar implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return new Fecon().apply(x, arg)[0];
        }
    }

    public static final class Fieqcon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {x[0] - x[1] - 1};
        }
    }

    public static final class FprimeIeqcon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {1, -1};
        }
    }

    public static final class Fieqcon2 implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return x;
        }
    }

    public static final class FprimeIeqcon2 implements Vector2MatrixFunc
    {
        public double[][] apply(double[] x, double... arg)
        {
            return identity(x.length);
        }

        private double[][] identity(int length)
        {
            final double[][] identityMat = new double[length][length];
            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < length; j++)
                {
                    if (i == j)
                    {
                        identityMat[i][j] = 1;
                    }
                    else
                    {
                        identityMat[i][j] = 0;
                    }
                }
            }
            return identityMat;
        }
    }
}
