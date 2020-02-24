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

import com.skew.slsqp4j.functions.Vector2ScalarFunc;
import com.skew.slsqp4j.functions.Vector2VectorFunc;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class JacobianTest
{
    @Test
    public void testJacobian()
    {
        final TestFunc testFunc = new TestFunc();

        final double[] inputArr = new double[]{-1, -1};
        final double[] fprime = Jacobian.approxJacobian(inputArr, testFunc);

        final double[] exactJac = testFunc.exactJac(inputArr);
        for (int i = 0; i < inputArr.length; i++)
        {
            assertTrue(Math.abs(fprime[i] - exactJac[i]) < TestUtil.ERROR);
        }
    }

    @Test
    public void testIdentity()
    {
        final Vector2VectorFunc testFunc = (x, arg) -> x;

        final double[] inputArr = new double[]{-1, -1, 2, 4, 5, 6, 8};
        final double[][] fprime = Jacobian.approxJacobian(inputArr, testFunc);

        for (int i = 0; i < fprime.length; i++)
        {
            for (int j = 0; j < fprime[0].length; j++)
            {
                if (i == j)
                {
                    assertTrue(Math.abs(fprime[i][j] - 1) < TestUtil.ERROR);
                }
                else
                {
                    assertTrue(Math.abs(fprime[i][j] - 0) < TestUtil.ERROR);
                }
            }
        }
    }

    private static final class TestFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            final double x1 = x[0];
            final double x2 = x[1];
            return 2 * x1 * x2 + 2 * x1 - Math.pow(x1, 2) - 2 * Math.pow(x2, 2);
        }

        // returns the analytical jacobian of func
        double[] exactJac(double[] x)
        {
            final double x1 = x[0];
            final double x2 = x[1];
            return new double[] {-2 * x1 + 2 * x2 + 2, 2 * x1 - 4 * x2};
        }
    }
}
