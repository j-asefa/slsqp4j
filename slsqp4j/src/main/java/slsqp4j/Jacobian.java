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

import slsqp4j.functions.Vector2ScalarFunc;
import slsqp4j.functions.Vector2VectorFunc;

/**
 * Utility class for calculating the approximate jacobian of {@link Vector2ScalarFunc} functions and
 * {@link Vector2VectorFunc} functions using forward differences.
 */
public class Jacobian
{
    public static double epsilon = Math.sqrt(Math.ulp((double)1)); // 1.4901161193847656e-08

    /**
     * Approximate the Jacobian of <code>func</code> using forward differences.
     *
     * @param x the input to <code>func</code> at which to calculate the Jacobian.
     * @param func the function for which to calculate the Jacobian.
     * @param arg any arguments the function takes.
     * @return a numerical approximation of the Jacobian of <code>func</code>, using forward differences.
     */
    public static double[][] approxJacobian(double[] x, Vector2VectorFunc func, double... arg)
    {
        final int n = x.length;
        final double[] f0;
        if (arg != null && arg.length > 0)
        {
            f0 = func.apply(x, arg);
        }
        else
        {
            f0 = func.apply(x);
        }

        final double[][] jac = new double[n][f0.length];
        final double[] dx = new double[n];

        for (int i = 0; i < n; i++)
        {
            dx[i] = Jacobian.epsilon;

            final double[] add = new double[n];
            for (int j = 0; j < n; j++)
            {
                add[j] = x[j] + dx[j];
            }
            for (int j = 0; j < f0.length; j++)
            {
                if (arg != null && arg.length > 0)
                {
                    jac[i][j] = (func.apply(add, arg)[j] - f0[j]) / Jacobian.epsilon;
                }
                else
                {
                    jac[i][j] = (func.apply(add)[j] - f0[j]) / Jacobian.epsilon;
                }
            }
            dx[i] = 0;
        }
        return jac;
    }

    /**
     * Approximate the Jacobian of <code>func</code> using forward differences.
     *
     * @param x the input to <code>func</code> at which to calculate the Jacobian.
     * @param func the function for which to calculate the Jacobian.
     * @param arg any arguments the function takes.
     * @return a numerical approximation of the Jacobian of <code>func</code>, using forward differences.
     */
    public static double[] approxJacobian(double[] x, Vector2ScalarFunc func, double... arg)
    {
        final int n = x.length;
        final double f0;
        if (arg != null && arg.length > 0)
        {
            f0 = func.apply(x, arg);
        }
        else
        {
            f0 = func.apply(x);
        }

        final double[] jac = new double[n];
        final double[] dx = new double[n];
        for (int i = 0; i < n; i++)
        {
            dx[i] = Jacobian.epsilon;

            final double[] add = new double[n];
            for (int j = 0; j < n; j++)
            {
                add[j] = x[j] + dx[j];
            }

            if (arg != null && arg.length > 0)
            {
                jac[i] = (func.apply(add, arg) - f0) / Jacobian.epsilon;
            }
            else
            {
                jac[i] = (func.apply(add) - f0) / Jacobian.epsilon;
            }

            dx[i] = 0;
        }
        return jac;
    }

    /**
     * Transpose the given matrix.
     *
     * @param mat the matrix to transpose.
     * @return the transpose of mat.
     */
    public static double[][] transpose(double[][] mat)
    {
        final double[][] temp = new double[mat[0].length][mat.length];
        for (int i = 0; i < mat.length; i++)
        {
            for (int j = 0; j < mat[0].length; j++)
            {
                temp[j][i] = mat[i][j];
            }
        }
        return temp;
    }
}
