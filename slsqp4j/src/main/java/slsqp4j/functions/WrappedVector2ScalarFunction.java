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
package slsqp4j.functions;

import slsqp4j.Jacobian;

/**
 * A wrapper around a {@link Vector2ScalarFunc} that provides access to the function and its
 * jacobian.
 */
public class WrappedVector2ScalarFunction
{
    private final Vector2ScalarFunc func;
    private final double[] arg;
    private final Vector2VectorFunc jacobian;

    public WrappedVector2ScalarFunction(Vector2ScalarFunc func, Vector2VectorFunc jacobian, double... arg)
    {
        this.func = func;
        this.jacobian = jacobian;
        this.arg = arg;
    }

    public double apply(double[] x)
    {
        return func.apply(x, arg);
    }

    public double[] getJacobian(double[] x)
    {
        if (jacobian == null)
        {
            return Jacobian.approxJacobian(x, this.func, arg);
        }
        else
        {
            return jacobian.apply(x, arg);
        }
    }
}
