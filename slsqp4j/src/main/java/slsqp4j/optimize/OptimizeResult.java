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
package slsqp4j.optimize;

/**
 * An {@link OptimizeResult} represents the result of applying {@link Slsqp#minimize(double[])}. This class contains
 * public variables that provide information about the current state of the Slsqp optimizer.
 */
public class OptimizeResult
{
    private final double[] resultVec;
    private final int numIters;
    private final int exitMode;
    private final boolean success;

    OptimizeResult(
        double[] resultVec,
        int numIters,
        int exitMode,
        boolean success
    )
    {
        this.resultVec = resultVec;
        this.numIters = numIters;
        this.exitMode = exitMode;
        this.success = success;
    }

    /**
     * Get the point at which the objective function is minimized, if one exists.
     * Must check that {@link #success()} returns true.
     *
     * @return minima, if it exists, of the objective function minimized by the Slsqp solver.
     */
    public double[] resultVec()
    {
        return resultVec;
    }

    /**
     * Get the exit mode of the Slsqp solver. Useful for determining the state of the solver and for debugging.
     *
     *     Exit modes are defined as follows ::
     *         -1 : Gradient evaluation required (g & a)
     *          0 : Optimization terminated successfully
     *          1 : Function evaluation required (f & c)
     *          2 : More equality constraints than independent variables
     *          3 : More than 3*n iterations in LSQ subproblem
     *          4 : Inequality constraints incompatible
     *          5 : Singular matrix E in LSQ subproblem
     *          6 : Singular matrix C in LSQ subproblem
     *          7 : Rank-deficient equality constraint subproblem HFTI
     *          8 : Positive directional derivative for linesearch
     *          9 : Iteration limit exceeded
     *
     * @return the exit mode of the Slsqp solver.
     */
    public int exitMode()
    {
        return exitMode;
    }

    /**
     * Get the number of iterations the Slsqp solver has performed.
     *
     * @return number of iterations the solver has performed.
     */
    public int numIters()
    {
        return numIters;
    }

    /**
     * Did the optimizer complete successfully? A value of true corresponds to an {@link #exitMode()} return value of 0.
     *
     * @return true if the Slsqp solver completed successfully false otherwise.
     */
    public boolean success()
    {
        return success;
    }
}
