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
package slsqp4j.util;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Utility class for interacting with shared native libraries
 */
public final class NativeUtils
{

    private NativeUtils()
    {

    }

    // Native C method that we want to be able to call from Java.
    public static native int slsqp(
        int m, // standard int
        int meq, // standard int
        int la, // la = len(c). check len(c) >= la
        int n, // n = len(x). check len(x) >= n
        double[] x, // array of length n    ---  value is returned to caller
        double[] xl, // array of length n
        double[] xu, // array of length n
        double f, // standard double
        double[] c, // array of length la
        double[] g, // array of length n + 1
        double[][] a, // matrix of dims (la, n + 1)
        double[] acc, // standard double --- value is returned to caller
        int[] iter, // standard int -- value is returned to caller
        int[] mode, // standard int -- value is returned to caller
        double[] w, // array of length l_w
        int lenW, // standard int. check len(w) >= l_w
        int[] jw, // array of length l_jw
        int lenJw, // standard int. check len(jw) >= l_jw
        double[] alpha, // standard double  -- value is returned to caller
        double[] f0, // standard double -- value is returned to caller
        double[] gs, // standard double -- value is returned to caller
        double[] h1, // standard double -- value is returned to caller
        double[] h2, // standard double -- value is returned to caller
        double[] h3, // standard double -- value is returned to caller
        double[] h4, // standard double -- value is returned to caller
        double[] t, // standard double -- value is returned to caller
        double[] t0, // standard double -- value is returned to caller
        double[] tol, // standard double -- value is returned to caller
        int[] iexact, // standard int -- value is returned to caller
        int[] incons, // standard int -- value is returned to caller
        int[] ireset, // standard int -- value is returned to caller
        int[] itermx, // standard int -- value is returned to caller
        int[] line, // standard int -- value is returned to caller
        int[] n1, // standard int -- value is returned to caller
        int[] n2, // standard int -- value is returned to caller
        int[] n3 // standard int -- value is returned to caller
    );

    private static void loadLib(String lib)
    {
        try (InputStream is = NativeUtils.class.getResourceAsStream(lib))
        {
            File tempLib = null;
            final int dot = lib.indexOf('.');
            tempLib = File.createTempFile(lib.substring(0, dot), lib.substring(dot));
            try (FileOutputStream out = new FileOutputStream(tempLib))
            {
                final byte[] buf = new byte[1 << 18];
                while (true)
                {
                    final int read = is.read(buf);
                    if (read == -1)
                    {
                        break;
                    }
                    out.write(buf, 0, read);
                }
            }
            finally
            {
                tempLib.deleteOnExit();
            }
            System.load(tempLib.getAbsolutePath());
        }
        catch (IOException e)
        {
            throw new Error("Internal error: cannot find " + lib + ", broken package?");
        }
    }

    static
    {
        loadLib("/libslsqp.so");
    }

    public static void slsqp(
        int m,
        int meq,
        int la,
        double[] x,
        double[] xl,
        double[] xu,
        double fx,
        double[] c,
        double[] g,
        double[][] a,
        double[] acc,
        int[] majiter,
        int[] mode,
        double[] w,
        int[] jw,
        double[] alpha,
        double[] f0,
        double[] gs,
        double[] h1,
        double[] h2,
        double[] h3,
        double[] h4,
        double[] t,
        double[] t0,
        double[] tol,
        int[] iexact,
        int[] incons,
        int[] ireset,
        int[] itermx,
        int[] line,
        int[] n1,
        int[] n2,
        int[] n3)
    {
        slsqp(m, meq, la, x.length, x, xl, xu, fx, c, g, a, acc, majiter, mode, w, w.length, jw, jw.length,
            alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, iexact, incons, ireset, itermx, line, n1, n2, n3);
    }
}