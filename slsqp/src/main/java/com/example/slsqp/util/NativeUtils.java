package com.example.slsqp.util;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

public class NativeUtils
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
        int l_w, // standard int. check len(w) >= l_w
        int[] jw, // array of length l_jw
        int l_jw, // standard int. check len(jw) >= l_jw
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
        loadLib("/libslsqp_solver.so");
    }

    public static void slsqp(
        int m,
        int meq,
        int la,
        double[] x,
        double[] xl,
        double[] xu,
        double[] fx,
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
        slsqp(m, meq, la, x.length, x, xl, xu, fx[0], c, g, a, acc, majiter, mode, w, w.length, jw, jw.length,
            alpha, f0, gs, h1, h2, h3, h4, t, t0, tol, iexact, incons, ireset, itermx, line, n1, n2, n3);
    }
}