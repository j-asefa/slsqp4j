package com.example.stl;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

public class Stl {

    private Stl()
    {

    }

    public static native void square_cube(int i, int[] isquare, int[] icube);

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
        try (InputStream is = Stl.class.getResourceAsStream(lib))
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
        loadLib("/libsquare_cube.so");
    }
}