package com.example.square_cube;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

public class SquareCube
{

    private SquareCube()
    {

    }

    public static native void square_cube(int i, int[] isquare, int[] icube);

    private static void loadLib(String lib)
    {
        try (InputStream is = SquareCube.class.getResourceAsStream(lib))
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