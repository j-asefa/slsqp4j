package com.example.slsqp;

public class Main
{
    public static void main(String[] args)
    {
        int m = 1; // standard int
        int meq = 1; // standard int
        int la = 2; // la = len(c). check len(c) >= la
        int n = 1; // n = len(x). check len(x) >= n
        double[] x = new double[]{1}; // array of length n    ---  value is returned to caller
        double[] xl = new double[]{1}; // array of length n
        double[] xu = new double[]{1}; // array of length n
        double f = 1; // standard double
        double[] c = new double[]{1, 1}; // array of length la
        double[] g = new double[]{1,1}; // array of length n + 1

        double[][] a = new double[la][n + 1]; // matrix of dims (la = ; n + 1)
        for (int i = 0; i < la; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                if (i == j)
                {
                    a[i][j] = 1;
                }
                else
                {
                    a[i][j] = 0;
                }
            }
        }

        double[] acc = new double[]{1}; // standard double --- value is returned to caller
        int[] iter = new int[]{1}; // standard int -- value is returned to caller
        int[] mode = new int[]{0}; // standard int -- value is returned to caller
        double[] w = new double[]{1}; // array of length l_w
        int l_w = 1; // standard int. check len(w) >= l_w
        int[] jw = new int[]{1}; // array of length l_jw
        int l_jw = 1; // standard int. check len(jw) >= l_jw
        double[] alpha = new double[]{1}; // standard double  -- value is returned to caller
        double[] f0 = new double[]{1}; // standard double -- value is returned to caller
        double[] gs = new double[]{1}; // standard double -- value is returned to caller
        double[] h1 = new double[]{1}; // standard double -- value is returned to caller
        double[] h2 = new double[]{1}; // standard double -- value is returned to caller
        double[] h3 = new double[]{1}; // standard double -- value is returned to caller
        double[] h4 = new double[]{1}; // standard double -- value is returned to caller
        double[] t = new double[]{1}; // standard double -- value is returned to caller
        double[] t0 = new double[]{1}; // standard double -- value is returned to caller
        double[] tol = new double[]{1}; // standard double -- value is returned to caller
        int[] iexact = new int[]{1}; // standard int -- value is returned to caller
        int[] incons = new int[]{1}; // standard int -- value is returned to caller
        int[] ireset = new int[]{1}; // standard int -- value is returned to caller
        int[] itermx = new int[]{1}; // standard int -- value is returned to caller
        int[] line = new int[]{1}; // standard int -- value is returned to caller
        int[] n1 = new int[]{1}; // standard int -- value is returned to caller
        int[] n2 = new int[]{1}; // standard int -- value is returned to caller
        int[] n3 = new int[]{1};// standard int -- value is returned to caller

        final int result = Slsqp.slsqp(
            m,
            meq,
            la,
            n,
            x,
            xl,
            xu,
            f,
            c,
            g,
            a,
            acc,
            iter,
            mode,
            w,
            l_w,
            jw,
            l_jw,
            alpha,
            f0,
            gs,
            h1,
            h2,
            h3,
            h4,
            t,
            t0,
            tol,
            iexact,
            incons,
            ireset,
            itermx,
            line,
            n1,
            n2,
            n3
        );
        System.out.println("result = " + result);
    }
}
