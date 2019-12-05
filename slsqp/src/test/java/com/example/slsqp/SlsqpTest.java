package com.example.slsqp;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class SlsqpTest
{
    int meq = 1; // number of equality constraints
    int la = 1; // la = len(c). check len(c) >= la
    int l_w; // length of w array
    int l_jw; // length of jw array

    int[] iter = new int[]{100}; // standard int -- value is returned to caller
    int[] mode = new int[]{0}; // standard int -- value is returned to caller

    double[] alpha = new double[]{0}; // standard double  -- value is returned to caller
    double[] f0 = new double[]{0}; // standard double -- value is returned to caller
    double[] gs = new double[]{0}; // standard double -- value is returned to caller
    double[] h1 = new double[]{0}; // standard double -- value is returned to caller
    double[] h2 = new double[]{0}; // standard double -- value is returned to caller
    double[] h3 = new double[]{0}; // standard double -- value is returned to caller
    double[] h4 = new double[]{0}; // standard double -- value is returned to caller
    double[] t = new double[]{0}; // standard double -- value is returned to caller
    double[] t0 = new double[]{0}; // standard double -- value is returned to caller
    double[] tol = new double[]{0}; // standard double -- value is returned to caller
    int[] iexact = new int[]{0}; // standard int -- value is returned to caller
    int[] incons = new int[]{0}; // standard int -- value is returned to caller
    int[] ireset = new int[]{0}; // standard int -- value is returned to caller
    int[] itermx = new int[]{0}; // standard int -- value is returned to caller
    int[] line = new int[]{0}; // standard int -- value is returned to caller
    int[] n1 = new int[]{0}; // standard int -- value is returned to caller
    int[] n2 = new int[]{0}; // standard int -- value is returned to caller
    int[] n3 = new int[]{0};// standard int -- value is returned to caller


    @BeforeEach
    public void before()
    {

    }

    @Test
    public void testSlsqp()
    {
        double[] xl = new double[]{0, 100};
        double[] xu = new double[]{0, 5};

        double[] x = new double[]{0, 0};
        int n = x.length;

        int m = meq;
        int n_1 = n + 1;
        int mineq = m - meq+ n_1 + n_1;
        l_jw = mineq;
        l_w = (3*n_1+m)*(n_1+1)+(n_1-meq+1)*(mineq+2) + 2*mineq+(n_1+mineq)*(n_1-meq) + 2*meq + n_1 + ((n+1)*n);

        double[] w = new double[l_w]; // array of length l_w
        int[] jw = new int[l_jw]; // array of length l_jw

        double[][] a = new double[la][n + 1]; // matrix of dims (la = ; n + 1)
        for (int i = 0; i < la; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                if (i == j)
                {
                    a[i][j] = 0.5;
                }
                else
                {
                    a[i][j] = 0.7;
                }
            }
        }

        for (int i = 0; i < l_w; i++)
        {
            w[i] = i;
        }

        for (int i = 0; i < l_jw; i++)
        {
            jw[i] = i;
        }

        final Func constraintFunc = new ConstraintFunc();
        final Func inputFunc = new InputFunc();
        double ftol = 1.0E-6;
        double[] acc = new double[]{ftol};

        final Cjac constraintCjac = new Cjac(constraintFunc);
        final Cjac inputCjac = new Cjac(inputFunc);

        double[] c = new double[la];
        double[] g = new double[n + 1];

        double fx = inputFunc.func(x);
        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                fx = inputFunc.func(x);
                double ceq = constraintFunc.func(x);
                c = new double[] {ceq};
            }
            if (mode[0] == 0 || mode[0] == -1)
            {
                double[][] fprime = inputCjac.jacobian(x);

                for (int i = 0; i < n; i++)
                {
                    g[i] = fprime[0][i];
                }
                g[n] = 0;


                double ceq = constraintFunc.func(x);
                c = new double[] {ceq};

                a = constraintCjac.jacobian(x);
            }

            System.out.println(" ************* BEFORE ******************");
            System.out.println("m = " + m);
            System.out.println("meq = " + meq);
            System.out.println("la = " + la);
            System.out.println("n = " + n);
            System.out.println("x = " + x[0]);
            System.out.println("xl = " + xl[0]);
            System.out.println("xu = " + xu[0]);
            System.out.println("f = " + fx);
            System.out.println("c = " + c[0]);
            System.out.println("g = " + g[0] + ", " + g[1]);
            for (int i = 0; i < la; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    System.out.print("a[" + i + "]" + "[" + j + "] = " + a[i][j] + " ");
                }
                System.out.println();
            }
            System.out.println("iter = " + iter[0]);
            System.out.println("mode = " + mode[0]);
            System.out.println("w = " + w[0]);
            System.out.println("l_w = " + l_w);
            System.out.println("jw = " + jw[0]);
            System.out.println("l_jw = " + l_jw);
            System.out.println("alpha = " + alpha[0]);
            System.out.println("f0 = " + f0[0]);
            System.out.println("gs = " + gs[0]);
            System.out.println("h1 = " + h1[0]);
            System.out.println("h2 = " + h2[0]);
            System.out.println("h3 = " + h3[0]);
            System.out.println("h4 = " + h4[0]);
            System.out.println("t = " + t[0]);
            System.out.println("t0 = " + t0[0]);
            System.out.println("tol = " + tol[0]);
            System.out.println("iexact = " + iexact[0]);
            System.out.println("incons = " + incons[0]);
            System.out.println("ireset = " + ireset[0]);
            System.out.println("itermx = " + itermx[0]);
            System.out.println("line = " + line[0]);
            System.out.println("n1 = " + n1[0]);
            System.out.println("n2 = " + n2[0]);
            System.out.println("n3 = " + n3[0]);

            final int result = Slsqp.slsqp(
                m,
                meq,
                la,
                n,
                x,
                xl,
                xu,
                fx,
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

            System.out.println(" ************* AFTER ******************");
            System.out.println("m = " + m);
            System.out.println("meq = " + meq);
            System.out.println("la = " + la);
            System.out.println("n = " + n);
            System.out.println("x = " + x[0]);
            System.out.println("xl = " + xl[0]);
            System.out.println("xu = " + xu[0]);
            System.out.println("f = " + fx);
            System.out.println("c = " + c[0]);
            System.out.println("g = " + g[0] + ", " + g[1]);
            for (int i = 0; i < la; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    System.out.print("a[" + i + "]" + "[" + j + "] = " + a[i][j] + " ");
                }
                System.out.println();
            }
            System.out.println("iter = " + iter[0]);
            System.out.println("mode = " + mode[0]);
            System.out.println("w = " + w[0]);
            System.out.println("l_w = " + l_w);
            System.out.println("jw = " + jw[0]);
            System.out.println("l_jw = " + l_jw);
            System.out.println("alpha = " + alpha[0]);
            System.out.println("f0 = " + f0[0]);
            System.out.println("gs = " + gs[0]);
            System.out.println("h1 = " + h1[0]);
            System.out.println("h2 = " + h2[0]);
            System.out.println("h3 = " + h3[0]);
            System.out.println("h4 = " + h4[0]);
            System.out.println("t = " + t[0]);
            System.out.println("t0 = " + t0[0]);
            System.out.println("tol = " + tol[0]);
            System.out.println("iexact = " + iexact[0]);
            System.out.println("incons = " + incons[0]);
            System.out.println("ireset = " + ireset[0]);
            System.out.println("itermx = " + itermx[0]);
            System.out.println("line = " + line[0]);
            System.out.println("n1 = " + n1[0]);
            System.out.println("n2 = " + n2[0]);
            System.out.println("n3 = " + n3[0]);

            if (Math.abs(mode[0]) != 1)
            {
                break;
            }
        }
    }

    private void minimize(
        Func inputFuncInter,
        double[] x0,
        double[][] bounds,
        Constraint constraint,
        int maxiter,
        double ftol,
        int iprint)
    {

    }

    private Cjac cjacFactory(Func func)
    {
        return new Cjac(func);
    }

    private double lambda(double[] x)
    {
        return x[0] + x[1] - 5;
    }

    /*@Test
    public void testSlsqp()
    {
        int m = 1; // standard int
        int meq = 1; // standard int
        int la = 2; // la = len(c). check len(c) >= la
        int n = 1; // n = len(x). check len(x) >= n

        int n_1 = n + 1;
        int l_jw = m - meq + n_1 + n_1;
        int l_w = (3*n_1+m)*(n_1+1)+(n_1-meq+1)*(l_jw+2) + 2*l_jw+(n_1+l_jw)*(n_1-meq)
            + 2*meq + n_1 + ((n+1)*n);

        double[] x = new double[]{0.953}; // array of length n    ---  value is returned to caller
        double[] xl = new double[]{1}; // array of length n
        double[] xu = new double[]{1}; // array of length n
        double f = 1; // standard double
        double[] c = new double[]{0.5, 0.5}; // array of length la
        double[] g = new double[]{1,1}; // array of length n + 1

        double[][] a = new double[la][n + 1]; // matrix of dims (la = ; n + 1)
        for (int i = 0; i < la; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                if (i == j)
                {
                    a[i][j] = 0.5;
                }
                else
                {
                    a[i][j] = 0.7;
                }
            }
        }

        double[] acc = new double[]{-0.5}; // standard double --- value is returned to caller
        int[] iter = new int[]{50}; // standard int -- value is returned to caller
        int[] mode = new int[]{0}; // standard int -- value is returned to caller
        double[] w = new double[l_w]; // array of length l_w
        for (int i = 0; i < l_w; i++)
        {
            w[i] = i;
        }
        int[] jw = new int[l_jw]; // array of length l_jw
        for (int i = 0; i < l_jw; i++)
        {
            jw[i] = i;
        }
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

        System.out.println(" ************* BEFORE ******************");
        System.out.println("m = " + m);
        System.out.println("meq = " + meq);
        System.out.println("la = " + la);
        System.out.println("n = " + n);
        System.out.println("x = " + x[0]);
        System.out.println("xl = " + xl[0]);
        System.out.println("xu = " + xu[0]);
        System.out.println("f = " + f);
        System.out.println("c = " + c[0] + ", " + c[1]);
        System.out.println("g = " + g[0] + ", " + g[1]);
        for (int i = 0; i < la; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                System.out.print("a[" + i + "]" + "[" + j + "] = " + a[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println("iter = " + iter[0]);
        System.out.println("mode = " + mode[0]);
        System.out.println("w = " + w[0]);
        System.out.println("l_w = " + l_w);
        System.out.println("jw = " + jw[0]);
        System.out.println("l_jw = " + l_jw);
        System.out.println("alpha = " + alpha[0]);
        System.out.println("f0 = " + f0[0]);
        System.out.println("gs = " + gs[0]);
        System.out.println("h1 = " + h1[0]);
        System.out.println("h2 = " + h2[0]);
        System.out.println("h3 = " + h3[0]);
        System.out.println("h4 = " + h4[0]);
        System.out.println("t = " + t[0]);
        System.out.println("t0 = " + t0[0]);
        System.out.println("tol = " + tol[0]);
        System.out.println("iexact = " + iexact[0]);
        System.out.println("incons = " + incons[0]);
        System.out.println("ireset = " + ireset[0]);
        System.out.println("itermx = " + itermx[0]);
        System.out.println("line = " + line[0]);
        System.out.println("n1 = " + n1[0]);
        System.out.println("n2 = " + n2[0]);
        System.out.println("n3 = " + n3[0]);

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

        System.out.println(" ************* AFTER *******************");
        System.out.println("result = " + result);
        System.out.println("m = " + m);
        System.out.println("meq = " + meq);
        System.out.println("la = " + la);
        System.out.println("n = " + n);
        System.out.println("x = " + x[0]);
        System.out.println("xl = " + xl[0]);
        System.out.println("xu = " + xu[0]);
        System.out.println("f = " + f);
        System.out.println("c = " + c[0] + ", " + c[1]);
        System.out.println("g = " + g[0] + ", " + g[1]);
        for (int i = 0; i < la; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                System.out.print("a[" + i + "]" + "[" +j + "] = " + a[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println("iter = " + iter[0]);
        System.out.println("mode = " + mode[0]);
        System.out.println("w = " + w[0]);
        System.out.println("l_w = " + l_w);
        System.out.println("jw = " + jw[0]);
        System.out.println("l_jw = " + l_jw);
        System.out.println("alpha = " + alpha[0]);
        System.out.println("f0 = " + f0[0]);
        System.out.println("gs = " + gs[0]);
        System.out.println("h1 = " + h1[0]);
        System.out.println("h2 = " + h2[0]);
        System.out.println("h3 = " + h3[0]);
        System.out.println("h4 = " + h4[0]);
        System.out.println("t = " + t[0]);
        System.out.println("t0 = " + t0[0]);
        System.out.println("tol = " + tol[0]);
        System.out.println("iexact = " + iexact[0]);
        System.out.println("incons = " + incons[0]);
        System.out.println("ireset = " + ireset[0]);
        System.out.println("itermx = " + itermx[0]);
        System.out.println("line = " + line[0]);
        System.out.println("n1 = " + n1[0]);
        System.out.println("n2 = " + n2[0]);
        System.out.println("n3 = " + n3[0]);
    }*/

}
