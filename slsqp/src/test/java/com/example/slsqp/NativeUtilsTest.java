package com.example.slsqp;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNotEquals;

public class NativeUtilsTest
{
    private int[] iter = new int[]{100};
    private int[] mode = new int[]{0};

    private double[] alpha = new double[]{0};
    private double[] f0 = new double[]{0};
    private double[] gs = new double[]{0};
    private double[] h1 = new double[]{0};
    private double[] h2 = new double[]{0};
    private double[] h3 = new double[]{0};
    private double[] h4 = new double[]{0};
    private double[] t = new double[]{0};
    private double[] t0 = new double[]{0};
    private double[] tol = new double[]{0};
    private int[] iexact = new int[]{0};
    private int[] incons = new int[]{0};
    private int[] ireset = new int[]{0};
    private int[] itermx = new int[]{0};
    private int[] line = new int[]{0};
    private int[] n1 = new int[]{0};
    private int[] n2 = new int[]{0};
    private int[] n3 = new int[]{0};


    @BeforeEach
    public void before()
    {

    }

    @Test
    public void testSlsqp()
    {
        double[] xl = new double[]{0, 0};
        double[] xu = new double[]{100, 5};

        double[] x = new double[]{0, 0};
        int n = x.length;

        int meq = 1;
        int m = meq;
        int n_1 = n + 1;
        int mineq = m - meq + n_1 + n_1;

        int l_w = (3 * n_1 + m) * (n_1 + 1) + (n_1 - meq + 1) * (mineq + 2) + 2 * mineq + (n_1 + mineq) * (n_1 - meq) +
            2 * meq + n_1 * n / 2 + 2 * m + 3 * n + 4 * n_1 + 1;
        int l_jw = Math.max(mineq, n_1 - meq);

        double[] w = new double[l_w];
        int[] jw = new int[l_jw];

        int la = 1;
        double[][] a = new double[la][n + 1];
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

        final Vector2VectorFunc constraintFunc = new ConstraintFunc();
        final Vector2ScalarFunc inputFunc = new InputFunc();
        double ftol = 1.0E-6;
        double[] acc = new double[]{ftol};

        double[] c = new double[la];
        double[] g = new double[n + 1];

        final Slsqp slsqp = new Slsqp(
            m,
            meq,
            la,
            n,
            x,
            xl,
            xu,
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
            n3,
            constraintFunc,
            inputFunc
        );
        slsqp.solveSlsqp();
        assertNotEquals(Math.abs(mode[0]), 1);
    }

    private static class ConstraintFunc implements Vector2VectorFunc
    {
        @Override
        public double[] func(double[] x)
        {
            return new double[]{x[0] + x[1] - 5};
        }
    }

    private static class InputFunc implements Vector2ScalarFunc
    {
        @Override
        public double func(double[] x)
        {
            return x[0] * x[1];
        }
    }
}
