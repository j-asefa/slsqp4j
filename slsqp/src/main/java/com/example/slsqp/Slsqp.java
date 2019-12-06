package com.example.slsqp;

public class Slsqp
{
    private int m; // standard int
    private int meq; // standard int
    private int la; // la = len(c). check len(c) >= la
    private int n; // n = len(x). check len(x) >= n
    private double[] x; // array of length n    ---  value is returned to caller
    private double[] xl; // array of length n
    private double[] xu; // array of length n
    private double[] c; // array of length la
    private double[] g; // array of length n + 1
    private double[][] a; // matrix of dims (la; n + 1)
    private double[] acc; // standard double --- value is returned to caller
    private int[] iter; // standard int -- value is returned to caller
    private int[] mode; // standard int -- value is returned to caller
    private double[] w; // array of length l_w
    private int l_w; // standard int. check len(w) >= l_w
    private int[] jw; // array of length l_jw
    private int l_jw; // standard int. check len(jw) >= l_jw
    private double[] alpha; // standard double  -- value is returned to caller
    private double[] f0; // standard double -- value is returned to caller
    private double[] gs; // standard double -- value is returned to caller
    private double[] h1; // standard double -- value is returned to caller
    private double[] h2; // standard double -- value is returned to caller
    private double[] h3; // standard double -- value is returned to caller
    private double[] h4; // standard double -- value is returned to caller
    private double[] t; // standard double -- value is returned to caller
    private double[] t0; // standard double -- value is returned to caller
    private double[] tol; // standard double -- value is returned to caller
    private int[] iexact; // standard int -- value is returned to caller
    private int[] incons; // standard int -- value is returned to caller
    private int[] ireset; // standard int -- value is returned to caller
    private int[] itermx; // standard int -- value is returned to caller
    private int[] line; // standard int -- value is returned to caller
    private int[] n1; // standard int -- value is returned to caller
    private int[] n2; // standard int -- value is returned to caller
    private int[] n3; // standard int -- value is returned to caller
    private Vector2VectorFunc constraintFunc;
    private Vector2ScalarFunc inputFunc;

    private double fx;

    public Slsqp(
        int m,
        int meq,
        int la,
        int n,
        double[] x,
        double[] xl,
        double[] xu,
        double[] c,
        double[] g,
        double[][] a,
        double[] acc,
        int[] iter,
        int[] mode,
        double[] w,
        int l_w,
        int[] jw,
        int l_jw,
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
        int[] n3,
        Vector2VectorFunc constraintFunc,
        Vector2ScalarFunc inputFunc
    )
    {
        this.m = m;
        this.meq = meq;
        this.la = la;
        this.n = n;
        this.x = x;
        this.xl = xl;
        this.xu = xu;
        this.c = c;
        this.g = g;
        this.a = a;
        this.acc = acc;
        this.iter = iter;
        this.mode = mode;
        this.w = w;
        this.l_w = l_w;
        this.jw = jw;
        this.l_jw = l_jw;
        this.alpha = alpha;
        this.f0 = f0;
        this.gs = gs;
        this.h1 = h1;
        this.h2 = h2;
        this.h3 = h3;
        this.h4 = h4;
        this.t = t;
        this.t0 = t0;
        this.tol = tol;
        this.iexact = iexact;
        this.incons = incons;
        this.ireset = ireset;
        this.itermx = itermx;
        this.line = line;
        this.n1 = n1;
        this.n2 = n2;
        this.n3 = n3;
        this.constraintFunc = constraintFunc;
        this.inputFunc = inputFunc;
    }

    public void solveSlsqp()
    {
        while (true)
        {
            if (mode[0] == 0 || mode[0] == 1)
            {
                fx = inputFunc.func(x);
                c = constraintFunc.func(x);
            }
            if (mode[0] == 0 || mode[0] == -1)
            {
                double[][] fprime = Jacobian.approx_jacobian(x, inputFunc);

                System.arraycopy(fprime[0], 0, g, 0, n);
                g[n] = 0;


                c = constraintFunc.func(x);

                a = Jacobian.approx_jacobian(x, constraintFunc);
            }

            System.out.println(" ************* BEFORE ******************");
            System.out.println("m = " + m);
            System.out.println("meq = " + meq);
            System.out.println("la = " + la);
            System.out.println("n = " + n);
            System.out.println("x = " + x[0] + ", " + x[1]);
            System.out.println("xl = " + xl[0] + ", " + xl[1]);
            System.out.println("xu = " + xu[0] + ", " + xu[0]);
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

            final int result = NativeUtils.slsqp(
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
            System.out.println("x = " + x[0] + ", " + x[1]);
            System.out.println("xl = " + xl[0] + ", " + xl[1]);
            System.out.println("xu = " + xu[0] + ", " + xu[0]);
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
}
