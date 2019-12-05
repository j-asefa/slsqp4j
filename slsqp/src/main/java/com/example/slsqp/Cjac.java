package com.example.slsqp;

public class Cjac
{
    private double epsilon = Math.sqrt(Math.ulp(1));
    private Func func;

    public Cjac(Func func)
    {
        this.func = func;
    }

    public Func getFunc()
    {
        return func;
    }

    public double[][] jacobian(double[] x)
    {
        return approx_jacobian(x, func);
    }

    private double[][] approx_jacobian(double[] x, Func func)
    {
        double[] f0 = new double[]{func.func(x)};

        double[][] jac = new double[x.length][f0.length];
        double[] dx = new double[x.length];
        for (int i = 0; i < x.length; i++)
        {
            dx[i] = epsilon;

            double[] add = new double[x.length];
            for (int j = 0; j < x.length; j++)
            {
                add[j] = x[j] + dx[j];
            }
            jac[i] = new double[] {(func.func(add) - f0[0]) / epsilon};
            dx[i] = 0;
        }
        return transpose(jac);
    }

    private double[][] transpose(double[][] arr)
    {
        double[][] temp = new double[arr[0].length][arr.length];
        for (int i = 0; i < arr.length; i++)
            for (int j = 0; j < arr[0].length; j++)
                temp[j][i] = arr[i][j];
        return temp;
    }

}
