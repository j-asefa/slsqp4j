package com.example.slsqp;

public class Jacobian
{
    public static double epsilon = Math.sqrt(Math.ulp((double)1)); // 1.4901161193847656e-08

    public static double[][] approx_jacobian(double[] x, Vector2VectorFunc func)
    {
        double[] f0 = func.func(x);

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
            for (int j = 0; j < f0.length; j++)
            {
                jac[i][j] = (func.func(add)[j] - f0[j]) / epsilon;
            }
            dx[i] = 0;
        }
        return transpose(jac);
    }

    public static double[][] approx_jacobian(double[] x, Vector2ScalarFunc func)
    {
        double f0 = func.func(x);

        double[][] jac = new double[x.length][1];
        double[] dx = new double[x.length];
        for (int i = 0; i < x.length; i++)
        {
            dx[i] = epsilon;

            double[] add = new double[x.length];
            for (int j = 0; j < x.length; j++)
            {
                add[j] = x[j] + dx[j];
            }

            jac[i][0] = (func.func(add) - f0) / epsilon;
            dx[i] = 0;
        }
        return transpose(jac);
    }

    private static double[][] transpose(double[][] arr)
    {
        double[][] temp = new double[arr[0].length][arr.length];
        for (int i = 0; i < arr.length; i++)
            for (int j = 0; j < arr[0].length; j++)
                temp[j][i] = arr[i][j];
        return temp;
    }

}
