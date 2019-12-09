package com.example.slsqp;

public class Jacobian
{
    public static double epsilon = Math.sqrt(Math.ulp((double)1)); // 1.4901161193847656e-08

    public static double[][] approx_jacobian(double[] x, Vector2VectorFunc func)
    {
        final int n = x.length;
        double[] f0 = func.apply(x);
        double[][] jac = new double[n][f0.length];
        double[] dx = new double[n];

        for (int i = 0; i < n; i++)
        {
            dx[i] = Jacobian.epsilon;

            double[] add = new double[n];
            for (int j = 0; j < n; j++)
            {
                add[j] = x[j] + dx[j];
            }
            for (int j = 0; j < f0.length; j++)
            {
                jac[i][j] = (func.apply(add)[j] - f0[j]) / Jacobian.epsilon;
            }
            dx[i] = 0;
        }
        return jac;
    }

    // jacobian of a scalar valued function is n-dimensional vector
    public static double[] approx_jacobian(double[] x, Vector2ScalarFunc func)
    {
        final int n = x.length;
        double f0 = func.apply(x);

        double[] jac = new double[n];
        double[] dx = new double[n];
        for (int i = 0; i < n; i++)
        {
            dx[i] = Jacobian.epsilon;

            double[] add = new double[n];
            for (int j = 0; j < n; j++)
            {
                add[j] = x[j] + dx[j];
            }
            jac[i] = (func.apply(add) - f0) / Jacobian.epsilon;
            dx[i] = 0;
        }
        return jac;
    }

    public static double[] approx_jacobian(double[] x, Vector2ScalarFunc func, double arg)
    {
        final int n = x.length;
        double f0 = func.apply(x);

        double[] jac = new double[n];
        double[] dx = new double[n];
        for (int i = 0; i < n; i++)
        {
            dx[i] = Jacobian.epsilon;

            double[] add = new double[n];
            for (int j = 0; j < n; j++)
            {
                add[j] = x[j] + arg + dx[j];
            }
            jac[i] = (func.apply(add) - f0) / Jacobian.epsilon;
            dx[i] = 0;
        }
        return jac;
    }

    public static double[][] transpose(double[][] arr)
    {
        double[][] temp = new double[arr[0].length][arr.length];
        for (int i = 0; i < arr.length; i++)
            for (int j = 0; j < arr[0].length; j++)
                temp[j][i] = arr[i][j];
        return temp;
    }

}
