package com.example.slsqp;

public class WrappedVectorFunction
{
    private final Vector2VectorFunc func;
    private double arg = Double.MIN_VALUE;

    public WrappedVectorFunction(Vector2VectorFunc func)
    {
        this.func = func;
    }

    public WrappedVectorFunction(Vector2VectorFunc func, double arg)
    {
        this.func = func;
        this.arg = arg;
    }

    public double[][] approx_jacobian(double[] x)
    {
        final int n = x.length;
        double[] x0 = new double[n];

        if (arg != Double.MIN_VALUE)
        {
            for (int i = 0; i < n; i++)
            {
                x0[i] = x[i] + arg;
            }
        }

        double[] f0 = func.func(x0);

        double[][] jac = new double[n][f0.length];
        double[] dx = new double[n];
        for (int i = 0; i < n; i++)
        {
            dx[i] = Jacobian.epsilon;

            double[] add = new double[n];
            for (int j = 0; j < n; j++)
            {
                add[j] = x0[j] + dx[j];
            }
            for (int j = 0; j < f0.length; j++)
            {
                jac[i][j] = (func.func(add)[j] - f0[j]) / Jacobian.epsilon;
            }
            dx[i] = 0;
        }
        return Jacobian.transpose(jac);
    }
}
