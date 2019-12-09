package com.example.slsqp.functions;

import com.example.slsqp.Jacobian;

public class WrappedScalarFunction
{
    private final Vector2ScalarFunc func;
    private double arg = Double.MIN_VALUE;

    public WrappedScalarFunction(Vector2ScalarFunc func)
    {
        this.func = func;
    }

    public WrappedScalarFunction(Vector2ScalarFunc func, double arg)
    {
        this.func = func;
        this.arg = arg;
    }

    public double[] approx_jacobian(double[] x)
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

        return Jacobian.approx_jacobian(x0, this.func);
    }
}
