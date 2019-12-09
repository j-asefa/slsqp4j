package com.example.slsqp.functions;

import com.example.slsqp.Jacobian;

public class WrappedScalarFunction
{
    private final Vector2ScalarFunc func;

    public WrappedScalarFunction(Vector2ScalarFunc func)
    {
        this.func = func;
    }

    public double[] approx_jacobian(double[] x)
    {
        return Jacobian.approx_jacobian(x, this.func);
    }
}
