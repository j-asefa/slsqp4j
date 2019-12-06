package com.example.slsqp;

public class InputFunc implements Vector2ScalarFunc
{
    @Override
    public double func(double[] x)
    {
        return x[0] * x[1];
    }
}
