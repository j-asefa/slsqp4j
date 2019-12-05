package com.example.slsqp;

public class InputFunc implements Func
{
    @Override
    public double func(double[] x)
    {
        return x[0] * x[1];
    }
}
