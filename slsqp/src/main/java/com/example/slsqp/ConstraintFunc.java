package com.example.slsqp;

public class ConstraintFunc implements Func
{
    @Override
    public double func(double[] x)
    {
        return x[0] + x[1] - 5;
    }
}
