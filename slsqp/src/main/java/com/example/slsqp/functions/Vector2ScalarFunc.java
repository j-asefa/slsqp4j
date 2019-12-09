package com.example.slsqp.functions;

@FunctionalInterface
public interface Vector2ScalarFunc
{
    double apply(double[] x, double... arg);
}
