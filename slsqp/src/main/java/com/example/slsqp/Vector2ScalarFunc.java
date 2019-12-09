package com.example.slsqp;

@FunctionalInterface
public interface Vector2ScalarFunc
{
    double apply(double[] x, double... arg);
}
