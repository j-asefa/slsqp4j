package com.example.slsqp;

@FunctionalInterface
public interface Vector2ScalarFunc
{
    double func(double[] x, double... arg);
}
