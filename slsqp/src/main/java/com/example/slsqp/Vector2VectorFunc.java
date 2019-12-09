package com.example.slsqp;

@FunctionalInterface
public interface Vector2VectorFunc
{
    double[] apply(double[] x, double... arg);
}
