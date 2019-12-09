package com.example.slsqp.functions;

@FunctionalInterface
public interface Vector2VectorFunc
{
    double[] apply(double[] x, double... arg);
}
