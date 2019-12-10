package com.example.slsqp.functions;

@FunctionalInterface
public interface Vector2MatrixFunc
{
    double[][] apply(double[] x, double... arg);
}
