package com.example.stl;

public class Main
{
    public static void main(String[] args)
    {
        final int result = Stl.slsqp(
            1,
            2,
            1,
            1,
            new double[]{1},
            new double[]{1},
            new double[]{1},
            0.1,
            new double[]{1},
            new double[]{1},
            new double[][]{new double[]{1}, new double[]{2}},
            new double[]{0.2},
            new int[]{1},
            new int[]{-1},
            new double[]{0.1},
            1,
            new int[]{1,2},
            2,
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new double[]{0.1},
            new int[]{1},
            new int[]{3},
            new int[]{3},
            new int[]{2},
            new int[]{1},
            new int[]{2},
            new int[]{1},
            new int[]{0}
        );
    }
}
