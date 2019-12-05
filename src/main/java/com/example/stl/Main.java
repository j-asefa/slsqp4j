package com.example.stl;

public class Main
{
    public static void main(String[] args)
    {
        int i = 2;
        int[] isquare = {0};
        int[] icube = {0};
        System.out.println("i = " + i);
        Stl.square_cube(i, isquare, icube);

        System.out.println("isquare = " + isquare[0]);
        System.out.println("icube = " + icube[0]);
        /*final int result = Stl.slsqp(
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
        );*/
    }
}
