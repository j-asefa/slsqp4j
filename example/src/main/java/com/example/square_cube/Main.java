package com.example.square_cube;

public class Main
{
    public static void main(String[] args)
    {
        int i = 2;
        int[] isquare = {0};
        int[] icube = {0};
        System.out.println("i = " + i);
        SquareCube.square_cube(i, isquare, icube);
        System.out.println("isquare = " + isquare[0]);
        System.out.println("icube = " + icube[0]);
    }
}
