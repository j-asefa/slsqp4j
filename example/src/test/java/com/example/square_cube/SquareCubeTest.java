package com.example.square_cube;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SquareCubeTest
{
    @Test
    public void testSquareCube()
    {
        int i = 2;
        int[] isquare = {0};
        int[] icube = {0};

        SquareCube.square_cube(i, isquare, icube);

        assertEquals(i, 2);

        assertEquals(isquare.length, 1);
        assertEquals(isquare[0], 4);

        assertEquals(icube.length, 1);
        assertEquals(icube[0], 8);
    }
}
