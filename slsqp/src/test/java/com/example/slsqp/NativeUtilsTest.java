package com.example.slsqp;

import com.example.slsqp.constraints.ConstraintType;
import com.example.slsqp.constraints.ScalarConstraint;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

public class NativeUtilsTest
{
    /*@Test
    public void testSlsqp()
    {
        double[] xl = new double[]{0, 0};
        double[] xu = new double[]{100, 5};

        final double[][] bounds = new double[][] {xl, xu};

        double[] x = new double[]{0, 0};
        final ScalarConstraint scalarConstraint = new ScalarConstraint(
            ConstraintType.EQ,
            new ConstraintFunc(),
            null
        );
        final List<ScalarConstraint> constraintList = new ArrayList<>();
        constraintList.add(scalarConstraint);
        final Vector2ScalarFunc inputFunc = new InputFunc();
        final Slsqp slsqp = new Slsqp(
            inputFunc,
            null,
            x
        );
        double tolerance = 1.0E-6;
        int maxIter = 100;
        slsqp.minimize_slsqp_with_scalar_constraints(
            bounds,
            0,
            constraintList,
            tolerance,
            maxIter,
            null
        );
    }*/

    private static class ConstraintFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return x[0] + x[1] - 5;
        }
    }

    private static class InputFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return x[0] * x[1];
        }
    }
}
