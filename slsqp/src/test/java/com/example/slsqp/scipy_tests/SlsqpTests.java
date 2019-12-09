package com.example.slsqp.scipy_tests;

import com.example.slsqp.NoOpSlsqp;
import com.example.slsqp.OptimizeResult;
import com.example.slsqp.Slsqp;
import com.example.slsqp.functions.Vector2ScalarFunc;
import com.example.slsqp.constraints.ConstraintType;
import com.example.slsqp.constraints.ScalarConstraint;
import com.example.slsqp.constraints.VectorConstraint;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class SlsqpTests
{
    final double defaultTol = 1.0E-6;
    final int defaultMaxIter = 100;

    // this test case and the one below it are taken from the example at
    // https://stackoverflow.com/questions/26882087/python-scipy-optimization-minimize-using-slsqp-showing-maximized-results
    /*@Test
    public void testSymmetricInput()
    {
        double[] xl = new double[]{0, 0};
        double[] xu = new double[]{100, 5};

        final double[][] bounds = new double[][] {xl, xu};

        double[] x = new double[]{0, 0};
        final ScalarConstraint scalarConstraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.TestConstraintFunc(),
            null
        );
        final List<ScalarConstraint> constraintList = new ArrayList<>();
        constraintList.add(scalarConstraint);
        final Vector2ScalarFunc inputFunc = new TestUtil.TestInputFunc();
        final Slsqp slsqp = new Slsqp(
            inputFunc,
            null,
            x
        );
        double tolerance = 1.0E-6;
        int maxIter = 100;
        final OptimizeResult result = slsqp.minimize_slsqp_with_scalar_constraints(
            bounds,
            0,
            constraintList,
            tolerance,
            maxIter,
            null
        );

        assertTrue(Math.abs(result.x[0] - 2.5) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - 2.5) < TestUtil.ERROR);
    }

    @Test
    public void testASymmetricInput()
    {
        double[] xl = new double[]{0, 0};
        double[] xu = new double[]{100, 5};

        final double[][] bounds = new double[][] {xl, xu};

        double[] x = new double[]{1, 0};
        final ScalarConstraint scalarConstraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.TestConstraintFunc(),
            null
        );
        final List<ScalarConstraint> constraintList = new ArrayList<>();
        constraintList.add(scalarConstraint);
        final Vector2ScalarFunc inputFunc = new TestUtil.TestInputFunc();
        final Slsqp slsqp = new Slsqp(
            inputFunc,
            null,
            x
        );
        double tolerance = 1.0E-6;
        int maxIter = 100;
        final OptimizeResult result = slsqp.minimize_slsqp_with_scalar_constraints(
            bounds,
            0,
            constraintList,
            tolerance,
            maxIter,
            null
        );

        assertTrue(Math.abs(result.x[0] - 5) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - 3.18634e-14) < TestUtil.ERROR);
    }*/

    @Test
    public void test_minimize_equality_approximated()
    {
        final Slsqp slsqp = new Slsqp(new TestUtil.Fun(), null, new double[]{-1, -1});
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null);
        constraints.add(constraint);
        final OptimizeResult result = slsqp.minimize_slsqp_with_vector_constraints(
            null,
            -1,
            constraints,
            defaultTol,
            defaultMaxIter,
            null
        );
        final double[] resX = result.x;
        final double[] expected = {1, 1};
        assertArrayEquals(resX, expected);
        assertTrue(result.success);
    }

    @Test
    public void test_no_op_slsqp()
    {
        final double[] x = new double[]{-1, -1};
        final NoOpSlsqp slsqp = new NoOpSlsqp(new TestUtil.Fun(), null, x);
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint eqCon = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null);
        final VectorConstraint ieqCon = new VectorConstraint(ConstraintType.INEQ, new TestUtil.Fieqcon(), null);
        constraints.add(eqCon);
        constraints.add(ieqCon);
        final OptimizeResult result = slsqp.minimize_slsqp_with_vector_constraints(
            null,
            -1,
            constraints,
            defaultTol,
            defaultMaxIter,
            null
        );
        final double[][] a = slsqp.getA();
        final double[] resX = result.x;
        assertArrayEquals(resX, x);
        for (int i = 0; i < a.length; i++)
        {
            for (int j = 0; j < a[0].length; j++)
            {
                assertTrue(Math.abs(a[i][j] - result.a[i][j]) < TestUtil.ERROR);
            }
        }
    }
}
