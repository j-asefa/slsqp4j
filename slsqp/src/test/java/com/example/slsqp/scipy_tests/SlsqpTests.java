package com.example.slsqp.scipy_tests;

import com.example.slsqp.NoOpSlsqp;
import com.example.slsqp.OptimizeResult;
import com.example.slsqp.Slsqp;
import com.example.slsqp.constraints.ConstraintType;
import com.example.slsqp.constraints.VectorConstraint;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

public class SlsqpTests
{
    final double defaultTol = 1.0E-6;
    final int defaultMaxIter = 100;

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
                assertEquals(a[i][j], result.a[i][j]);
            }
        }
    }
}
