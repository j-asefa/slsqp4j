package com.example.slsqp.scipy_tests;

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
}
