package com.example.slsqp;

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

    // Scalar tests

    // this test case and the one below it are taken from the example at
    // https://stackoverflow.com/questions/26882087/python-scipy-optimization-minimize-using-slsqp-showing-maximized-results
    @Test
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
            constraintList,
            tolerance,
            maxIter,
            null,
            1
        );

        assertTrue(Math.abs(result.x[0] - 2.5) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - 2.5) < TestUtil.ERROR);
        assertTrue(result.success);
    }


    @Test
    public void testASymmetricInput()
    {
        double[] xl = new double[]{0, 0};
        double[] xu = new double[]{100, 5};

        final double[][] bounds = new double[][] {xl, xu};

        double[] x = new double[] {0.2, 0.9};
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
            constraintList,
            tolerance,
            maxIter,
            null,
            1
        );

        assertTrue(Math.abs(result.x[0] - 1.99840144e-14) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - 5) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void test_minimize_equality_given_cons_scalar()
    {
        final double[] x = new double[] {-1, 1};
        final Slsqp slsqp = new Slsqp(new TestUtil.Fun(), new TestUtil.Jac(), x);
        final List<ScalarConstraint> constraints = new ArrayList<>();
        final ScalarConstraint constraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.FeconScalar(),
            new TestUtil.FprimeEconScalar().apply(x));
        constraints.add(constraint);
        final OptimizeResult result = slsqp.minimize_slsqp_with_scalar_constraints(
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] expected = {1, 1};
        assertTrue(Math.abs(result.x[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void test_scalar_constraints()
    {
        final double[] x = new double[] {3};

        final Vector2ScalarFunc inputFunc = (x1, arg) -> Math.pow(x1[0], 2);

        final Vector2ScalarFunc constraintFunc = (x1, arg) -> x1[0] - 1;

        final Slsqp slsqp = new Slsqp(inputFunc, null, x);

        final List<ScalarConstraint> constraints = new ArrayList<>();
        final ScalarConstraint constraint = new ScalarConstraint(
            ConstraintType.INEQ,
            constraintFunc,
            null);
        constraints.add(constraint);
        final OptimizeResult result = slsqp.minimize_slsqp_with_scalar_constraints(
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            1
        );
        final double[] resX = result.x;
        final double[] expected = {1};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void test_infeasible_initial()
    {
        double[] x = new double[] {10};

        final Vector2ScalarFunc inputFunc = (x1, arg) -> Math.pow(x1[0], 2) - 2*x1[0] + 1;

        final Vector2ScalarFunc constraintFunc1 = (x1, arg) -> 0 - x1[0];
        final ScalarConstraint constraint1 = new ScalarConstraint(ConstraintType.INEQ, constraintFunc1, null);
        final List<ScalarConstraint> constraintList1 = new ArrayList<>();
        constraintList1.add(constraint1);

        final Vector2ScalarFunc constraintFunc2 = (x1, arg) -> x1[0] - 2;
        final ScalarConstraint constraint2 = new ScalarConstraint(ConstraintType.INEQ, constraintFunc2, null);
        final List<ScalarConstraint> constraintList2 = new ArrayList<>();
        constraintList2.add(constraint2);

        final Vector2ScalarFunc constraintFunc3 = (x1, arg) -> x1[0] + 1;
        final ScalarConstraint constraint3 = new ScalarConstraint(ConstraintType.INEQ, constraintFunc3, null);
        final List<ScalarConstraint> constraintList3 = new ArrayList<>();
        constraintList3.add(constraint3);
        constraintList3.add(constraint1);

        Slsqp slsqp = new Slsqp(inputFunc, null, x);
        OptimizeResult result = slsqp.minimize_slsqp_with_scalar_constraints(
            null,
            constraintList1,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        double[] resX = result.x;
        double[] expected = {0};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(result.success);

        x = new double[] {-10};
        slsqp = new Slsqp(inputFunc, null, x);
        result = slsqp.minimize_slsqp_with_scalar_constraints(
            null,
            constraintList2,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        resX = result.x;
        expected = new double[]{2};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    // Vector tests

    @Test
    public void test_minimize_equality_approximated()
    {
        final Slsqp slsqp = new Slsqp(new TestUtil.Fun(), null, new double[]{-1, 1});
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null, -1);
        constraints.add(constraint);
        final OptimizeResult result = slsqp.minimize_slsqp_with_vector_constraints(
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {1, 1};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void test_minimize_equality_given()
    {
        final double[] x = new double[] {-1, 1};
        final Slsqp slsqp = new Slsqp(new TestUtil.Fun(), new TestUtil.Jac(), x);
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null, -1);
        constraints.add(constraint);
        final OptimizeResult result = slsqp.minimize_slsqp_with_vector_constraints(
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {1, 1};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(result.success);
    }
}
