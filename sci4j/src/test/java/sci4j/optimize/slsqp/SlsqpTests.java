package sci4j.optimize.slsqp;

import sci4j.optimize.slsqp.functions.Vector2MatrixFunc;
import sci4j.optimize.slsqp.functions.Vector2ScalarFunc;
import sci4j.optimize.slsqp.constraints.ConstraintType;
import sci4j.optimize.slsqp.constraints.ScalarConstraint;
import sci4j.optimize.slsqp.constraints.VectorConstraint;
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
        final double[] xl = new double[]{0, 0};
        final double[] xu = new double[]{100, 5};

        final double[][] bounds = new double[][] {xl, xu};

        final double[] x = new double[]{0, 0};
        final ScalarConstraint scalarConstraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.TestConstraintFunc(),
            null
        );
        final List<ScalarConstraint> constraintList = new ArrayList<>();
        constraintList.add(scalarConstraint);
        final Vector2ScalarFunc inputFunc = new TestUtil.TestInputFunc();
        final double tolerance = 1.0E-6;
        final int maxIter = 100;
        final OptimizeResult result = Slsqp.minimizeSlsqpWithScalarConstraints(
            inputFunc,
            null,
            x,
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
        final double[] xl = new double[]{0, 0};
        final double[] xu = new double[]{100, 5};

        final double[][] bounds = new double[][] {xl, xu};

        final double[] x = new double[] {0.2, 0.9};
        final ScalarConstraint scalarConstraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.TestConstraintFunc(),
            null
        );
        final List<ScalarConstraint> constraintList = new ArrayList<>();
        constraintList.add(scalarConstraint);

        final Vector2ScalarFunc inputFunc = new TestUtil.TestInputFunc();
        final double tolerance = 1.0E-6;
        final int maxIter = 100;
        final OptimizeResult result = Slsqp.minimizeSlsqpWithScalarConstraints(
            inputFunc,
            null,
            x,
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
    public void testMinimizeEqualityGivenConsScalar()
    {
        final double[] x = new double[] {-1, 1};
        final List<ScalarConstraint> constraints = new ArrayList<>();
        final ScalarConstraint constraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.FeconScalar(),
            new TestUtil.FprimeEconScalar());
        constraints.add(constraint);
        final OptimizeResult result = Slsqp.minimizeSlsqpWithScalarConstraints(
            new TestUtil.Fun(),
            new TestUtil.Jac(),
            x,
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
    public void testScalarConstraints()
    {
        final double[] x = new double[] {3};

        final Vector2ScalarFunc inputFunc = (x1, arg) -> Math.pow(x1[0], 2);

        final Vector2ScalarFunc constraintFunc = (x1, arg) -> x1[0] - 1;

        final List<ScalarConstraint> constraints = new ArrayList<>();
        final ScalarConstraint constraint = new ScalarConstraint(
            ConstraintType.INEQ,
            constraintFunc,
            null);
        constraints.add(constraint);
        final OptimizeResult result = Slsqp.minimizeSlsqpWithScalarConstraints(
            inputFunc,
            null,
            x,
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
    public void testInfeasibleInitial()
    {
        double[] x = new double[] {10};

        final Vector2ScalarFunc inputFunc = (x1, arg) -> Math.pow(x1[0], 2) - 2 * x1[0] + 1;

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

        OptimizeResult result = Slsqp.minimizeSlsqpWithScalarConstraints(
            inputFunc,
            null,
            x,
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
        result = Slsqp.minimizeSlsqpWithScalarConstraints(
            inputFunc,
            null,
            x,
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
    public void testMinimizeEqualityApproximated()
    {
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null, -1);
        constraints.add(constraint);
        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            new TestUtil.Fun(),
            null,
            new double[]{-1, 1},
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
    public void testMinimizeEqualityGiven()
    {
        final double[] x = new double[] {-1, 1};
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null, -1);
        constraints.add(constraint);
        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            new TestUtil.Fun(),
            new TestUtil.Jac(),
            x,
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
    public void testMinimizeInequalityGiven()
    {
        final double[] x = new double[]{-1, 1};
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.INEQ, new TestUtil.Fieqcon(), null, -1);
        constraints.add(constraint);
        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            new TestUtil.Fun(),
            new TestUtil.Jac(),
            x,
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {2, 1};
        assertTrue(Math.abs(resX[0] - expected[0]) < 1.0E-3);
        assertTrue(Math.abs(resX[1] - expected[1]) < 1.0E-3);
        assertTrue(result.success);
    }

    @Test
    public void testMinimizeInequalityGivenVectorConstraints()
    {
        final double[] x = new double[]{-1, 1};
        final List<VectorConstraint> constraints = new ArrayList<>();
        final VectorConstraint constraint = new VectorConstraint(
            ConstraintType.INEQ,
            new TestUtil.Fieqcon2(),
            new TestUtil.FprimeIeqcon2(),
            -1);
        constraints.add(constraint);
        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            new TestUtil.Fun(),
            new TestUtil.Jac(),
            x,
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {2, 1};
        assertTrue(Math.abs(resX[0] - expected[0]) < 1.0E-3);
        assertTrue(Math.abs(resX[1] - expected[1]) < 1.0E-3);
        assertTrue(result.success);
    }

    @Test
    public void testMinimizeBoundEqualityGiven2()
    {
        final double[] x = new double[]{-1, 1};

        final double[] lowerBounds = new double[] {-0.8, -1};
        final double[] upperBounds = new double[] {1, 0.8};
        final double[][] bounds = new double[][] {lowerBounds, upperBounds};

        final List<VectorConstraint> constraints = new ArrayList<>();
        final Vector2MatrixFunc constraintJac = (x1, arg) ->
            Jacobian.transpose(new double[][] {new TestUtil.FprimeEcon().apply(x1, arg)});

        final VectorConstraint constraint = new VectorConstraint(
            ConstraintType.EQ,
            new TestUtil.Fecon(),
            constraintJac,
            -1);
        constraints.add(constraint);

        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            new TestUtil.Fun(),
            new TestUtil.Jac(),
            x,
            bounds,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {0.8, 0.8};
        assertTrue(Math.abs(resX[0] - expected[0]) < 1.0E-3);
        assertTrue(Math.abs(resX[1] - expected[1]) < 1.0E-3);
        assertTrue(-0.8 <= resX[0] && resX[0] <= 1);
        assertTrue(-1 <= resX[1] && resX[1] <= 0.8);
        assertTrue(result.success);
    }
}
