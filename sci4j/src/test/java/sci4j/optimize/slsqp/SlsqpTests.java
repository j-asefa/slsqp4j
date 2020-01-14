package sci4j.optimize.slsqp;

import sci4j.optimize.slsqp.functions.Vector2MatrixFunc;
import sci4j.optimize.slsqp.functions.Vector2ScalarFunc;
import sci4j.optimize.slsqp.constraints.ConstraintType;
import sci4j.optimize.slsqp.constraints.ScalarConstraint;
import sci4j.optimize.slsqp.constraints.VectorConstraint;
import org.junit.jupiter.api.Test;
import sci4j.optimize.slsqp.functions.Vector2VectorFunc;

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
        final ScalarConstraint[] constraintList = new ScalarConstraint[1];
        constraintList[0] = scalarConstraint;
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
        final ScalarConstraint[] constraintList = new ScalarConstraint[1];
        constraintList[0] = scalarConstraint;

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
        final ScalarConstraint[] constraints = new ScalarConstraint[1];
        final ScalarConstraint constraint = new ScalarConstraint(
            ConstraintType.EQ,
            new TestUtil.FeconScalar(),
            new TestUtil.FprimeEconScalar());
        constraints[0] = constraint;

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

        final ScalarConstraint[] constraints = new ScalarConstraint[1];
        final ScalarConstraint constraint = new ScalarConstraint(
            ConstraintType.INEQ,
            constraintFunc,
            null);
        constraints[0] = constraint;

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
        final ScalarConstraint[] constraintList1 = new ScalarConstraint[1];
        constraintList1[0] = constraint1;

        final Vector2ScalarFunc constraintFunc2 = (x1, arg) -> x1[0] - 2;
        final ScalarConstraint constraint2 = new ScalarConstraint(ConstraintType.INEQ, constraintFunc2, null);
        final ScalarConstraint[] constraintList2 = new ScalarConstraint[1];
        constraintList2[0] = constraint2;

        final Vector2ScalarFunc constraintFunc3 = (x1, arg) -> x1[0] + 1;
        final ScalarConstraint constraint3 = new ScalarConstraint(ConstraintType.INEQ, constraintFunc3, null);
        final ScalarConstraint[] constraintList3 = new ScalarConstraint[2];
        constraintList3[0] = constraint3;
        constraintList3[1] = constraint1;

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
        final VectorConstraint[] constraints = new VectorConstraint[1];
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null, -1);
        constraints[0] = constraint;
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
        final VectorConstraint[] constraints = new VectorConstraint[1];
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.EQ, new TestUtil.Fecon(), null, -1);
        constraints[0] = constraint;

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
        final VectorConstraint[] constraints = new VectorConstraint[1];
        final VectorConstraint constraint = new VectorConstraint(ConstraintType.INEQ, new TestUtil.Fieqcon(), null, -1);
        constraints[0] = constraint;

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
        final VectorConstraint[] constraints = new VectorConstraint[1];
        final VectorConstraint constraint = new VectorConstraint(
            ConstraintType.INEQ,
            new TestUtil.Fieqcon2(),
            new TestUtil.FprimeIeqcon2(),
            -1);
        constraints[0] = constraint;

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

        final VectorConstraint[] constraints = new VectorConstraint[1];
        final Vector2MatrixFunc constraintJac = (x1, arg) ->
            Jacobian.transpose(new double[][] {new TestUtil.FprimeEcon().apply(x1, arg)});

        final VectorConstraint constraint = new VectorConstraint(
            ConstraintType.EQ,
            new TestUtil.Fecon(),
            constraintJac,
            -1);
        constraints[0] = constraint;

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

    @Test
    public void testLargerMatrix()
    {
        final double[] x = new double[]{0.84, 1.3, -0.992, 0.18};

        final VectorConstraint[] constraints = new VectorConstraint[1];
        final Vector2VectorFunc constraintFunc = (x1, arg) -> x1;

        final VectorConstraint constraint = new VectorConstraint(
            ConstraintType.INEQ,
            constraintFunc,
            null);
        constraints[0] = constraint;

        final Vector2ScalarFunc objectiveFunction = (x12, arg) ->
        {
            final double a = x12[0];
            final double b = x12[1];
            final double c = x12[2];
            final double d = x12[3];
            int sign = 1;
            if (arg != null && arg.length > 0)
            {
                sign = (int)arg[0];
            }
            return sign * (a * b * c * d + 2 * a - 2 * b + Math.pow(c, 2) + Math.pow(d, 2));
        };
        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            objectiveFunction,
            null,
            x,
            null,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {
            81869346598383174108963135153684611072.0000000000000000000000000000000000000,
            57919931473235792828070634541023232.0000000000000000000000000000000000000000000000000000000,
            53311942694712538609603367273781788672.0000000000000000000000000000000000000000000000000000000,
            3329199916413834719641054178443264.0000000000000000000000000000000000000000000000000000000};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[2] - expected[2]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[3] - expected[3]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void testLargerMatrixWithBounds()
    {
        final double[] x = new double[]{0.84, 1.3, -0.992, 0.18};

        final double[] lowerBounds = new double[] {0.5, -1, -1.4, -2.2};
        final double[] upperBounds = new double[] {1, 1.9, 1.3, 0.8};
        final double[][] bounds = new double[][] {lowerBounds, upperBounds};


        final VectorConstraint[] constraints = new VectorConstraint[1];
        final Vector2VectorFunc constraintFunc = (x1, arg) -> x1;

        final VectorConstraint constraint = new VectorConstraint(
            ConstraintType.INEQ,
            constraintFunc,
            null);
        constraints[0] = constraint;

        final Vector2ScalarFunc objectiveFunction = (x12, arg) ->
        {
            final double a = x12[0];
            final double b = x12[1];
            final double c = x12[2];
            final double d = x12[3];
            int sign = 1;
            if (arg != null && arg.length > 0)
            {
                sign = (int)arg[0];
            }
            return sign * (a * b * c * d + 2 * a - 2 * b + Math.pow(c, 2) + Math.pow(d, 2));
        };
        final OptimizeResult result = Slsqp.minimizeSlsqpWithVectorConstraints(
            objectiveFunction,
            null,
            x,
            bounds,
            constraints,
            defaultTol,
            defaultMaxIter,
            null,
            -1
        );
        final double[] resX = result.x;
        final double[] expected = {1.00000000e+00, -1.33226763e-15, -1.33226763e-15, -4.99600361e-16};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[2] - expected[2]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[3] - expected[3]) < TestUtil.ERROR);
        assertTrue(result.success);
    }
}
