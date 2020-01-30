package slsqp.optimize;

import org.junit.jupiter.api.Test;
import slsqp.optimize.constraints.ConstraintType;
import slsqp.optimize.constraints.ScalarConstraint;
import slsqp.optimize.constraints.VectorConstraint;
import slsqp.optimize.functions.Vector2MatrixFunc;
import slsqp.optimize.functions.Vector2ScalarFunc;
import slsqp.optimize.functions.Vector2VectorFunc;

import static org.junit.jupiter.api.Assertions.assertTrue;

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

        final ScalarConstraint constraint = new ScalarConstraint.ScalarConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.TestConstraintFunc())
            .build();

        final Vector2ScalarFunc objectiveFunc = new TestUtil.TestInputFunc();
        final double tolerance = 1.0E-6;
        final int maxIter = 100;
        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunc)
            .withBounds(bounds)
            .addConstraint(constraint)
            .withTolerance(tolerance)
            .withMaxIterations(maxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);

        assertTrue(Math.abs(result.x[0] - 2.5) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - 2.5) < TestUtil.ERROR);
        assertTrue(result.success);
    }


    @Test
    public void testASymmetricInput()
    {
        final double[] xl = new double[]{0, 0};
        final double[] xu = new double[]{100, 5};
        final double[] x = new double[] {0.2, 0.9};
        final double[][] bounds = new double[][] {xl, xu};

        final ScalarConstraint constraint = new ScalarConstraint.ScalarConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.TestConstraintFunc())
            .build();

        final Vector2ScalarFunc objectiveFunc = new TestUtil.TestInputFunc();
        final double tolerance = 1.0E-6;
        final int maxIter = 100;

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunc)
            .withBounds(bounds)
            .addConstraint(constraint)
            .withTolerance(tolerance)
            .withMaxIterations(maxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        assertTrue(Math.abs(result.x[0] - 1.99840144e-14) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - 5) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void testMinimizeEqualityGivenConsScalar()
    {
        final double[] x = new double[] {-1, 1};
        final ScalarConstraint constraint = new ScalarConstraint.ScalarConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.FeconScalar())
            .withJacobian(new TestUtil.FprimeEconScalar())
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] expected = {1, 1};
        assertTrue(Math.abs(result.x[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(result.x[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void testScalarConstraints()
    {
        final double[] x = new double[] {3};
        final Vector2ScalarFunc objectiveFunction = (x1, arg) -> Math.pow(x1[0], 2);
        final Vector2ScalarFunc constraintFunc = (x1, arg) -> x1[0] - 1;

        final ScalarConstraint constraint = new ScalarConstraint.ScalarConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(constraintFunc)
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunction)
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {1};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void testInfeasibleInitial()
    {
        double[] x = new double[] {10};

        final Vector2ScalarFunc objectiveFunction = (x1, arg) -> Math.pow(x1[0], 2) - 2 * x1[0] + 1;

        final Vector2ScalarFunc constraintFunc1 = (x1, arg) -> 0 - x1[0];
        final ScalarConstraint constraint1 = new ScalarConstraint.ScalarConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(constraintFunc1)
            .build();

        final Vector2ScalarFunc constraintFunc2 = (x1, arg) -> x1[0] - 2;
        final ScalarConstraint constraint2 = new ScalarConstraint.ScalarConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(constraintFunc2)
            .build();

        final Slsqp slsqp1 = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunction, -1)
            .addConstraint(constraint1)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        OptimizeResult result = slsqp1.minimize(x);

        double[] resX = result.x;
        double[] expected = {0};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(result.success);

        x = new double[] {-10};

        final Slsqp slsqp2 = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunction, -1)
            .addConstraint(constraint2)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();
        result = slsqp2.minimize(x);
        resX = result.x;
        expected = new double[]{2};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    // Vector tests

    @Test
    public void testMinimizeEqualityApproximated()
    {
        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.Fecon(), -1)
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(new double[] {-1, 1});
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
        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.Fecon(), -1)
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
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
        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(new TestUtil.Fieqcon(), -1)
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {2, 1};
        assertTrue(Math.abs(resX[0] - expected[0]) < 1.0E-3);
        assertTrue(Math.abs(resX[1] - expected[1]) < 1.0E-3);
        assertTrue(result.success);
    }

    @Test
    public void testMinimizeInequalityGiven2()
    {
        final double[] x = new double[]{-1, 1};
        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(new TestUtil.Fieqcon2(), -1)
            .withJacobian(new TestUtil.FprimeIeqcon2())
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {2, 1};
        assertTrue(Math.abs(resX[0] - expected[0]) < 1.0E-3);
        assertTrue(Math.abs(resX[1] - expected[1]) < 1.0E-3);
        assertTrue(result.success);
    }

    @Test
    public void testMinimizeBoundedEqualityConstraint()
    {
        final double[] x = new double[]{-1, 1};

        final double[] lowerBounds = new double[] {-0.8, -1};
        final double[] upperBounds = new double[] {1, 0.8};
        final double[][] bounds = new double[][] {lowerBounds, upperBounds};

        final Vector2MatrixFunc constraintJac = (x1, arg) ->
            Jacobian.transpose(new double[][] {new TestUtil.FprimeEcon().apply(x1, arg)});

        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.Fecon(), -1)
            .withJacobian(constraintJac)
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .withBounds(bounds)
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {0.8, 0.8};
        assertTrue(Math.abs(resX[0] - expected[0]) < 1.0E-3);
        assertTrue(Math.abs(resX[1] - expected[1]) < 1.0E-3);
        assertTrue(-0.8 <= resX[0] && resX[0] <= 1);
        assertTrue(-1 <= resX[1] && resX[1] <= 0.8);
        assertTrue(result.success);
    }

    @Test
    public void test4DVector()
    {
        final double[] x = new double[]{0.84, 1.3, -0.992, 0.18};
        final Vector2VectorFunc constraintFunc = (x1, arg) -> x1;

        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(constraintFunc)
            .build();

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

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunction, -1)
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);

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
    public void test4DVectorWithBounds()
    {
        final double[] x = new double[]{0.84, 1.3, -0.992, 0.18};

        final double[] lowerBounds = new double[] {0.5, -1, -1.4, -2.2};
        final double[] upperBounds = new double[] {1, 1.9, 1.3, 0.8};
        final double[][] bounds = new double[][] {lowerBounds, upperBounds};

        final Vector2VectorFunc constraintFunc = (x1, arg) -> x1;

        final VectorConstraint constraint = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(constraintFunc)
            .build();

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

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(objectiveFunction, -1)
            .withBounds(bounds)
            .addConstraint(constraint)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {1.00000000e+00, -1.33226763e-15, -1.33226763e-15, -4.99600361e-16};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[2] - expected[2]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[3] - expected[3]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void testMultipleConstraints()
    {
        final double[] x = new double[]{-1.4, 0.9};

        final double[] lowerBounds = new double[] {-0.8, -1};
        final double[] upperBounds = new double[] {1, 0.8};
        final double[][] bounds = new double[][] {lowerBounds, upperBounds};

        final Vector2VectorFunc constraintFunc = (x1, arg) -> new double[]{2 * x1[0] - 3 * x1[1]};

        final VectorConstraint constraint1 = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(constraintFunc)
            .build();

        final VectorConstraint constraint2 = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.EQ)
            .withConstraintFunction(new TestUtil.Fecon())
            .build();

        final VectorConstraint constraint3 = new VectorConstraint.VectorConstraintBuilder()
            .withConstraintType(ConstraintType.INEQ)
            .withConstraintFunction(new TestUtil.Fieqcon2())
            .build();

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .withBounds(bounds)
            .addConstraint(constraint1)
            .addConstraint(constraint2)
            .addConstraint(constraint3)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {-2.38418578e-08, -2.38418574e-08};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(result.success);
    }

    @Test
    public void testUnconstrained()
    {
        final double[] x = new double[]{-1.4, 0.9};

        final double[] lowerBounds = new double[] {-0.8, -1};
        final double[] upperBounds = new double[] {1, 0.8};
        final double[][] bounds = new double[][] {lowerBounds, upperBounds};

        final Slsqp slsqp = new Slsqp.SlsqpBuilder()
            .withObjectiveFunction(new TestUtil.Fun(), -1)
            .withJacobian(new TestUtil.Jac())
            .withBounds(bounds)
            .withTolerance(defaultTol)
            .withMaxIterations(defaultMaxIter)
            .build();

        final OptimizeResult result = slsqp.minimize(x);
        final double[] resX = result.x;
        final double[] expected = {1, 0.5};
        assertTrue(Math.abs(resX[0] - expected[0]) < TestUtil.ERROR);
        assertTrue(Math.abs(resX[1] - expected[1]) < TestUtil.ERROR);
        assertTrue(result.success);
    }
}
