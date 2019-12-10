package sci4j.optimize.slsqp;

import sci4j.optimize.slsqp.functions.Vector2ScalarFunc;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class JacobianTest
{
    @Test
    public void testJacobian()
    {
        final TestFunc testFunc = new TestFunc();

        final double[] inputArr = new double[]{-1, -1};
        final double[] fprime = Jacobian.approxJacobian(inputArr, testFunc);

        final double[] exactJac = testFunc.exactJac(inputArr);
        for (int i = 0; i < inputArr.length; i++)
        {
            assertTrue(Math.abs(fprime[i] - exactJac[i]) < TestUtil.ERROR);
        }
    }

    private static final class TestFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            final double x1 = x[0];
            final double x2 = x[1];
            return 2 * x1 * x2 + 2 * x1 - Math.pow(x1, 2) - 2 * Math.pow(x2, 2);
        }

        // returns the analytical jacobian of func
        double[] exactJac(double[] x)
        {
            final double x1 = x[0];
            final double x2 = x[1];
            return new double[] {-2 * x1 + 2 * x2 + 2, 2 * x1 - 4 * x2};
        }
    }
}
