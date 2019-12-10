package sci4j.optimize.slsqp;

import sci4j.optimize.slsqp.functions.Vector2MatrixFunc;
import sci4j.optimize.slsqp.functions.Vector2ScalarFunc;
import sci4j.optimize.slsqp.functions.Vector2VectorFunc;

public class TestUtil
{
    public static final double ERROR = 1.0E-6;

    // this constraint and input function is taken from the example at
    // https://stackoverflow.com/questions/26882087/python-scipy-optimization-minimize-using-slsqp-showing-maximized-results
    public static class TestConstraintFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return x[0] + x[1] - 5;
        }
    }

    public static class TestInputFunc implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return x[0] * x[1];
        }
    }

    public static final class Fun implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] arr, double... arg)
        {
            if (arg.length > 1)
            {
                throw new IllegalArgumentException("optional argument must be one of {1, -1}");
            }

            final double x = arr[0];
            final double y = arr[1];
            double sign = 1;
            if (arg.length > 0)
            {
                sign = arg[0];
            }

            return sign * (2 * x * y + 2 * x - Math.pow(x, 2) - 2 * Math.pow(y, 2));
        }
    }

    public static final class Jac implements Vector2VectorFunc
    {
        // returns the analytical jacobian of Fun above.
        @Override
        public double[] apply(double[] arr, double... arg)
        {
            if (arg.length > 1)
            {
                throw new IllegalArgumentException("optional argument must be one of {1, -1}");
            }
            final double x = arr[0];
            final double y = arr[1];

            double sign = 1;
            if (arg.length > 0)
            {
                sign = arg[0];
            }
            final double dfdx = sign * (-2 * x + 2 * y + 2);
            final double dfdy = sign * (2 * x - 4 * y);
            return new double[] {dfdx, dfdy};
        }
    }

    public static final class Fecon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {x[0] - x[1]};
        }
    }

    public static final class FprimeEcon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {1, -1};
        }
    }

    public static final class FprimeEconScalar implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new FprimeEcon().apply(x, arg);
        }
    }

    public static final class FeconScalar implements Vector2ScalarFunc
    {
        @Override
        public double apply(double[] x, double... arg)
        {
            return new Fecon().apply(x, arg)[0];
        }
    }

    public static final class Fieqcon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {x[0] - x[1] - 1};
        }
    }

    public static final class FprimeIeqcon implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return new double[] {1, -1};
        }
    }

    public static final class Fieqcon2 implements Vector2VectorFunc
    {
        @Override
        public double[] apply(double[] x, double... arg)
        {
            return x;
        }
    }

    public static final class FprimeIeqcon2 implements Vector2MatrixFunc
    {
        public double[][] apply(double[] x, double... arg)
        {
            return identity(x.length);
        }

        private double[][] identity(int length)
        {
            final double[][] identityMat = new double[length][length];
            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < length; j++)
                {
                    if (i == j)
                    {
                        identityMat[i][j] = 1;
                    }
                    else
                    {
                        identityMat[i][j] = 0;
                    }
                }
            }
            return identityMat;
        }
    }
}
