package sci4j.optimize.slsqp.functions;

import sci4j.optimize.slsqp.Jacobian;

public class WrappedVector2ScalarFunction
{
    private final Vector2ScalarFunc func;
    private final double[] arg;
    private final Vector2VectorFunc jacobian;

    public WrappedVector2ScalarFunction(Vector2ScalarFunc func, Vector2VectorFunc jacobian, double... arg)
    {
        this.func = func;
        this.jacobian = jacobian;
        this.arg = arg;
    }

    public double apply(double[] x)
    {
        return func.apply(x, arg);
    }

    public double[] getJacobian(double[] x)
    {
        if (jacobian == null)
        {
            return Jacobian.approxJacobian(x, this.func, arg);
        }
        else
        {
            return jacobian.apply(x, arg);
        }
    }
}
