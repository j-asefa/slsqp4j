package sci4j.optimize.slsqp.functions;

import sci4j.optimize.slsqp.Jacobian;

public class WrappedScalarFunction
{
    private final Vector2ScalarFunc func;
    private final double[] arg;

    public WrappedScalarFunction(Vector2ScalarFunc func, double... arg)
    {
        this.func = func;
        this.arg = arg;
    }

    public double apply(double[] x)
    {
        return func.apply(x, arg);
    }

    public double[] approxJacobian(double[] x)
    {
        return Jacobian.approxJacobian(x, this.func, arg);
    }
}
