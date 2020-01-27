package slsqp.optimize.functions;

@FunctionalInterface
public interface CallBackFunc
{
    void callback(double[] x);
}
