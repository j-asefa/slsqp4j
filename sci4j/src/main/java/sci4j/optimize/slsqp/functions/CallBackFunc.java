package sci4j.optimize.slsqp.functions;

@FunctionalInterface
public interface CallBackFunc
{
    void callback(double[] x);
}
