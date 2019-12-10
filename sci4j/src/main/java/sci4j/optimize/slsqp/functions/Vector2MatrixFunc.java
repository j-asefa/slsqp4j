package sci4j.optimize.slsqp.functions;

@FunctionalInterface
public interface Vector2MatrixFunc
{
    double[][] apply(double[] x, double... arg);
}
