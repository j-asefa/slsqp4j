package slsqp.optimize.functions;

@FunctionalInterface
public interface Vector2VectorFunc
{
    double[] apply(double[] x, double... arg);
}
