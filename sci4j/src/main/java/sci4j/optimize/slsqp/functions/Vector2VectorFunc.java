package sci4j.optimize.slsqp.functions;

@FunctionalInterface
public interface Vector2VectorFunc
{
    double[] apply(double[] x, double... arg);
}
