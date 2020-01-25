package sci4j.optimize.slsqp;

public class OptimizeResult
{
    public final double[] x;
    public final double fx;
    public final double[] jac;
    public final int numIters;
    public final int exitMode;
    public final boolean success;
    public double[][] a;

    public OptimizeResult(
        double[] x,
        double fx,
        double[] jac,
        int numIters,
        int exitMode,
        boolean success,
        double[][] a
    )
    {
        this.x = x;
        this.fx = fx;
        this.jac = jac;
        this.numIters = numIters;
        this.exitMode = exitMode;
        this.success = success;
        this.a = a;
    }
}
