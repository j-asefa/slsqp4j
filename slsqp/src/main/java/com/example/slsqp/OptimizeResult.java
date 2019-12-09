package com.example.slsqp;

public class OptimizeResult
{
    public final double[] x;
    public final double fx;
    public final double[] jac;
    public final int status;
    public final int nfev;
    public final int njev;
    public final int exit_mode;
    public final boolean success;

    public OptimizeResult(
        double[] x,
        double fx,
        double[] jac,
        int status,
        int nfev,
        int njev,
        int exit_mode,
        boolean success
    )
    {
        this.x = x;
        this.fx = fx;
        this.jac = jac;
        this.status = status;
        this.nfev = nfev;
        this.njev = njev;
        this.exit_mode = exit_mode;
        this.success = success;
    }
}
