package com.example.slsqp.constraints;

import com.example.slsqp.Jacobian;
import com.example.slsqp.functions.Vector2ScalarFunc;

public class ScalarConstraint
{
    private ConstraintType constraintType;
    private Vector2ScalarFunc constraintFunc;
    private double[] jacobian;
    private double[] arg;

    public ScalarConstraint(
        ConstraintType constraintType,
        Vector2ScalarFunc constraintFunc,
        double[] jacobian,
        double... arg)
    {
        this.constraintType = constraintType;
        this.constraintFunc = constraintFunc;
        this.jacobian = jacobian;
        this.arg = arg;
    }

    public ConstraintType getConstraintType()
    {
        return constraintType;
    }

    public double[] getJacobian(double[] x)
    {
        if (jacobian == null)
        {
            return Jacobian.approx_jacobian(x, constraintFunc, arg);
        }
        else
        {
            return jacobian;
        }
    }

    public double apply(double[] x)
    {
        return constraintFunc.apply(x, arg);
    }
}
