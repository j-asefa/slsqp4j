package com.example.slsqp.constraints;

import com.example.slsqp.Jacobian;
import com.example.slsqp.Vector2ScalarFunc;

public class ScalarConstraint
{
    private ConstraintType constraintType;
    private Vector2ScalarFunc constraintFunc;
    private double[] jacobian;
    private double arg = 0;

    public ScalarConstraint(
        ConstraintType constraintType,
        Vector2ScalarFunc constraintFunc,
        double[] jacobian,
        double arg)
    {
        this(constraintType, constraintFunc, jacobian);
        this.arg = arg;
    }

    public ScalarConstraint(
        ConstraintType constraintType,
        Vector2ScalarFunc constraintFunc,
        double[] jacobian)
    {
        this.constraintType = constraintType;
        this.constraintFunc = constraintFunc;
        this.jacobian = jacobian;
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

    public Vector2ScalarFunc getConstraintFunc()
    {
        return constraintFunc;
    }
}
