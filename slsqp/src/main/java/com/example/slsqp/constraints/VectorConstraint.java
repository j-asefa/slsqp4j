package com.example.slsqp.constraints;

import com.example.slsqp.Jacobian;
import com.example.slsqp.functions.Vector2MatrixFunc;
import com.example.slsqp.functions.Vector2VectorFunc;

public class VectorConstraint
{
    private double[] arg;
    private ConstraintType constraintType;
    private Vector2VectorFunc constraintFunc;
    private Vector2MatrixFunc jacobian;

    public VectorConstraint(
        ConstraintType constraintType,
        Vector2VectorFunc constraintFunc,
        Vector2MatrixFunc jacobian,
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

    public double[][] getJacobian(double[] x)
    {
        if (jacobian == null)
        {
            return Jacobian.approx_jacobian(x, constraintFunc, arg);
        }
        else
        {
            return jacobian.apply(x, arg);
        }
    }

    public double[] apply(double[] x)
    {
        return constraintFunc.apply(x, arg);
    }
}
