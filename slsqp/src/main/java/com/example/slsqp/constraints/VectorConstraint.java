package com.example.slsqp.constraints;

import com.example.slsqp.Jacobian;
import com.example.slsqp.Vector2VectorFunc;

public class VectorConstraint
{

    private ConstraintType constraintType;
    private Vector2VectorFunc constraintFunc;
    private double[][] jacobian;

    public VectorConstraint(
        ConstraintType constraintType,
        Vector2VectorFunc constraintFunc,
        double[][] jacobian)
    {
        this.constraintType = constraintType;
        this.constraintFunc = constraintFunc;
        this.jacobian = jacobian;
    }

    public ConstraintType getConstraintType()
    {
        return constraintType;
    }

    public double[][] getJacobian(double[] x)
    {
        if (jacobian == null)
        {
            return Jacobian.approx_jacobian(x, constraintFunc);
        }
        else
        {
            return jacobian;
        }
    }

    public Vector2VectorFunc getConstraintFunc()
    {
        return constraintFunc;
    }
}
