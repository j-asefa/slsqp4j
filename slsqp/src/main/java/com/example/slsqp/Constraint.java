package com.example.slsqp;

public class Constraint
{
    private Vector2ScalarFunc func;
    private ConstraintType constraintType;

    public Constraint(Vector2ScalarFunc func, ConstraintType constraintType)
    {
        this.func = func;
        this.constraintType = constraintType;
    }

    public Vector2ScalarFunc getFunc()
    {
        return func;
    }

    public ConstraintType getConstraintType()
    {
        return constraintType;
    }
}
