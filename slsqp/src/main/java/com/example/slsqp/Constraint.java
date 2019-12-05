package com.example.slsqp;

public class Constraint
{
    private Func func;
    private ConstraintType constraintType;

    public Constraint(Func func, ConstraintType constraintType)
    {
        this.func = func;
        this.constraintType = constraintType;
    }

    public Func getFunc()
    {
        return func;
    }

    public ConstraintType getConstraintType()
    {
        return constraintType;
    }
}
