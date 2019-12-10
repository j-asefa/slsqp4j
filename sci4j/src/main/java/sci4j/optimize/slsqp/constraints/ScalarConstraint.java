package sci4j.optimize.slsqp.constraints;

import sci4j.optimize.slsqp.Jacobian;
import sci4j.optimize.slsqp.functions.Vector2ScalarFunc;
import sci4j.optimize.slsqp.functions.Vector2VectorFunc;

public class ScalarConstraint
{
    private ConstraintType constraintType;
    private Vector2ScalarFunc constraintFunc;
    private Vector2VectorFunc jacobian;
    private double[] arg;

    public ScalarConstraint(
        ConstraintType constraintType,
        Vector2ScalarFunc constraintFunc,
        Vector2VectorFunc jacobian,
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
            return Jacobian.approxJacobian(x, constraintFunc, arg);
        }
        else
        {
            return jacobian.apply(x, arg);
        }
    }

    public double apply(double[] x)
    {
        return constraintFunc.apply(x, arg);
    }
}
