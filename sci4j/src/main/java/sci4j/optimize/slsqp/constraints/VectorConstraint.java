package sci4j.optimize.slsqp.constraints;

import sci4j.optimize.slsqp.Jacobian;
import sci4j.optimize.slsqp.functions.Vector2MatrixFunc;
import sci4j.optimize.slsqp.functions.Vector2VectorFunc;

public class VectorConstraint
{
    private double[] arg;
    private ConstraintType constraintType;
    private Vector2VectorFunc constraintFunc;
    private Vector2MatrixFunc jacobian;

    /**
     * Constructs a Vector Constraint to be used by the SLSQP solver
     *
     * @param constraintType type of this constraint (EQ, INEQ).
     * @param constraintFunc vector to vector function representing the constraint function to apply
     * @param jacobian function outputting the jacobian of the constraint function. Note that the solver expects the
     *                 Jacobian to be in column-major order, so a call to {@link Jacobian#transpose(double[][])}
     *                 should be made before passing the jacobian to this function.
     * @param arg optional arguments to the constraint function
     */
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

    /**
     * Constructs a VectorConstraint that does not take a Jacobian for the constraint function, thus the numerical
     * approximation for the Jacobian is used.
     *
     * @param constraintType type of this constraint (EQ, INEQ).
     * @param constraintFunc vector to vector function representing the constraint function to apply
     * @param arg optional arguments to the constraint function
     */
    public VectorConstraint(
        ConstraintType constraintType,
        Vector2VectorFunc constraintFunc,
        double... arg)
    {
        this.constraintType = constraintType;
        this.constraintFunc = constraintFunc;
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
            return Jacobian.approxJacobian(x, constraintFunc, arg);
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
