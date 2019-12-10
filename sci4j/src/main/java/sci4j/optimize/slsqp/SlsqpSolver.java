package sci4j.optimize.slsqp;

import sci4j.optimize.slsqp.constraints.ScalarConstraint;
import sci4j.optimize.slsqp.constraints.VectorConstraint;
import sci4j.optimize.slsqp.functions.CallBackFunc;

import java.util.List;

public interface SlsqpSolver
{
    OptimizeResult minimize_slsqp_with_scalar_constraints(
        double[][] bounds,
        List<ScalarConstraint> scalarConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc,
        double... objectiveFunctionArgs
        );

    OptimizeResult minimize_slsqp_with_vector_constraints(
        double[][] bounds,
        List<VectorConstraint> vectorConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc,
        double... objectiveFunctionArgs
        );
}
