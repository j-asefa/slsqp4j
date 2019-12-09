package com.example.slsqp;

import com.example.slsqp.constraints.ScalarConstraint;
import com.example.slsqp.constraints.VectorConstraint;
import com.example.slsqp.functions.CallBackFunc;

import java.util.List;

public interface SlsqpSolver
{
    OptimizeResult minimize_slsqp_with_scalar_constraints(
        double[][] bounds,
        double sign,
        List<ScalarConstraint> scalarConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc);

    OptimizeResult minimize_slsqp_with_vector_constraints(
        double[][] bounds,
        double sign,
        List<VectorConstraint> vectorConstraintList,
        double tolerance,
        int maxIterations,
        CallBackFunc callBackFunc);
}
