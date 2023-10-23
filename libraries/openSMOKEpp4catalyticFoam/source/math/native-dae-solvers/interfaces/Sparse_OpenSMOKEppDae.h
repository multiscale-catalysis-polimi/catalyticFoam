/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2015, 2014, 2013  Alberto Cuoci                          |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_Solve_Sparse_OpenSMOKEppDae_H
#define OpenSMOKE_Solve_Sparse_OpenSMOKEppDae_H

#include "math/native-dae-solvers/MultiValueSolver"

namespace DaeSMOKE
{
	template<typename Object, typename System>
	int Solve_Sparse_OpenSMOKEppDae(Object* object, const DaeSMOKE::DaeSolver_Parameters& parameters, const double t0, const double tEnd)
	{
		std::cout << "Sparse DAE solution (OpenSMOKE++)..." << std::endl;

		const unsigned int neq = object->NumberOfEquations();

		// Initial conditions
		Eigen::VectorXd yInitial(neq);
		object->UnknownsVector(yInitial.data());

		// Differential vs Algebraic
		std::vector<OpenSMOKE::EquationType> equation_type(neq);
		{
			Eigen::VectorXd temp(neq);
			object->AlgebraicDifferentialVector(temp.data());
			for (unsigned int i = 0; i < neq; i++)
			if (int(temp(i)) == 1)
				equation_type[i] = OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL;
			else
				equation_type[i] = OpenSMOKE::EQUATION_TYPE_ALGEBRAIC;
		}

		// Recognize the sparsity pattern
		std::vector<unsigned int> rows;
		std::vector<unsigned int> cols;
		object->SparsityPattern(rows, cols);

		typedef DaeSMOKE::KernelSparse<System> sparseDae;
		typedef DaeSMOKE::MethodGear<sparseDae> methodGear;
		DaeSMOKE::MultiValueSolver<methodGear> dae_solver;

		// Set initial conditions
		dae_solver.assign(object);
		dae_solver.SetSparsityPattern(neq, rows, cols, false);
		dae_solver.SetLinearAlgebraSolver(parameters.jacobian_solver());
		dae_solver.SetPreconditioner(parameters.preconditioner());
		dae_solver.SetInitialConditions(t0, yInitial, equation_type);
		
		// Minimum Constraints
		if (parameters.minimum_constraints() == true)
		{
			Eigen::VectorXd xMin(neq);
			object->MinimumUnknownsVector(xMin.data());
			dae_solver.SetMinimumValues(xMin);
		}

		// Maximum constraints
		if (parameters.maximum_constraints() == true)
		{
			Eigen::VectorXd xMax(neq);
			object->MaximumUnknownsVector(xMax.data());
			dae_solver.SetMaximumValues(xMax);
		}

		// Relative tolerance
		if (parameters.relative_tolerances().size() == neq)
			dae_solver.SetRelativeTolerances(parameters.relative_tolerances());
		else
			dae_solver.SetRelativeTolerances(parameters.relative_tolerance());

		// Absolute tolerance
		if (parameters.absolute_tolerances().size() == neq)
			dae_solver.SetAbsoluteTolerances(parameters.absolute_tolerances());
		else
			dae_solver.SetAbsoluteTolerances(parameters.absolute_tolerance());

		// Maximum order
		if (parameters.maximum_order() > 0)
			dae_solver.SetMaximumOrder(parameters.maximum_order());

		// Initial step
		if (parameters.initial_step() > 0.)
			dae_solver.SetFirstStepSize(parameters.initial_step());

		// Maximum number of steps
		if (parameters.maximum_number_of_steps() > 0)
			dae_solver.SetMaximumNumberOfSteps(parameters.maximum_number_of_steps());

		// Minimum step
		if (parameters.minimum_step() > 0.)
			dae_solver.SetMinimumStepSize(parameters.minimum_step());

		// Maximum step
		if (parameters.maximum_step() > 0.)
			dae_solver.SetMaximumStepSize(parameters.maximum_step());

		// Minimum sum of yp (mean)
		if (parameters.minimum_yp() > 0.)
			dae_solver.SetStopConditionMaximumYPrimeNorm1(parameters.minimum_yp()*neq);

		// Verbose
		if (parameters.verbosity_level() > 0)
			dae_solver.SetPrint(true);
		else
			dae_solver.SetPrint(false);

		// Solve the system
		double timeStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		DaeSMOKE::DaeStatus status = dae_solver.Solve(tEnd);
		double timeEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		// Check the solution
		if (status > 0)
		{
			std::string message("OpenSMOKE++ Dae System successfully solved: ");

			if (status == 1)		message += "INITIALIZATION_STATE";
			else if (status == 2)	message += "CONTINUATION_STATE";
			else if (status == 10)	message += "INTEGRATION_STOPPED_BEFORE_RECALCULATING_JACOBIAN";
			else if (status == 11)	message += "INTEGRATION_STOPPED_WHEN_SUM_ABS_Y1_IS_LESS_THAN";
			else if (status == 11)	message += "MAX_NUMBER_OF_STEPS_REACHED";

			std::cout << message << std::endl;

			Eigen::VectorXd yp(neq);
			dae_solver.FirstOrderDerivatives(yp);
			double sum_yp = 0.;
			for (unsigned int i = 0; i < neq; i++)
				sum_yp += std::fabs(yp(i));

			std::cout << std::endl;
			std::cout << " * CPU time (s):                   " << timeEnd - timeStart << std::endl;
			std::cout << " * number of steps:                " << dae_solver.numberOfSteps() << std::endl;
			std::cout << " * number of functions:            " << dae_solver.numberOfFunctionCalls() << std::endl;
			std::cout << " * number of solutions:            " << dae_solver.numberOfLinearSystemSolutions() << std::endl;
			std::cout << " * number of Jacobians:            " << dae_solver.numberOfJacobianEvaluations() << std::endl;
			std::cout << " * number of factorizations:       " << dae_solver.numberOfMatrixFactorizations() << std::endl;
			std::cout << " * number of functions (Jacobian): " << dae_solver.numberOfFunctionCallsForJacobian() << std::endl;
			std::cout << " * last order:                     " << dae_solver.lastOrderUsed() << std::endl;
			std::cout << " * last step size:                 " << std::scientific << dae_solver.lastStepUsed() << std::endl;
			std::cout << " * mean y':                        " << std::scientific << sum_yp / double(neq) << std::endl;
			std::cout << std::endl;

			Eigen::VectorXd ysol(neq);
			dae_solver.Solution(ysol);
			object->CorrectedUnknownsVector(ysol.data());
		}
		else
		{
			std::string message("OpenSMOKE++ Dae Solver Error: ");
			if (status == -2)	message += "TOO_STRICT_TOLERANCES";
			else if (status == -3)	message += "ILLEGAL_MAX_INDEPENDENT_VARIABLE";
			else if (status == -4)	message += "MAX_NUMBER_ERRORTEST_FAILURES";
			else if (status == -5)	message += "MAX_NUMBER_CONVERGENCETEST_FAILURES";
			else if (status == -6)	message += "TOO_SMALL_STEP_SIZE";
			else if (status == -7)	message += "YOU_MUST_USE_TCRITIC_STATE";
			else if (status == -8)	message += "ILLEGAL_CONTINUATION_REQUEST";
			else if (status == -9)	message += "ILLEGAL_CONSTRAINTS";
			else if (status == -10)	message += "EXCEPTION_HANDLING_STOP";
			else if (status == -11)	message += "DEINITIALIZE_STATE";
			else if (status == -12)	message += "YOU_CANNOT_OVERSHOOT_TCRITIC";

			std::cout << message << std::endl;

			object->CorrectedUnknownsVector(yInitial.data());

			if (status == -10)
			{
				std::cout << std::endl;
				std::cout << " * CPU time (s):                   " << timeEnd - timeStart << std::endl;
				std::cout << " * number of steps:                " << dae_solver.numberOfSteps() << std::endl;
				std::cout << " * number of functions:            " << dae_solver.numberOfFunctionCalls() << std::endl;
				std::cout << " * number of solutions:            " << dae_solver.numberOfLinearSystemSolutions() << std::endl;
				std::cout << " * number of Jacobians:            " << dae_solver.numberOfJacobianEvaluations() << std::endl;
				std::cout << " * number of factorizations:       " << dae_solver.numberOfMatrixFactorizations() << std::endl;
				std::cout << " * number of functions (Jacobian): " << dae_solver.numberOfFunctionCallsForJacobian() << std::endl;
				std::cout << " * last order:                     " << dae_solver.lastOrderUsed() << std::endl;
				std::cout << " * last step size:                 " << std::scientific << dae_solver.lastStepUsed() << std::endl;
				std::cout << std::endl;
			}
		}

		return status;
	}
}

#endif // OpenSMOKE_Solve_Sparse_OpenSMOKEppDae_H
