/*----------------------------------------------------------------------*\
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

namespace DaeSMOKE
{
	template<typename Object, typename System>
	int Solve_TridiagonalBlock_BzzDae(Object* object, BzzDaeSparseObject& dae_object, const DaeSMOKE::DaeSolver_Parameters& parameters, const double t0, const double tEnd)
	{
		std::cout << "Tridiagonal DAE solution (BzzDae)..." << std::endl;

		const unsigned int neq = object->NumberOfEquations();

		// Initial conditions
		BzzVector yInitial(neq);
		object->UnknownsVector(yInitial.GetHandle());

		// Differential vs Algebraic
		BzzVectorInt inDerAlg(neq);
		{
			BzzVector temp(neq);
			object->AlgebraicDifferentialVector(temp.GetHandle());
			for (int i = 1; i <= object->NumberOfEquations(); i++)
				inDerAlg[i] = int(temp[i]);
		}

		System dae_system;
		dae_system.assign(object);
		dae_object(yInitial, t0, inDerAlg, &dae_system, object->BlockDimensions());

		// Minimum Constraints
		if (parameters.minimum_constraints() == true)
		{
			BzzVector xMin(neq);
			object->MinimumUnknownsVector(xMin.GetHandle());
			dae_object.SetMinimumConstraints(xMin);
		}

		// Maximum constraints
		if (parameters.maximum_constraints() == true)
		{
			BzzVector xMax(neq);
			object->MaximumUnknownsVector(xMax.GetHandle());
			dae_object.SetMaximumConstraints(xMax);
		}

		// Relative tolerances
		if (parameters.relative_tolerances().size() == neq)
		{
			BzzVector reltol(neq);
			for (unsigned int i = 1; i <= neq; i++)
				reltol[i] = parameters.relative_tolerances()(i - 1);
			dae_object.SetTolRel(reltol);
		}
		else
		{
			dae_object.SetTolRel(parameters.relative_tolerance());
		}

		// Absolute tolerances
		if (parameters.absolute_tolerances().size() == neq)
		{
			BzzVector abstol(neq);
			for (unsigned int i = 1; i <= neq; i++)
				abstol[i] = parameters.absolute_tolerances()(i - 1);
			dae_object.SetTolAbs(abstol);
		}
		else
		{
			dae_object.SetTolAbs(parameters.absolute_tolerance());
		}

		// Maximum order
		if (parameters.maximum_order() > 0)
			dae_object.SetMaxOrder(parameters.maximum_order());

		// Initial step
		if (parameters.initial_step() > 0.)
			dae_object.SetH0(parameters.initial_step());

		// Maximum number of steps
		if (parameters.maximum_number_of_steps() > 0)
			dae_object.SetMaxStep(parameters.maximum_number_of_steps());

		// Minimum step
		if (parameters.minimum_step() > 0.)
			dae_object.SetHMin(parameters.minimum_step());

		// Maximum step
		if (parameters.maximum_step() > 0.)
			dae_object.SetHMax(parameters.maximum_step());

		// Max number Jacobians
		if (parameters.maximum_number_jacobians() > 0)
			dae_object.StopIntegrationBeforeRecalcuatingJacobian(parameters.maximum_number_jacobians());

		// Minimum sum of yp (mean)
		if (parameters.minimum_yp() > 0.)
			dae_object.StopIntegrationWhenSumAbsY1IsLessThan(parameters.minimum_yp()*neq);

		// Verbose
		if (parameters.verbosity_level() > 0)
			dae_object.StepPrint(DaePrint);

		// Solving the DAE system
		double timeStart = BzzGetCpuTime();
		dae_object(tEnd);
		double timeEnd = BzzGetCpuTime();

		const int status = dae_object.GetCalculationState();

		if (status > 0)
		{
			std::string message("BzzDae System successfully solved: ");

			if (status == 1)		message += "INITIALIZATION_STATE";
			else if (status == 2)	message += "CONTINUATION_STATE";
			else if (status == 10)	message += "INTEGRATION_STOPPED_BEFORE_RECALCULATING_JACOBIAN";
			else if (status == 11)	message += "INTEGRATION_STOPPED_WHEN_SUM_ABS_Y1_IS_LESS_THAN";

			std::cout << message << std::endl;

			BzzVector yp = dae_object.GetY1InMeshPoint();
			const double sum_yp = yp.GetSumAbsElements();

			std::cout << std::endl;
			std::cout << " * CPU time (s):                   " << timeEnd - timeStart << std::endl;
			std::cout << " * number of steps:                " << dae_object.GetNumStep() << std::endl;
			std::cout << " * number of functions:            " << dae_object.GetNumFunction() << std::endl;
			std::cout << " * number of solutions:            " << dae_object.GetNumSolution() << std::endl;
			std::cout << " * number of Jacobians:            " << dae_object.GetNumNumericalJacobian() << std::endl;
			std::cout << " * number of factorizations:       " << dae_object.GetNumFactorization() << std::endl;
			std::cout << " * number of functions (Jacobian): " << dae_object.GetNumFunctionForJacobian() << std::endl;
			std::cout << " * last order:                     " << dae_object.GetOrderUsed() << std::endl;
			std::cout << " * last step size:                 " << std::scientific << dae_object.GetHUsed() << std::endl;
			std::cout << " * mean y':                        " << std::scientific << sum_yp / double(neq) << std::endl;
			std::cout << std::endl;

			object->CorrectedUnknownsVector(dae_object.GetYInMeshPoint().GetHandle());
		}
		else
		{
			std::string message("BzzDae Solver Error: ");
			if (status == -1)		message += "EXCESSIVE_WORK_STATE";
			else if (status == -2)	message += "TOO_MUCH_ACCURACY_REQUESTED_STATE";
			else if (status == -3)	message += "ILLEGAL_VALUE_OF_TOUT_STATEE";
			else if (status == -4)	message += "REPEATED_ERROR_TEST_FAILURE_STATE";
			else if (status == -5)	message += "CONVERGENCE_TEST_FAILURE_STATE";
			else if (status == -6)	message += "H_EQUAL_HMIN_STATE";
			else if (status == -7)	message += "YOU_MUST_USE_TCRITIC_STATE";
			else if (status == -8)	message += "TCRITIC_IS_BEHIND_TOUT";
			else if (status == -9)	message += "ILLEGAL_CONSTRAINTS";
			else if (status == -10)	message += "EXCEPTION_HANDLING_STOP";
			else if (status == -11)	message += "DEINITIALIZE_STATE";
			else if (status == -12)	message += "YOU_CANNOT_OVERSHOOT_TCRITIC";

			std::cout << message << std::endl;

			object->CorrectedUnknownsVector(yInitial.GetHandle());
		}

		return 0;
	}
}

