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

#ifndef OpenSMOKE_Solve_Band_OpenSMOKEppNls_H
#define OpenSMOKE_Solve_Band_OpenSMOKEppNls_H

#include "math/nls-solvers/NonLinearSystemSolver"

namespace NlsSMOKE
{
	template<typename Object, typename System>
	int Solve_Band_OpenSMOKEppNls(Object* object, const NlsSMOKE::NonLinearSolver_Parameters& parameters)
	{
		std::cout << "Band NLS solution (OpenSMOKE++)..." << std::endl;

		const unsigned int neq = object->NumberOfEquations();

		// Initial conditions
		Eigen::VectorXd yInitial(neq);
		object->UnknownsVector(yInitial.data());

		typedef NlsSMOKE::KernelBand<System> kernel;
		NlsSMOKE::NonLinearSolver<kernel> nls_solver;

		// Set initial conditions
		nls_solver.assign(object);
		nls_solver.SetFirstGuessSolution(yInitial);
		nls_solver.SetBandSizes(object->UpperBand(), object->LowerBand());

		// Minimum Constraints
		if (parameters.minimum_constraints() == true)
		{
			Eigen::VectorXd xMin(neq);
			object->MinimumUnknownsVector(xMin.data());
			nls_solver.SetMinimumValues(xMin);
		}

		// Maximum constraints
		if (parameters.maximum_constraints() == true)
		{
			Eigen::VectorXd xMax(neq);
			object->MaximumUnknownsVector(xMax.data());
			nls_solver.SetMaximumValues(xMax);
		}

		// Maximum number of iterations
		if (parameters.maximum_number_iterations() > 0.)
			nls_solver.SetMaximumNumberOfIterations(parameters.maximum_number_iterations());

		// Relative tolerance
		if (parameters.tolerance_relative() > 0.)
			nls_solver.SetRelativeTolerances(parameters.tolerance_relative());

		// Absolute tolerance
		if (parameters.tolerance_absolute() > 0.)
			nls_solver.SetAbsoluteTolerances(parameters.tolerance_absolute());

		// Verbose
		//if (parameters.verbosity_level() > 0)
			nls_solver.SetPrint(true);
		//else
		//	nls_solver.SetPrint(false);

		// Scaling factors (ONE means no scaling)
		if (parameters.scaling_policy() != 0)
		{
			Eigen::VectorXd w(neq);
			for (unsigned int i = 0; i < neq; i++)
				w(i) = 1. / std::max(1.e-2, yInitial(i));
			nls_solver.SetWeights(w);
		}

		// Solve the non linear system
		double timeStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		const int status = nls_solver();
		double timeEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		// Analyze the solution
		if (status >= 0)
		{
			std::string message("Nls System successfully solved: ");
			if (status == 0)		message += "Start conditions";
			else if (status == 1)	message += "The maximum number of functions calls has been performed";
			else if (status == 2)	message += "The Newton correction has reached the required precision";
			else if (status == 3)	message += "The Quasi Newton correction has reached the required precision";
			else if (status == 4)	message += "The Gradient has reached the required precision";
			else if (status == 5)	message += "The Objective function phiNew has reached the required precision";
			else if (status == 6)	message += "The Objective function phiW has reached the required precision";
			else if (status == 7)	message += "The Objective function phiNew has reached the required precision but solution is dubious";
			else if (status == 8)	message += "Reached the assigned max value for Newton calls";

			std::cout << message << std::endl;

			Eigen::VectorXd y(neq), f(neq);
			nls_solver.Solution(y, f);

			//if (parameters.verbosity_level() > 0)
			nls_solver.NlsSummary(std::cout);

			object->CorrectedUnknownsVector(y.data());
		}
		else
		{
			std::string message("Nls Solver Error: ");
			if (status == -1)		message += "It has been impossible to reach the solution";
			else if (status == -2)	message += "The search has been stopped";
			else if (status == -3)	message += "The object is not initialized";
			else if (status == -4)	message += "It has been impossible to reach the solution in Restart";

			std::cout << message << std::endl;

			object->CorrectedUnknownsVector(yInitial.data());
		}

		return status;
	}
}

#endif // OpenSMOKE_Solve_Band_OpenSMOKEppNls_H


