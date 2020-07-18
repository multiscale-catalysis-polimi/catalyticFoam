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

#ifndef OpenSMOKE_Solve_Band_BzzNls_H
#define OpenSMOKE_Solve_Band_BzzNls_H

#include "math/nls-solvers/NonLinearSystemSolver"

namespace NlsSMOKE
{
	template<typename Object, typename System>
	int Solve_Band_BzzNls(Object* object, BzzNonLinearSystemSparseObject& nls_object, const NlsSMOKE::NonLinearSolver_Parameters& parameters)
	{
		std::cout << "Band NLS solution (BzzNls)..." << std::endl;

		const unsigned int neq = object->NumberOfEquations();

		// First guess solution
		BzzVector yInitial(neq);
		object->UnknownsVector(yInitial.GetHandle());

		System nls_system;
		nls_system.assign(object);
		nls_object(yInitial, &nls_system, object->LowerBand(), object->UpperBand());

		// Minimum Constraints
		if (parameters.minimum_constraints() == true)
		{
			BzzVector xMin(neq);
			object->MinimumUnknownsVector(xMin.GetHandle());
			nls_object.SetMinimumConstraints(xMin);
		}

		// Maximum constraints
		if (parameters.maximum_constraints() == true)
		{
			BzzVector xMax(neq);
			object->MaximumUnknownsVector(xMax.GetHandle());
			nls_object.SetMaximumConstraints(xMax);
		}

		// Maximum step
		if (parameters.maximum_number_iterations() > 0)
			nls_object.SetMaxNewtonCallsNumber(parameters.maximum_number_iterations());

		// Tolerances
		if (parameters.tolerance_absolute() > 0. && parameters.tolerance_relative() > 0.)
			nls_object.SetTolerance(parameters.tolerance_absolute(), parameters.tolerance_relative());
		else if (parameters.tolerance_absolute() > 0. && parameters.tolerance_relative() <= 0.)
			nls_object.SetTolerance(parameters.tolerance_absolute(), std::sqrt(MachEpsFloat()));
		else if (parameters.tolerance_absolute() <= 0. && parameters.tolerance_relative() > 0.)
			nls_object.SetTolerance(1.e-10, parameters.tolerance_relative());

		// Scaling factors (ONE means no scaling)
		if (parameters.scaling_policy() != 0)
		{
			BzzVector w(neq);
			for (unsigned int i = 1; i <= neq; i++)
				w[i] = 1. / std::max(1.e-2, yInitial[i]);
			nls_object.SetWeights(w);
		}

		// Verbose
		// if (parameters.verbose() == true)
		//	 nls_object.StepPrint(NlsPrint);

		// Solving the system
		double timeStart = BzzGetCpuTime();
		const int status = nls_object();
		double timeEnd = BzzGetCpuTime();

		// Analyze the result
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

			BzzVector y(neq), f(neq);
			nls_object.GetSolution(&y, &f);

			std::cout << std::endl;
			std::cout << " * CPU time (s):                   " << timeEnd - timeStart << std::endl;
			std::cout << " * number of iterations:           " << nls_object.IterationCounter() << std::endl;
			std::cout << " * number of functions (gradient): " << nls_object.NumFunctionsForGradient() << std::endl;
			std::cout << " * number of Newtons:              " << nls_object.NumNewtons() << std::endl;
			std::cout << " * number of Jacobians:            " << nls_object.NumNumericalJacobians() << std::endl;
			std::cout << " * number of functions (Jacobian): " << nls_object.NumFunctionsForNumericalJacobian() << std::endl;
			std::cout << std::endl;

			object->CorrectedUnknownsVector(y.GetHandle());
		}
		else
		{
			std::string message("Nls Solver Error: ");
			if (status == -1)		message += "It has been impossible to reach the solution";
			else if (status == -2)	message += "The search has been stopped";
			else if (status == -3)	message += "The object is not initialized";
			else if (status == -4)	message += "It has been impossible to reach the solution in Restart";

			std::cout << message << std::endl;

			object->CorrectedUnknownsVector(yInitial.GetHandle());
		}

		return status;
	}
}

#endif // OpenSMOKE_Solve_Band_BzzNls_H


