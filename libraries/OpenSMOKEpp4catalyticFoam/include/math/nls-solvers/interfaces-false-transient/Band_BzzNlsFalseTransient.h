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

#ifndef OpenSMOKE_Solve_Band_BzzNlsFalseTransient_H
#define OpenSMOKE_Solve_Band_BzzNlsFalseTransient_H

#include "math/nls-solvers/NonLinearSystemSolver"
#include "FalseTransient_Utilities.h"

namespace NlsSMOKE
{
	template<typename Object, typename System>
	int loop_BzzNls(int& global_count, Object* object, System& nls_system, BzzNonLinearSystemSparseObject& solver, const NlsSMOKE::FalseTransientSolver_Parameters& parameters);

	template<typename Object, typename System>
	int Solve_Band_BzzNlsFalseTransient(Object* object, BzzNonLinearSystemSparseObject& nls_object, const NlsSMOKE::FalseTransientSolver_Parameters& parameters)
	{
		std::cout << "Band False Transient solution (BzzNls)..." << std::endl;

		const unsigned int neq = object->NumberOfEquations();

		// First guess solution
		BzzVector yInitial(neq);
		object->UnknownsVector(yInitial.GetHandle());

		// Differential vs Algebraic
		BzzVectorInt indices_differential_algebraic(neq);
		{
			BzzVector temp(neq);
			object->AlgebraicDifferentialVector(temp.GetHandle());
			for (int i = 1; i <= object->NumberOfEquations(); i++)
				indices_differential_algebraic[i] = int(temp[i]);
		}

		System nls_system;
		nls_system.assign(object);
		nls_system.SetDifferentialAlgebraic(indices_differential_algebraic);
		nls_system.SetInitialConditions(yInitial);
		nls_system.SetTimeStep(parameters.initial_step());
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

		// Maximum number of iterations
		if (parameters.maximum_number_iterations() > 0.)
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

		double timeStart = BzzGetCpuTime();

		// Call main solver
		int global_count = 0;
		int flag_solver;
		for (;;)
		{
			int loop_result = loop_BzzNls(global_count, object, nls_system, nls_object, parameters);

			// The loop was successfully accomplished
			if (loop_result == 0)
			{
				flag_solver = 0;
				break;
			}
			// The time step was too large
			else
			{
				nls_system.SetTimeStep(nls_system.deltat() / parameters.decrement_factor());
				nls_system.SetTimeStep(std::max(nls_system.deltat(), parameters.minimum_step()));

				std::cout << "Decrease time step..." << std::endl;
				std::cout << "dt=" << nls_system.deltat() << " alfa=" << parameters.decrement_factor() << " n=" << global_count << std::endl;

				flag_solver = -1;
			}
		}

		double timeEnd = BzzGetCpuTime();

		return flag_solver;
	}

	template<typename Object, typename System>
	int loop_BzzNls(int& global_count, Object* object, System& nls_system, BzzNonLinearSystemSparseObject& solver, const NlsSMOKE::FalseTransientSolver_Parameters& parameters)
	{
		unsigned int count = 0;

		do
		{
			// Update the Jacobian at first call: true
			if (count == 0)
			{
			}

			if (count == parameters.steps_before_increasing())
			{
				nls_system.SetTimeStep(nls_system.deltat() * parameters.increment_factor());
				nls_system.SetTimeStep(std::min(nls_system.deltat(), parameters.maximum_step()));
				count = 0;

				std::cout << "dt=" << nls_system.deltat() << " beta=" << parameters.increment_factor() << " n=" << global_count << std::endl;
			}

			if (count == 0)
			{
				//solver.Restart(nls_system.InitialConditions(), &nls_system, J, 0, object->BlockDimensions());
			}
			else
			{
				int status = solver();
			}

			BzzVector y(object->NumberOfEquations()), f(object->NumberOfEquations());
			solver.GetSolution(&y, &f);

			object->CorrectedUnknownsVector(y.GetHandle());
			object->norm();

			nls_system.SetInitialConditions(y);
			solver.Restart(y);

			count++;
			global_count++;

		} 
		while (global_count <= parameters.minimum_number_steps());

		return 0;
	}
}

#endif // OpenSMOKE_Solve_Band_BzzNlsFalseTransient_H


