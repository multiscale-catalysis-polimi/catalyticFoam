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

#ifndef OpenSMOKE_Solve_Sparse_OpenSMOKEppFalseTransient_H
#define OpenSMOKE_Solve_Sparse_OpenSMOKEppFalseTransient_H

#include "math/nls-solvers/NonLinearSystemSolver"
#include "FalseTransient_Utilities.h"

namespace NlsSMOKE
{
	template<typename Object, typename System>
	int Solve_Sparse_OpenSMOKEppFalseTransient(Object* object, const NlsSMOKE::FalseTransientSolver_Parameters& parameters)
	{
		std::cout << "Sparse False Transient solution (OpenSMOKE++)..." << std::endl;

		// Total number of equations
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

		typedef NlsSMOKE::KernelSparse<System> kernel;
		NlsSMOKE::NonLinearSolver<kernel> nls_solver;

		// Set initial conditions
		nls_solver.assign(object);
		nls_solver.SetSparsityPattern(rows, cols, false);
		nls_solver.SetLinearAlgebraSolver(parameters.jacobian_solver());
		nls_solver.SetPreconditioner(parameters.preconditioner());
		nls_solver.SetDifferentialAlgebraic(equation_type);
		nls_solver.SetTimeStep(parameters.initial_step());
		nls_solver.SetInitialConditions(yInitial);
		nls_solver.SetFirstGuessSolution(yInitial);

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

		// No verbose output
		nls_solver.SetPrint(false);

		// Scaling factors (ONE means no scaling)
		if (parameters.scaling_policy() != 0)
		{
			Eigen::VectorXd w(neq);
			for (unsigned int i = 0; i < neq; i++)
				w(i) = 1. / std::max(1.e-2, yInitial(i));
			nls_solver.SetWeights(w);
		}

		// Call main solver
		int global_count = 0;
		int flag_solver = 0;

		double timeStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		for (;;)
		{
			const int loop_result = loop_OpenSMOKEpp(global_count, object, nls_solver, parameters);

			// The loop was successfully accomplished
			if (loop_result == 0)
			{
				flag_solver = 0;
				break;
			}
			// The time step was too large
			else
			{
				nls_solver.SetTimeStep(nls_solver.deltat() / parameters.decrement_factor());
				nls_solver.SetTimeStep(std::max(nls_solver.deltat(), parameters.minimum_step()));

				std::cout << "Decrease time step..." << std::endl;
				std::cout << "dt=" << nls_solver.deltat() << " alfa=" << parameters.decrement_factor() << " n=" << global_count << std::endl;

				flag_solver = -1;
			}
		}

		double timeEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		return flag_solver;
	}
}

#endif // OpenSMOKE_Solve_Sparse_OpenSMOKEppFalseTransient_H


