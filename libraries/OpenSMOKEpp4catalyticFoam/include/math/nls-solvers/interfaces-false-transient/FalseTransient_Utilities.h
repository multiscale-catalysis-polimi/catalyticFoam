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

#ifndef OpenSMOKE_FalseTransient_Utilities_H
#define OpenSMOKE_FalseTransient_Utilities_H

#include "math/nls-solvers/NonLinearSystemSolver"

namespace NlsSMOKE
{
	template<typename Object, typename Solver>
	int loop_OpenSMOKEpp(int& global_count, Object* object, Solver& nls_solver, const NlsSMOKE::FalseTransientSolver_Parameters& parameters)
	{
		unsigned int count = 0;
		do
		{
			// The time step is kept constant for a specific number of times
			if (count == parameters.steps_before_increasing())
			{
				nls_solver.SetTimeStep(nls_solver.deltat() * parameters.increment_factor());
				nls_solver.SetTimeStep(std::min(nls_solver.deltat(), parameters.maximum_step()));
				count = 0;

				std::cout << "dt=" << nls_solver.deltat() << " beta=" << parameters.increment_factor() << " n=" << global_count << std::endl;
			}

			// Solve the system of equations
			if (count > 0)
				nls_solver();

			// Extract the solution
			const unsigned int neq = object->NumberOfEquations();
			Eigen::VectorXd y(neq), f(neq);
			nls_solver.Solution(y, f);

			// Move the solution to the flame
			object->CorrectedUnknownsVector(y.data());
			object->norm();

			// Prepare the solver with a new set of initial conditions (new step)
			nls_solver.SetInitialConditions(y);
			nls_solver.Reset(y);

			// Update counters
			count++;
			global_count++;

		} while (global_count <= parameters.minimum_number_steps());

		return 0;
	}
}

#endif // OpenSMOKE_FalseTransient_Utilities_H


