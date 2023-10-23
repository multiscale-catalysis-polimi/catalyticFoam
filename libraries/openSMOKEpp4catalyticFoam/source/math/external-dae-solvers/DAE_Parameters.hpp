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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

namespace OpenSMOKE
{
	DAE_Parameters::DAE_Parameters()
	{
		// Default values
		type_ = DAE_INTEGRATOR_OPENSMOKE;
		linear_algebra_ = "Eigen";
		absolute_tolerance_ = 1.e-10;
		relative_tolerance_ = 100.*MachEpsFloat();
		number_of_steps_ = 0;
		minimum_step_ = -1.;
		maximum_step_ = -1.;
		maximum_number_of_steps_ = 5000000;
		maximum_err_test_fails_ = 50;
		maximum_conv_fails_ = 50;
		maximum_nl_iter_ = 50;
		initial_step_ = -1.;
		maximum_order_ = -1;
		full_pivoting_ = false;
		
		// Reset the counters
		time_spent_to_factorize_ = 0.;
		time_spent_to_evaluate_jacobian_ = 0.;
		time_spent_to_solve_linear_system_ = 0.;
		cpu_time_ = 0.;
		number_of_function_calls_ = 0;
		number_of_function_calls_jacobian_ = 0;
		number_of_jacobians_ = 0;
		number_of_factorizations_ = 0;
		number_of_steps_ = 0;
		last_order_used_ = 0;
		last_step_used_ = 0.;
		max_order_used_ = 0;
		max_step_used_ = 0.;
		min_step_used_ = 0.;
		number_of_nonlinear_iterations_ = 0;
		number_of_convergence_failures_ = 0;
		number_of_error_test_failures_ = 0;

		// Names of DAE solvers
		list_names_.push_back("OpenSMOKE");
		list_names_.push_back("IDA");
		list_names_.push_back("DASPK");
		list_names_.push_back("BzzDae");
	}

	void DAE_Parameters::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_DAEParameters_Options grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@DaeSolver") == true)
		{
			std::string name;
			dictionary.ReadString("@DaeSolver", name);
			
			if (name == "OpenSMOKE") 
			{
				type_ = DAE_INTEGRATOR_OPENSMOKE;
			}
			else if (name == "BzzDae")
			{
				#if OPENSMOKE_USE_BZZMATH == 0
					FatalErrorMessage("OpenSMOKE++ was built without the BzzMath support. Please select a different DAE solver");
				#endif
				type_ = DAE_INTEGRATOR_BZZDAE;
			}
			else if (name == "IDA")
			{
				#if OPENSMOKE_USE_SUNDIALS == 0
					FatalErrorMessage("OpenSMOKE++ was built without the Sundials (IDA) support. Please select a different DAE solver");
				#endif
				type_ = DAE_INTEGRATOR_IDA;
			}
			else if (name == "DASPK")
			{
				#if OPENSMOKE_USE_DASPK == 0
					FatalErrorMessage("OpenSMOKE++ was built without the DASPK support. Please select a different DAE solver");
				#endif
				type_ = DAE_INTEGRATOR_DASPK;
			}

			else FatalErrorMessage("Unknown DAE Solver: " + name);
		}

		if (dictionary.CheckOption("@LinearAlgebra") == true)
		{
			std::string name;
			dictionary.ReadString("@LinearAlgebra", name);
				 if (name == "Eigen")    linear_algebra_ = "Eigen";
			else if (name == "BzzMath")  linear_algebra_ = "BzzMath";
			else if (name == "Plasma")   linear_algebra_ = "Plasma";
			else if (name == "Flame")    linear_algebra_ = "Flame";
			else FatalErrorMessage("Unknown Linear Algebra: " + name);
		}

		if (dictionary.CheckOption("@RelativeTolerance") == true)
			dictionary.ReadDouble("@RelativeTolerance", relative_tolerance_);

		if (dictionary.CheckOption("@AbsoluteTolerance") == true)
			dictionary.ReadDouble("@AbsoluteTolerance", absolute_tolerance_);

		if (dictionary.CheckOption("@MaximumOrder") == true)
			dictionary.ReadInt("@MaximumOrder", maximum_order_);

		if (dictionary.CheckOption("@MaximumStep") == true)
			dictionary.ReadDouble("@MaximumStep", maximum_step_);

		if (dictionary.CheckOption("@MinimumStep") == true)
			dictionary.ReadDouble("@MinimumStep", minimum_step_);

		if (dictionary.CheckOption("@InitialStep") == true)
			dictionary.ReadDouble("@InitialStep", initial_step_);

		if (dictionary.CheckOption("@MaximumNumberOfSteps") == true)
			dictionary.ReadInt("@MaximumNumberOfSteps", maximum_number_of_steps_);

		if (dictionary.CheckOption("@FullPivoting") == true)
			dictionary.ReadBool("@FullPivoting", full_pivoting_);
	}

	void DAE_Parameters::Status(std::ostream& fOut)
	{
		double other_cpu_time = cpu_time_ - (time_spent_to_evaluate_jacobian_ + time_spent_to_factorize_ + time_spent_to_solve_linear_system_);

		fOut << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << "   DAE Integration                                                           " << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
		fOut << " * Integrator:                          " << list_names_[type_] << std::endl;
		
		fOut << " * CPU Time (s):                        " << std::scientific << cpu_time_ << std::endl;
		fOut << " *   assembling Jacobian:               " << std::scientific << time_spent_to_evaluate_jacobian_;
		fOut                                               << std::fixed << " (" << time_spent_to_evaluate_jacobian_ / cpu_time_*100. << "%)" << std::endl;
		fOut << " *   LU decomposition:                  " << std::scientific << time_spent_to_factorize_;
		fOut                                               << std::fixed << " (" << time_spent_to_factorize_ / cpu_time_*100. << "%)" << std::endl;
		fOut << " *   linear system solution:            " << std::scientific << time_spent_to_solve_linear_system_;
		fOut                                               << std::fixed << " (" << time_spent_to_solve_linear_system_ / cpu_time_*100. << "%)" << std::endl;
		fOut << " *   other:                             " << std::scientific << other_cpu_time;
		fOut                                               << std::fixed << " (" << other_cpu_time / cpu_time_*100. << "%)" << std::endl;
		fOut << " * Relative tolerance:                  " << std::scientific << relative_tolerance_ << std::endl;
		fOut << " * Absolute tolerance:                  " << std::scientific << absolute_tolerance_ << std::endl;
		fOut << " * Number of steps:                     " << std::fixed << number_of_steps_ << std::endl;
		fOut << " * Number of function calls:            " << std::fixed << number_of_function_calls_ + number_of_function_calls_jacobian_ << std::endl;
		fOut << " * Number of function calls (Jacobian): " << std::fixed << number_of_function_calls_jacobian_ << " (" << double(number_of_function_calls_jacobian_) / double(number_of_function_calls_ + number_of_function_calls_jacobian_)*100. << "%)" << std::endl;
		fOut << " * Number of Jacobians:                 " << std::fixed << number_of_jacobians_ << std::endl;
		fOut << " * Number of factorizations:            " << std::fixed << number_of_factorizations_ << std::endl;
		fOut << " * Number of convergence failures:      " << std::fixed << number_of_convergence_failures_ << std::endl;
		fOut << " * Number of error test failures:       " << std::fixed << number_of_error_test_failures_ << std::endl;
		fOut << " * Last order used:                     " << std::fixed << last_order_used_ << std::endl;
		fOut << " * Last step used:                      " << std::scientific << last_step_used_ << std::endl;
		fOut << " * Max order used:                      " << std::fixed << max_order_used_ << std::endl;
		fOut << " * Max step used:                       " << std::scientific << max_step_used_ << std::endl;
		fOut << "-----------------------------------------------------------------------------" << std::endl;
	}

	template <typename DaeSolver>
	void DAE_Parameters::TransferDataFromDaeSolver(const DaeSolver& dae_solver, const double cpuTime)
	{
		SetCPUTime(cpuTime);
		SetNumberOfFunctionCalls(dae_solver.numberOfFunctionCalls());
		SetNumberOfJacobians(dae_solver.numberOfJacobianEvaluations());
		SetNumberOfFactorizations(dae_solver.numberOfMatrixFactorizations());
		SetNumberOfSteps(dae_solver.numberOfSteps());
		SetLastOrderUsed(dae_solver.lastOrderUsed());
		SetLastStepUsed(dae_solver.lastStepUsed());
		SetMaxOrderUsed(dae_solver.maxOrderUsed());
		SetMinStepUsed(dae_solver.minimumStepUsed());
		SetMaxStepUsed(dae_solver.maximumStepUsed());
		SetNumberOfFunctionCallsForJacobian(dae_solver.numberOfFunctionCallsForJacobian());
		SetTimeSpentToFactorize(dae_solver.cpuTimeToFactorize());
		SetTimeSpentToEvaluateJacobian(dae_solver.cpuTimeToAssembleJacobian());
		SetTimeSpentToSolveLinearSystem(dae_solver.cpuTimeToSolveLinearSystem());
		SetNumberOfConvergenceFailures(dae_solver.numberOfConvergenceFailure());
		SetNumberOfErrorTestFailures(dae_solver.numberOfErrorCheckFailure());
	}
}
