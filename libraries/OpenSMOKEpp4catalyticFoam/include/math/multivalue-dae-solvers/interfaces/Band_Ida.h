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

#include <boost/timer/timer.hpp>
#include <ida/ida_spgmr.h>
#include <ida/ida_spbcgs.h>
#include <ida/ida_sptfqmr.h>

namespace DaeSMOKE
{
	static void IDAStatistics(void *mem, realtype cpuTime, N_Vector y, N_Vector yp);
	
	void CheckSolver(const int flag_solver);

	template<typename Object>
	int Solve_Band_Ida(Object* object, const DaeSMOKE::DaeSolver_Parameters& parameters, const double t0, const double tEnd)
	{
		std::cout << "Band DAE solution (IDA)..." << std::endl;

		const int neq = object->NumberOfEquations();

		void *mem;
		int flag;

		N_Vector y;
		N_Vector yp;
		N_Vector res;
		N_Vector id;
		N_Vector yInitial;
		IDAUserData data;

		mem = NULL;
		data = NULL;

		y = NULL;
		yp = NULL;
		res = NULL;
		id = NULL;
		yInitial = NULL;

		// Memory allocation
		{
			y = N_VNew_Serial(neq);
			if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
			yp = N_VNew_Serial(neq);
			if (check_flag((void *)yp, "N_VNew_Serial", 0)) return(1);
			res = N_VNew_Serial(neq);
			if (check_flag((void *)res, "N_VNew_Serial", 0)) return(1);
			id = N_VNew_Serial(neq);
			if (check_flag((void *)id, "N_VNew_Serial", 0)) return(1);
			yInitial = N_VNew_Serial(neq);
			if (check_flag((void *)yInitial, "N_VNew_Serial", 0)) return(1);
		}


		// User data
		data = (IDAUserData)malloc(sizeof *data);

		//Initialize
		{
			object->UnknownsVector(NV_DATA_S(y));
			object->UnknownsVector(NV_DATA_S(yInitial));
		}

		// Initialize and allocate memory
		mem = IDACreate();
		if (check_flag((void *)mem, "IDACreate", 0)) return(1);

		// Set user data
		flag = IDASetUserData(mem, data);
		if (check_flag(&flag, "IDASetUserData", 1)) return(1);

		if (parameters.sparse_linear_algebra() == true)
		{
			data->J = N_VNew_Serial(neq);
			if (check_flag((void *)data->J, "N_VNew_Serial", 0)) return(1);
			data->invJ = N_VNew_Serial(neq);
			if (check_flag((void *)data->invJ, "N_VNew_Serial", 0)) return(1);
		}

		// Differential/Algebraic equations
		{
			object->AlgebraicDifferentialVector(NV_DATA_S(id));

			flag = IDASetId(mem, id);
			if (check_flag(&flag, "IDASetId", 1)) return(1);
		}

		// Assign the system
		{
			//Initialize yp vector to 0.
			ida_initial_derivatives(t0, y, yp, data);
			flag = IDAInit(mem, ida_equations, ida_print_solution, t0, y, yp);
			if (check_flag(&flag, "IDAInit", 1)) return(1);
		}


		// Set optional input
		{
			// Tolerances
			if (parameters.relative_tolerances().size() == neq)
			{
				// Allocate memory for absolute tolerances vector
				N_Vector abstol;
				abstol = N_VNew_Serial(neq);
				if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

				// Fill the vector
				for (int i=0;i<neq;i++)
					NV_Ith_S(abstol,i) = parameters.absolute_tolerances()(i);

				// Assign the vector
				flag = IDASVtolerances(mem, parameters.relative_tolerance(), abstol);
				if (check_flag(&flag, "IDASVtolerances", 1)) return(1);

				// Free the memory
				N_VDestroy_Serial(abstol);
			}
			else
			{
				flag = IDASStolerances(mem, parameters.relative_tolerance(), parameters.absolute_tolerance());
				if (check_flag(&flag, "IDASStolerances", 1)) return(1);
			}

			// Maximum number of steps
			flag = IDASetMaxNumSteps(mem, parameters.maximum_number_of_steps());
			if (check_flag(&flag, "IDASetMaxNumSteps", 1)) return(1);

			/*
			// Maximum number of error test failures
			flag = IDASetMaxErrTestFails(mem, parameters.maximum_err_test_fails());
			if (check_flag(&flag, "IDASetMaxErrTestFails", 1)) return(1);

			// Maximum number of non linear iterations
			flag = IDASetMaxNonlinIters(mem, parameters.maximum_nl_iter());
			if (check_flag(&flag, "IDASetMaxNonlinIters", 1)) return(1);

			// Maximum number of convergence failures
			flag = IDASetMaxConvFails(mem, parameters.maximum_conv_fails());
			if (check_flag(&flag, "IDASetMaxConvFails", 1)) return(1);

			// Maximum number of non linear convergence tests
			flag = IDASetNonlinConvCoef(mem, parameters.coefficient_nonlinear_convergence_test());
			if (check_flag(&flag, "IDASetNonlinConvCoef", 1)) return(1);

			// Maximum order
			flag = IDASetMaxOrd(mem, parameters.maximum_order());
			if (check_flag(&flag, "IDASetMaxOrd", 1)) return(1);

			// Initial step size
			if (parameters.initial_step() > 0.)
			{
			flag = IDASetInitStep(mem, parameters.initial_step());
			if (check_flag(&flag, "IDASetInitStep", 1)) return(1);
			}

			// Maximum step size
			if (parameters.maximum_step() > 0.)
			{
			flag = IDASetMaxStep(mem, parameters.maximum_step());
			if (check_flag(&flag, "IDASetInitStep", 1)) return(1);
			}

			*/
		}

		// Set minimum constraints
		if (parameters.non_negative_unknowns() == true)
		{
			N_Vector constraints;
			constraints = NULL;
			constraints = N_VNew_Serial(neq);
			if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);

			// All the unknowns will be constrained to be >= 0
			N_VConst(1., constraints);
			flag = KINSetConstraints(mem, constraints);
			if (check_flag(&flag, "KINSetConstraints", 1)) return(1);

			N_VDestroy_Serial(constraints);
		}

		// Set maximum constraints
		if (parameters.maximum_constraints() == true)
		{
			// Sundials IDA does not support this function
		}

		// Attach band linear solver
		{
			// Band solver
			if (parameters.sparse_linear_algebra() == false)
			{
				int mu = object->UpperBand();
				int ml = object->LowerBand();

				#if OPENSMOKE_USE_MKL == 1
				{
					std::cout << "IDA with Lapack solver..." << std::endl;
					flag = IDALapackBand(mem, neq, mu, ml);
					if (check_flag(&flag, "IDALapackBand", 1)) return(1);
				}
				#else
				{
					flag = IDABand(mem, neq, mu, ml);
					if (check_flag(&flag, "IDABand", 1)) return(1);
				}
				#endif
			}

			// Iterative (sparse) solver
			else
			{
				// Call IDASpgmr to specify the linear solver
				flag = IDASpgmr(mem, 0);
				if (check_flag(&flag, "IDASpgmr", 1)) return(1);

				// Specify preconditioner
				flag = IDASpilsSetPreconditioner(mem, ida_preconditioner_setup, ida_preconditioner_solution);
				if (check_flag(&flag, "IDASpilsSetPreconditioner", 1)) return(1);
			}
		}

		double increasing_factor = 10.;
		unsigned int n_steps = 7;
		double t_initial = t0;
		double dt0 = 0.;
		
		if (t0 == 0.)
			t_initial = tEnd / std::pow(increasing_factor, double(n_steps-1));
		else
		{
			double sum = 1.;
			for (unsigned int i = 1; i <= n_steps - 1; i++)
				sum += std::pow(increasing_factor, double(i));
			dt0 = (tEnd - t0) / sum;
			t_initial = t0 + dt0;
		}
		

		// Call IDACalcIC to correct the initial values.
		std::cout << "Initial conditions..." << t0 << " " << tEnd << std::endl;
		flag = IDACalcIC(mem, IDA_YA_YDP_INIT, t_initial);
		if (check_flag(&flag, "IDACalcIC", 1)) return(1);

		// Loop over output times, call IDASolve, and print results.	
		int flag_solver;
		boost::timer::cpu_timer timer;
		double tout = t_initial;
		for (unsigned int iout = 1; iout <= n_steps; iout++)
		{
			std::cout << "Integrating to: " << tout << std::endl;

			realtype tret;

			flag_solver = IDASolve(mem, tout, &tret, y, yp, IDA_NORMAL);
			if (check_flag(&flag_solver, "IDASolve", 1)) return(1);
			CheckSolver(flag_solver);

			if (flag_solver >= 0)
			{
				realtype sum_abs_yp = N_SumAbs(yp);
				realtype norm2_yp = N_Norm2(yp);

				// Get scaled norm of the system function
				std::cout << " IDA solution: " << "||y'||sum = " << sum_abs_yp << " ||y'||2 = " << norm2_yp << std::endl;

				// Check if the solution reached the steady state conditions
				if (sum_abs_yp < parameters.minimum_yp()*neq)
					break;
			}
			else
			{
				break;
			}

			if (t0 == 0.)	tout *= increasing_factor;
			else            tout += dt0*std::pow(increasing_factor, double(iout));
		}
		boost::timer::cpu_times elapsed = timer.elapsed();

		if (flag_solver >= 0)
		{
			// Give results back
			object->CorrectedUnknownsVector(NV_DATA_S(y));

			// Print remaining counters and free memory.
			double cpuTime = elapsed.wall / 1.e9;
			IDAStatistics(mem, cpuTime, y, yp);
		}
		else
		{
			// Recover initial solution
			object->CorrectedUnknownsVector(NV_DATA_S(yInitial));
		}

		// Free memory
		IDAFree(&mem);
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(yInitial);
		N_VDestroy_Serial(yp);
		N_VDestroy_Serial(id);
		N_VDestroy_Serial(res);
		free(data);

		return flag_solver;
	}


	static void IDAStatistics(void *mem, realtype cpuTime, N_Vector y, N_Vector yp)
	{
		int flag;
		realtype hused;
		long int nst, nni, nje, nre, nreLS;
		int kused;

		int neq = NV_LENGTH_S(y);

		check_flag(&flag, "IDAGetNumSteps", 1);
		flag = IDAGetNumNonlinSolvIters(mem, &nni);
		flag = IDAGetLastOrder(mem, &kused);
		check_flag(&flag, "IDAGetLastOrder", 1);
		flag = IDAGetLastStep(mem, &hused);
		check_flag(&flag, "IDAGetLastStep", 1);
		flag = IDADlsGetNumJacEvals(mem, &nje);
		check_flag(&flag, "IDADlsGetNumJacEvals", 1);
		flag = IDAGetNumResEvals(mem, &nre);
		check_flag(&flag, "IDAGetNumResEvals", 1);
		flag = IDADlsGetNumResEvals(mem, &nreLS);
		check_flag(&flag, "IDADlsGetNumResEvals", 1);
		flag = IDAGetNumSteps(mem, &nst);
		check_flag(&flag, "IDAGetNumNonlinSolvIters", 1);


		std::cout << std::endl;
		std::cout << " * CPU time (s):                   " << cpuTime << std::endl;
		std::cout << " * number of steps:                " << nni << std::endl;
		std::cout << " * number of functions:            " << nre << std::endl;
		std::cout << " * number of non linear iter.:     " << nst << std::endl;
		std::cout << " * number of Jacobians:            " << nje << std::endl;
		std::cout << " * dummy:                          " << 0 << std::endl;
		std::cout << " * number of functions (Jacobian): " << nreLS << std::endl;
		std::cout << " * last order:                     " << kused << std::endl;
		std::cout << " * last step size:                 " << std::scientific << hused << std::endl;
		std::cout << " * mean y':                        " << std::scientific << N_SumAbs(yp) / double(neq) << std::endl;
		std::cout << std::endl;
	}

	void CheckSolver(const int flag_solver)
	{
		if (flag_solver >= 0)
		{
			std::string message("IDA successfully solved: ");

			// IDASolve succeeded
			if (flag_solver == IDA_SUCCESS)				message += "IDA_SUCCESS";
			// IDASolve succeeded by reaching the stop point specified through the optional input function IDASetStopTime.
			else if (flag_solver == IDA_TSTOP_RETURN)	message += "IDA_TSTOP_RETURN";
			// IDASolve succeeded and found one or more roots.In this case,
			// tret is the location of the root.If nrtfn > 1, call IDAGetRootInfo
			// to see which gi were found to have a root.See §4.5.9.3 for more information.
			else if (flag_solver == IDA_ROOT_RETURN)	message += "IDA_ROOT_RETURN";
		}
		else
		{
			std::string message("IDA solver error: ");

			// The ida mem argument was NULL
			if (flag_solver == IDA_MEM_NULL)			message += "IDA_MEM_NULL";
			// One of the inputs to IDASolve was illegal, or some other input to the solver was either illegal or missing.The latter category
			// includes the following situations : (a)The tolerances have not been set. (b) A component of the error weight vector became zero during
			// internal time - stepping. (c)The linear solver initialization function (called by the user after calling IDACreate) failed to set the linear
			// solver - specific lsolve field in ida mem. (d)A root of one of the root functions was found both at a point t and also very near t.In
			// any case, the user should see the printed error message for details.
			else if (flag_solver == IDA_ILL_INPUT)		message += "IDA_ILL_INPUT";
			// The solver took mxstep internal steps but could not reach tout. 
			// The default value for mxstep is MXSTEP DEFAULT = 500.
			else if (flag_solver == IDA_TOO_MUCH_WORK)	message += "IDA_TOO_MUCH_WORK";
			// The solver could not satisfy the accuracy demanded by the user for some internal step.
			else if (flag_solver == IDA_TOO_MUCH_ACC)	message += "IDA_TOO_MUCH_ACC";
			//  Error test failures occurred too many times(MXNEF = 10) during
			// one internal time step or occurred with | h | = hmin.
			else if (flag_solver == IDA_ERR_FAIL)		message += "IDA_ERR_FAIL";
			//  Convergence test failures occurred too many times(MXNCF = 10)
			// during one internal time step or occurred with | h | = hmin.
			else if (flag_solver == IDA_CONV_FAIL)		message += "IDA_CONV_FAIL";
			// The linear solver’s initialization function failed.
			else if (flag_solver == IDA_LINIT_FAIL)		message += "IDA_LINIT_FAIL";
			// The linear solver’s setup function failed in an unrecoverable manner.
			else if (flag_solver == IDA_LSETUP_FAIL)		message += "IDA_LSETUP_FAIL";
			// The linear solver’s solve function failed in an unrecoverable manner.
			else if (flag_solver == IDA_LSOLVE_FAIL)		message += "IDA_LSOLVE_FAIL";
			// The inequality constraints were violated and the solver was unable to recover.
			else if (flag_solver == IDA_CONSTR_FAIL)		message += "IDA_CONSTR_FAIL";
			//  The user’s residual function repeatedly returned a recoverable error flag, but the solver was unable to recover.
			else if (flag_solver == IDA_REP_RES_ERR)		message += "IDA_REP_RES_ERR";
			// The user’s residual function returned a nonrecoverable error flag.
			else if (flag_solver == IDA_RES_FAIL)		message += "IDA_RES_FAIL";
			//  The rootfinding function failed.
			else if (flag_solver == IDA_RTFUNC_FAIL)		message += "IDA_RTFUNC_FAIL";

			std::cout << message << std::endl;
		}
	}
}