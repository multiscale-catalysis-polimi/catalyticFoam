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

#ifndef OpenSMOKE_Solve_Band_KinSol_H
#define OpenSMOKE_Solve_Band_KinSol_H

#include "math/nls-solvers/NonLinearSystemSolver"

namespace NlsSMOKE
{
	static void KinSolStatistics(void *kmem);

	template<typename Object>
	int Solve_Band_KinSol(Object* object, const NlsSMOKE::NonLinearSolver_Parameters& parameters)
	{
		std::cout << "Band NLS solution (KinSol)..." << std::endl;

		const int neq = object->NumberOfEquations();

		N_Vector y;
		N_Vector yInitial;
		N_Vector scale;
		KINSolUserData data;
	
		int flag;
		void *kmem;

		y = NULL;
		yInitial = NULL;
		scale = NULL;
		kmem = NULL;
		data = NULL;

		// Memory allocation
		{
			y = N_VNew_Serial(neq);
			if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

			yInitial = N_VNew_Serial(neq);
			if (check_flag((void *)yInitial, "N_VNew_Serial", 0)) return(1);

			scale = N_VNew_Serial(neq);
			if (check_flag((void *)scale, "N_VNew_Serial", 0)) return(1);
		}

		// User data
		{
			data = (KINSolUserData)malloc(sizeof *data);
			data->yInitial = NV_DATA_S(yInitial);
		}

		// Initialize and allocate memory for KINSOL
		kmem = KINCreate();
		if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

		// Set user data
		flag = KINSetUserData(kmem, data);
		if (check_flag(&flag, "KINSetUserData", 1)) return(1);

		// Assign the system
		flag = KINInit(kmem, kinsol_equations, y);
		if (check_flag(&flag, "KINInit", 1)) return(1);

		// Set optional input
		{
			// Function tolerance
			if (parameters.tolerance_function() > 0.)
			{
				flag = KINSetFuncNormTol(kmem, parameters.tolerance_function());
				if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
			}

			// Step tolerance
			if (parameters.tolerance_step() > 0.)
			{
				flag = KINSetScaledStepTol(kmem, parameters.tolerance_step());
				if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);
			}

			// Relative error
			if (parameters.relative_error() > 0.)
			{
				flag = KINSetRelErrFunc(kmem, parameters.relative_error());
				if (check_flag(&flag, "KINSetRelErrFunc", 1)) return(1);
			}
		
			// Force a Jacobian re-evaluation every mset iterations (default 10)
			flag = KINSetMaxSetupCalls(kmem, parameters.maximum_setup_calls());
			if (check_flag(&flag, "KINSetMaxSetupCalls", 1)) return(1);

			// Every msubset iterations, test if a Jacobian evaluation is necessary (default 5)
			flag = KINSetMaxSubSetupCalls(kmem, parameters.maximum_sub_setup_calls());
			if (check_flag(&flag, "KINSetMaxSubSetupCalls", 1)) return(1);

			// Maximum number of nonlinear iterations allowed
			flag = KINSetNumMaxIters(kmem, parameters.maximum_number_iterations());
			if (check_flag(&flag, "KINSetNumMaxIters", 1)) return(1);

			// Maximum beta failures in linesearch algorithm
			flag = KINSetMaxBetaFails(kmem, parameters.maximum_beta_fails());
			if (check_flag(&flag, "KINSetMaxBetaFails", 1)) return(1);

			// Maximum allowed scaled lenght of Newton step
			if (parameters.maximum_newton_step() > 0.)
			{
				flag = KINSetMaxNewtonStep(kmem, parameters.maximum_newton_step());
				if (check_flag(&flag, "KINSetMaxNewtonStep", 1)) return(1);
			}

			// Parameters omega_min and omega_max
			if (parameters.omega_min() > 0. || parameters.omega_max() > 0.)
			{
				flag = KINSetResMonParams(kmem, parameters.omega_min(), parameters.omega_max());
				if (check_flag(&flag, "KINSetResMonParams", 1)) return(1);
			}

			// Eta choice
			{
				if (parameters.eta_choice() == 0)	flag = KINSetEtaForm(kmem, KIN_ETACHOICE1);
				if (parameters.eta_choice() == 1)	flag = KINSetEtaForm(kmem, KIN_ETACHOICE2);
				if (parameters.eta_choice() == 2)	flag = KINSetEtaForm(kmem, KIN_ETACONSTANT);
				if (check_flag(&flag, "KINSetEtaForm", 1)) return(1);
			}

			// Print level
			flag = KINSetPrintLevel(kmem, parameters.verbosity_level());
			if (check_flag(&flag, "KINSetPrintLevel", 1)) return(1);
		}

		// Set minimum constraints
		if (parameters.non_negative_unknowns() == true)
		{
			N_Vector constraints;
			constraints = NULL;
			constraints = N_VNew_Serial(neq);
			if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);

			N_VConst(1., constraints);
			flag = KINSetConstraints(kmem, constraints);
			if (check_flag(&flag, "KINSetConstraints", 1)) return(1);

			N_VDestroy_Serial(constraints);
		}

		// Set maximum constraints (TODO)
		if (parameters.maximum_constraints() == true)
		{
		}

		// Attach band linear solver
		{
			const int mu = object->UpperBand();
			const int ml = object->LowerBand();
		
			#if OPENSMOKE_USE_MKL == 1
			{
				flag = KINLapackBand(kmem, neq, mu, ml);
				if (check_flag(&flag, "KINLapackBand", 1)) return(1);
			}
			#else
			{
				flag = KINBand(kmem, neq, mu, ml);
				if (check_flag(&flag, "KINBand", 1)) return(1);
			}
			#endif
		}

		// Initialize
		object->UnknownsVector(NV_DATA_S(y));
		object->UnknownsVector(NV_DATA_S(yInitial));

		// Scaling factors (ONE means no scaling)
		N_VConst_Serial(1., scale);
		if (parameters.scaling_policy() != 0)
		{
			realtype *pt_scale = NV_DATA_S(scale);
			realtype *pt_y = NV_DATA_S(y);
			for (int i = 0; i < neq; i++)
				pt_scale[i] = 1. / std::max(1.e-2, pt_y[i]);
		}

		// Call main solver
		int flag_solver = 0;
		if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_BASIC)
		{
			flag_solver = KINSol(kmem, y, KIN_NONE, scale, scale);
			if (check_flag(&flag_solver, "KINSol::KIN_NONE", 0)) return(1);
		}
		else if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_NEWTON_GLOBALIZATION)
		{
			flag_solver = KINSol(kmem, y, KIN_LINESEARCH, scale, scale);
			if (check_flag(&flag_solver, "KINSol::KIN_LINESEARCH", 0)) return(1);
		}
		else if (parameters.strategy() == NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_FIXED_POINT)
		{
			flag_solver = KINSol(kmem, y, KIN_FP, scale, scale);
			if (check_flag(&flag_solver, "KINSol::KIN_FP", 0)) return(1);
		}
		else if (parameters.strategy() ==NlsSMOKE::NonLinearSolver_Parameters::NLS_STRATEGY_PICARD)
		{
			flag_solver = KINSol(kmem, y, KIN_PICARD, scale, scale);
			if (check_flag(&flag_solver, "KINSol::KIN_PICARD", 0)) return(1);
		}

		// Check solution
		if (flag_solver >= 0)
		{
			std::string message("NLS successfully solved: ");

			// KINSol succeeded; the scaled norm of F(u) is less than fnormtol
			if (flag_solver == KIN_SUCCESS)				message += "KIN_SUCCESS";
			// The guess u = u0 satisfied the system F(u) = 0 within the tolerances specified.
			else if (flag_solver == KIN_INITIAL_GUESS_OK)	message += "KIN_INITIAL_GUESS_OK";
			// kinsol stopped based on scaled step length. This means that the current iterate may be an approximate
			// solution of the given nonlinear system, but it is also quite possible that the algorithm is “stalled”
			// (making insufficient progress) near an invalid solution, or that the scalar scsteptol is too large
			// (see KINSetScaledStepTol in §4.5.4 to change scsteptol from its default value).
			else if (flag_solver == KIN_STEP_LT_STPTOL)	message += "KIN_STEP_LT_STPTOL";

			std::cout << message << std::endl;

			// Get scaled norm of the system function
			double fnorm;
			flag = KINGetFuncNorm(kmem, &fnorm);
			if (check_flag(&flag, "KINGetfuncNorm", 1)) return(1);

			std::cout << " * computed solution ||F|| = " << std::scientific << fnorm << std::endl;;
		
			if (parameters.verbosity_level() > 0)
				KinSolStatistics(kmem);

			// Move the solution to the object
			object->CorrectedUnknownsVector(NV_DATA_S(y));
		}
		else
		{
			std::string message("NLS solver error: ");

			// The kinsol memory block pointer was NULL.
			if (flag_solver == KIN_MEM_NULL)					message += "KIN_MEM_NULL";
			// An input parameter was invalid.
			else if (flag_solver == KIN_ILL_INPUT)				message += "KIN_ILL_INPUT";
			// The kinsol memory was not allocated by a call to KINCreate.
			else if (flag_solver == KIN_NO_MALLOC)				message += "KIN_NO_MALLOC";
			// The line search algorithm was unable to find an iterate sufficiently distinct from the
			// current iterate, or could not find an iterate satisfying the sufficient decrease condition.
			// Failure to satisfy the sufficient decrease condition could mean the current iterate
			// is “close” to an approximate solution of the given nonlinear system, the difference
			// approximation of the matrix - vector product J(u)v is inaccurate, or the real scalar
			// scsteptol is too large.
			else if (flag_solver == KIN_LINESEARCH_NONCONV)	message += "KIN_LINESEARCH_NONCONV";
			// The maximum number of nonlinear iterations has been reached.
			else if (flag_solver == KIN_MAXITER_REACHED)		message += "KIN_MAXITER_REACHED";
			// Five consecutive steps have been taken that satisfy the inequality kDupkL2 > 0.99
			// mxnewtstep, where p denotes the current step and mxnewtstep is a scalar upper
			// bound on the scaled step length.Such a failure may mean that kDFF(u)kL2 asymptotes
			// from above to a positive value, or the real scalar mxnewtstep is too small.
			else if (flag_solver == KIN_MXNEWT_5X_EXCEEDED)		message += "KIN_MXNEWT_5X_EXCEEDED";
			// The line search algorithm was unable to satisfy the “beta - condition” for MXNBCF + 1
			// nonlinear iterations(not necessarily consecutive), which may indicate the algorithm
			// is making poor progress.
			else if (flag_solver == KIN_LINESEARCH_BCFAIL)		message += "KIN_LINESEARCH_BCFAIL";
			// The user - supplied routine psolve encountered a recoverable error, but the preconditioner is already current.
			else if (flag_solver == KIN_LINSOLV_NO_RECOVERY)		message += "KIN_LINSOLV_NO_RECOVERY";
			// The linear solver initialization routine(linit) encountered an error.
			else if (flag_solver == KIN_LINIT_FAIL)		message += "KIN_LINIT_FAIL";
			// The user - supplied routine pset(used to set up the preconditioner data) encountered an unrecoverable error.
			else if (flag_solver == KIN_LSETUP_FAIL)		message += "KIN_LSETUP_FAIL";
			// Either the user - supplied routine psolve(used to to solve the preconditioned linear system) 
			// encountered an unrecoverable error, or the linear solver routine(lsolve) encountered an error condition.
			else if (flag_solver == KIN_LSOLVE_FAIL)		message += "KIN_LSOLVE_FAIL";
			// The system function failed in an unrecoverable manner.
			else if (flag_solver == KIN_SYSFUNC_FAIL)		message += "KIN_SYSFUNC_FAIL";
			// The system function failed recoverably at the first call.
			else if (flag_solver == KIN_FIRST_SYSFUNC_ERR)		message += "KIN_FIRST_SYSFUNC_ERR";
			// The system function had repeated recoverable errors.No recovery is possible.
			else if (flag_solver == KIN_REPTD_SYSFUNC_ERR)		message += "KIN_REPTD_SYSFUNC_ERR";

			std::cout << message << std::endl;

			// Recover original starting point
			object->CorrectedUnknownsVector(NV_DATA_S(yInitial));

			//getchar();
		}

		// Free memory
		N_VDestroy_Serial(y);
		N_VDestroy_Serial(yInitial);
		N_VDestroy_Serial(scale);
		KINFree(&kmem);
		free(data);

		return flag_solver;
	}

	static void KinSolStatistics(void *kmem)
	{
		int flag;

		// Main solver statistics
		{
			long int nni, nfe;

			flag = KINGetNumNonlinSolvIters(kmem, &nni);
			check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
			flag = KINGetNumFuncEvals(kmem, &nfe);
			check_flag(&flag, "KINGetNumFuncEvals", 1);

			std::cout << "Main solver statistics..       " << std::endl;
			std::cout << "* # linear solver iterations = " << nni << std::endl;
			std::cout << "* # function evaluations     = " << nfe << std::endl;
			std::cout << std::endl;
		}

		// Linesearch statistics
		{
			long int nbcfails, nbacktr;

			flag = KINGetNumBetaCondFails(kmem, &nbcfails);
			check_flag(&flag, "KINGetNumBetacondFails", 1);
			flag = KINGetNumBacktrackOps(kmem, &nbacktr);
			check_flag(&flag, "KINGetNumBacktrackOps", 1);

			std::cout << "Line search statistics..       " << std::endl;
			std::cout << "* # beta failures            = " << nbcfails << std::endl;
			std::cout << "* # back track operations    = " << nbacktr << std::endl;
			std::cout << std::endl;
		}

		// Band solver statistics
		{
			long int nje, nfeD;
			flag = KINDlsGetNumJacEvals(kmem, &nje);
			check_flag(&flag, "KINDlsGetNumJacEvals", 1);
			flag = KINDlsGetNumFuncEvals(kmem, &nfeD);
			check_flag(&flag, "KINDlsGetNumFuncEvals", 1);

			std::cout << "Band linear solver statistics.. " << std::endl;
			std::cout << "* # function evaluations      = " << nfeD << std::endl;
			std::cout << "* # jacobians                 = " << nje << std::endl;
			std::cout << std::endl;
		}

		// Workspace sizes
		{
			long int lenrw, leniw, lenrwB, leniwB;

			// Main solver workspace size
			flag = KINGetWorkSpace(kmem, &lenrw, &leniw);
			check_flag(&flag, "KINGetWorkSpace", 1);

			// Band linear solver workspace size
			flag = KINDlsGetWorkSpace(kmem, &lenrwB, &leniwB);
			check_flag(&flag, "KINDlsGetWorkSpace", 1);

			std::cout << "Workspaces.. " << std::endl;
			std::cout << "* main     = " << lenrw << " = " << leniw << std::endl;
			std::cout << "* band     = " << lenrwB << " = " << leniwB << std::endl;
			std::cout << std::endl;
		}
	}
}

#endif // OpenSMOKE_Solve_Band_KinSol_H


