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

namespace DaeSMOKE
{
	#if defined(_WIN32) || defined(_WIN64) 
		extern "C" {	void DDASPK(void RES(double*, double*, double*, double*, double*, int*, double*, int*), int *NEQ, double *T, double *Y, double *YPRIME, double *TOUT,
			int *INFO, double *RTOL, double *ATOL, int *IDID, double *RWORK, int *LRW,
			int *IWORK, int *LIW, double *RPAR, int *IPAR,
			void JAC(double*, double*, double*, double*, double*, double*, int*),
			void PSOL(int*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*, double*, int*),
			void PRINT(double *T, double *Y)
			);
		}
	#else
		extern "C" {	void ddaspk_(void RES(double*, double*, double*, double*, double*, int*, double*, int*), int *NEQ, double *T, double *Y, double *YPRIME, double *TOUT,
			int *INFO, double *RTOL, double *ATOL, int *IDID, double *RWORK, int *LRW,
			int *IWORK, int *LIW, double *RPAR, int *IPAR,
			void JAC(double*, double*, double*, double*, double*, double*, int*),
			void PSOL(int*, double*, double*, double*, double*, double*, double*, double*, double*, int*, double*, double*, int*, double*, int*),
			void PRINT(double *T, double *Y)
			);
		}
	#endif

	template<typename Object>
	int Solve_Band_Daspk(Object* object, const DaeSMOKE::DaeSolver_Parameters& parameters, const double t0, const double tEnd)
	{
		std::cout << "Band DAE solution (DASPK)..." << std::endl;

		int neq = object->NumberOfEquations();

		// Memory allocation
		const int ml = object->LowerBand();
		const int mu = object->UpperBand();
		int lrw_ = 50 + 9 * neq + (2 * ml + mu + 1)*neq + 2 * (neq / (ml + mu + 1) + 1);
		int liw_ = 40 + neq + neq;

		double* y = new double[neq];
		double* yInitial = new double[neq];
		double* yp = new double[neq];
		double* rwork_ = new double[lrw_];
		int* iwork_ = new int[liw_];
		int* ipar_ = new int[1];
		double* rpar_ = nullptr;
		double* relTolerance_;
		double* absTolerance_;

		// Default options
		int info[30];
		for (int i = 0; i < 30; i++)
			info[i] = 0;

		// Tolerances
		if (parameters.scalar_absolute_tolerance() == true && 
			parameters.scalar_relative_tolerance() == true )
		{
			info[1] = 0;
			relTolerance_ = new double[1];
			absTolerance_ = new double[1];
			relTolerance_[0] = parameters.relative_tolerance();
			absTolerance_[0] = parameters.absolute_tolerance();

			std::cout << " * Absolute tolerance (scalar): " << std::scientific << std::setprecision(3) << absTolerance_[0] << std::endl;
			std::cout << " * Relative tolerance (scalar): " << std::scientific << std::setprecision(3) << relTolerance_[0] << std::endl;
		}
		else
		{
			info[1] = 1;
			relTolerance_ = new double[neq];
			absTolerance_ = new double[neq];

			if (parameters.scalar_absolute_tolerance() == false &&
				parameters.scalar_relative_tolerance() == false)
				{

					for (int i = 0; i < neq; i++)
					{
						relTolerance_[i] = parameters.relative_tolerances()(i);
						absTolerance_[i] = parameters.absolute_tolerances()(i);
					}
				}
			else if (parameters.scalar_absolute_tolerance() == true &&
				parameters.scalar_relative_tolerance() == false)
				{
					for (int i = 0; i < neq; i++)
					{
						relTolerance_[i] = parameters.relative_tolerances()(i);
						absTolerance_[i] = parameters.absolute_tolerance();
					}
				}
			else if (parameters.scalar_absolute_tolerance() == false &&
				parameters.scalar_relative_tolerance() == true)
				{
					for (int i = 0; i < neq; i++)
					{
						relTolerance_[i] = parameters.relative_tolerance();
						absTolerance_[i] = parameters.absolute_tolerances()(i);
					}
				}

			double meanAbsTolerance = 0.;
			double meanRelTolerance = 0.;
			for (int i = 0; i < neq; i++)
			{
				meanAbsTolerance += absTolerance_[i];
				meanRelTolerance += relTolerance_[i];
			}
			meanAbsTolerance /= double(neq);
			meanRelTolerance /= double(neq);

			std::cout << " * Absolute tolerance (mean): " << std::scientific << std::setprecision(3) << meanAbsTolerance << std::endl;
			std::cout << " * Relative tolerance (mean): " << std::scientific << std::setprecision(3) << meanRelTolerance << std::endl;
		}
		
		
		// Default values
		for (int i = 0; i < lrw_; i++)
			rwork_[i] = 0.0;
		for (int i = 0; i < liw_; i++)
			iwork_[i] = 0;

		// Set Banded Jacobian
		{
			info[5] = 1;
			iwork_[0] = ml;
			iwork_[1] = mu;
		}

		// Set nonnegative constraint
		{
			if (parameters.non_negative_unknowns() == true)
				info[9] = 2;
		}

		// Set algebraic/differential equations
		{
			std::cout << "Set differential algebraic equations..." << std::endl;
			info[10] = 1;
			double* v = new double[neq];
			DaspkAlgebraicDifferentialVector(v);
			for (int i = 0; i < neq; i++)
			{
				if (v[i] == 1.)	iwork_[40 + i] = 1;
				else            iwork_[40 + i] = -1;
			}
			delete[] v;
		}

		//Initialize
		{
			object->UnknownsVector(y);
			object->UnknownsVector(yInitial);
		}

		// Assign the system
		{
			//Initialize yp vector to 0.
			DaspkInitialDerivatives(t0, y, yp);
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

		// Loop over output times, call IDASolve, and print results.	
		int flag_solver;
		boost::timer::cpu_timer timer;
		double tstart = t0;
		double tout = t_initial;
		for (unsigned int iout = 1; iout <= n_steps; iout++)
		{
			std::cout << "Integrating to: " << tout << std::endl;

			#if defined(_WIN32) || defined(_WIN64) 
				DDASPK(DaspkEquations, &neq, &tstart, y, yp, &tout,
					info, relTolerance_, absTolerance_, &flag_solver, rwork_, &lrw_, iwork_, &liw_,
					rpar_, ipar_, DaspkAnalyticalJacobian, DaspkKrylovSolver, DaspkPrintSolution);
			#else
				ddaspk_(DaspkEquations, &neq, &tstart, y, yp, &tout,
					info, relTolerance_, absTolerance_, &flag_solver, rwork_, &lrw_, iwork_, &liw_,
					rpar_, ipar_, DaspkAnalyticalJacobian, DaspkKrylovSolver, DaspkPrintSolution);
			#endif

			if (flag_solver >= 0)
			{
				double sum_abs_yp = 0.;
				for (int i = 0; i < neq; i++)
					sum_abs_yp += std::fabs(yp[i]);

				double norm2_yp = 0.;
				for (int i = 0; i < neq; i++)
					norm2_yp += yp[i] * yp[i];

				// Get scaled norm of the system function
				std::cout << " DASPK solution: " << "||y'||sum = " << sum_abs_yp << " ||y'||2 = " << norm2_yp << std::endl;

				// Check if the solution reached the steady state conditions
				if (sum_abs_yp < parameters.minimum_yp()*neq)
					break;
			}
			else
			{
				break;
			}
			
			// New end time
			if (t0 == 0.)	tout *= increasing_factor;
			else            tout += dt0*std::pow(increasing_factor, double(iout));
		}
		boost::timer::cpu_times elapsed = timer.elapsed();
		 
		/*
		IDID = 1 --a step was successfully taken in the
		C                   intermediate - output mode.The code has not
		C                   yet reached TOUT.
		C
		C           IDID = 2 --the integration to TSTOP was successfully
		C                   completed(T = TSTOP) by stepping exactly to TSTOP.
		C
		C           IDID = 3 --the integration to TOUT was successfully
		C                   completed(T = TOUT) by stepping past TOUT.
		C                   Y(*) and YPRIME(*) are obtained by interpolation.
		C
		C           IDID = 4 --the initial condition calculation, with
		C                   INFO(11) > 0, was successful, and INFO(14) = 1.
		C                   No integration steps were taken, and the solution
		C                   is not considered to have been started.
		C
		C                    *** TASK INTERRUPTED ***
		C                Reported by negative values of IDID
		C
		C           IDID = -1 --a large amount of work has been expended
		C(about 500 steps).
		C
		C           IDID = -2 --the error tolerances are too stringent.
		C
		C           IDID = -3 --the local error test cannot be satisfied
		C                     because you specified a zero component in ATOL
		C                     and the corresponding computed solution component
		C                     is zero.Thus, a pure relative error test is
		C                     impossible for this component.
		C
		C           IDID = -5 --there were repeated failures in the evaluation
		C                     or processing of the preconditioner(in JAC).
		C
		C           IDID = -6 --DDASPK had repeated error test failures on the
		C                     last attempted step.
		C
		C           IDID = -7 --the nonlinear system solver in the time
		C                     integration could not converge.
		C
		C           IDID = -8 --the matrix of partial derivatives appears
		C                     to be singular(direct method).
		C
		C           IDID = -9 --the nonlinear system solver in the time integration
		C                     failed to achieve convergence, and there were repeated
		C                     error test failures in this step.
		C
		C           IDID = -10 --the nonlinear system solver in the time integration
		C                     failed to achieve convergence because IRES was equal
		C                     to - 1.
		C
		C           IDID = -11 --IRES = -2 was encountered and control is
		C                     being returned to the calling program.
		C
		C           IDID = -12 --DDASPK failed to compute the initial Y, YPRIME.
		C
		C           IDID = -13 --unrecoverable error encountered inside user's
		C                     PSOL routine, and control is being returned to
		C                     the calling program.
		C
		C           IDID = -14 --the Krylov linear system solver could not
		C                     achieve convergence.
		C
		C           IDID = -15, .., -32 --Not applicable for this code.
		C
		C                    *** TASK TERMINATED ***
		C                reported by the value of IDID = -33
		C
		C           IDID = -33 --the code has encountered trouble from which
		C                   it cannot recover.A message is output
		C                   explaining the trouble and control is returned
		C                   to the calling program.For example, this occurs
		C                   when invalid input is detected.
		*/
		if (flag_solver > 0)
		{
			// Give results back
			object->CorrectedUnknownsVector(y);

			// Print remaining counters and free memory.
			double cpuTime = elapsed.wall / 1.e9;
		}
		else
		{
			// Recover initial solution
			object->CorrectedUnknownsVector(yInitial);
		}

		// Free memory
		delete[] y;
		delete[] yInitial;
		delete[] yp;
		delete[] rwork_;
		delete[] iwork_;
		//	delete[] ipar_;
		//	delete rpar_;
		//	delete[] relTolerance_;
		//	delete[] absTolerance_;

		return 0;
	}
}
