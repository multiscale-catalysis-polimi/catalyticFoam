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

#include <typeinfo>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iomanip>

namespace OpenSMOKE
{
	template<typename T>
	OpenSMOKE_IDA_Sundials<T>::OpenSMOKE_IDA_Sundials(T* daeSystem)
	{
		firstCall_ = true;
		this->iUseLapack_ = true;
		calculateConsistentInitialConditions_ = true;

		this->SetDefaultValues();
		
		ySundials_ = NULL;
		ypSundials_ = NULL;
		y0Sundials_ = NULL;
		yp0Sundials_ = NULL;
		absToleranceSundials_ = NULL;
		algebraicSundials_ = NULL;
		ida_mem_ = NULL;

		A = NULL;
		LS = NULL;

		this->daeSystem_ = daeSystem;
	}

	template <typename T>
	void OpenSMOKE_IDA_Sundials<T>::SetDimensions(const int n)
	{
		MemoryAllocation(n);
	}

	template <typename T>
	void OpenSMOKE_IDA_Sundials<T>::SetLapackSolver(const bool flag)
	{
		this->iUseLapack_ = flag;
	}

	template <typename T>
	void OpenSMOKE_IDA_Sundials<T>::SetCalculateConsistentInitialConditions(const bool flag)
	{
		this->calculateConsistentInitialConditions_ = flag;
	}

	template<typename T>
	void OpenSMOKE_IDA_Sundials<T>::MemoryAllocation(const int n)
	{	
		this->n_		=	n;					// Number of equations

		this->y0_ 	= new double[this->n_];
		this->y_ 	= new double[this->n_];
		this->yp0_  = new double[this->n_];
		this->yp_   = new double[this->n_];

		y0Sundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)y0Sundials_, std::string("N_VNew_Serial"), 0))	exit(-1);

		yp0Sundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)yp0Sundials_, std::string("N_VNew_Serial"), 0))	exit(-1);

		ySundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)ySundials_, std::string("N_VNew_Serial"), 0))	exit(-1);

		ypSundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)ypSundials_, std::string("N_VNew_Serial"), 0))	exit(-1);

		this->relToleranceVector_ = new double[this->n_];
		this->absToleranceVector_ = new double[this->n_];

		absToleranceSundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)absToleranceSundials_, std::string("N_VNew_Serial"), 0))	exit(-1);

		this->algebraic_ = new bool[this->n_];
		algebraicSundials_ = N_VNew_Serial(this->n_);
		if (check_flag((void *)algebraicSundials_, std::string("N_VNew_Serial"), 0))	exit(-1);
	}

	template<typename T>
	void OpenSMOKE_IDA_Sundials<T>::AnalyzeUserOptions()
	{
		if (this->iSetMaximumNumberOfSteps_)
		{
			int flag = IDASetMaxNumSteps(ida_mem_, this->maximumNumberOfSteps_);
			if (check_flag(&flag, std::string("IDASetMaxNumSteps"), 1)) exit(-1);
		}
		else // Default value
		{
			int flag = IDASetMaxNumSteps(ida_mem_, 5000000);
			if (check_flag(&flag, std::string("IDASetMaxNumSteps"), 1)) exit(-1);
		}
		
		if (this->iSetMaximumErrorTestFailures_)
		{
			int flag = IDASetMaxErrTestFails(ida_mem_, this->maximumErrorTestFailures_);
			if (check_flag(&flag, std::string("IDASetMaxErrTestFails"), 1)) exit(-1);
		}
		
		if (this->iSetMaximumConvergenceFailures_)
		{
			int flag = IDASetMaxConvFails(ida_mem_, this->maximumConvergenceFailures_);
			if (check_flag(&flag, std::string("IDASetMaxConvFails"), 1)) exit(-1);
		}
		
		if (this->iSetMaximumNonLinearIterations_)
		{
			int flag = IDASetMaxNonlinIters(ida_mem_, this->maximumNonLinearIterations_);
			if (check_flag(&flag, std::string("IDASetMaxNonlinIters"), 1)) exit(-1);
		}

		if (this->iSetMaximumStep_)
		{
			int flag = IDASetMaxStep(ida_mem_, this->maximumStep_);
			if (check_flag(&flag, std::string("IDASetMaxStep"), 1)) exit(-1);
		}

		if (this->iSetMinimumStep_)
		{
			std::cout << "WARNING: IDA Solver does not accept a user-defined minimum step..." << std::endl;
		}

		if (this->iSetFirstStep_)
		{
			int flag = IDASetInitStep(ida_mem_, this->firstStep_);
			if (check_flag(&flag, std::string("IDASetInitStep"), 1)) exit(-1);
		}
	}

	template<typename T>
	void OpenSMOKE_IDA_Sundials<T>::Solve(const double xend, const bool* index)
	{
		for (int i = 0; i < this->n_; i++)
			this->algebraic_[i] = index[i];
		
		int flag;

		this->x_ = this->x0_;
		this->xend_ = xend;

		for (int i = 0; i<this->n_; i++)
			NV_Ith_S(y0Sundials_, i) = this->y0_[i];
		for (int i = 0; i<this->n_; i++)
			NV_Ith_S(yp0Sundials_, i) = this->yp0_[i];

		if (firstCall_ == true)
		{
			firstCall_ = false;

			/* Call IDACreate to create the solver memory and specify the 
			* Backward Differentiation Formula and the use of a Newton iteration */
			ida_mem_ = IDACreate();
			if (check_flag((void *)ida_mem_, std::string("IDACreate"), 0)) exit(-1);

			/* Call IDAInit to initialize the integrator memory and specify the
			* user's right hand side function in y'=f(t,y), the initial time t0, and
			* the initial dependent variable vector y0Sundials_. */
			flag = IDAInit(ida_mem_, this->daeSystem_->GetSystemFunctionsStatic, this->daeSystem_->GetWriteFunctionStatic, this->x0_, y0Sundials_, yp0Sundials_);
			if (check_flag(&flag, std::string("IDAInit"), 1)) exit(-1);

			/* Call IDASVtolerances to specify the scalar relative tolerance
			* and vector absolute tolerances */
			if (this->iAbsToleranceScalar_ == true)
			{
				flag = IDASStolerances(ida_mem_, this->relToleranceScalar_, this->absToleranceScalar_);
				if (check_flag(&flag, std::string("IDASVtolerances"), 1)) exit(-1);
			}
			else
			{
				for (int i = 0; i<this->n_; i++)
					NV_Ith_S(absToleranceSundials_, i) = this->absToleranceVector_[i];
				flag = IDASVtolerances(ida_mem_, this->relToleranceScalar_, absToleranceSundials_);
				if (check_flag(&flag, std::string("IDASVtolerances"), 1)) exit(-1);
			}
			
			if (this->iRelToleranceScalar_ == false)
			{
				FatalErrorMessage("The IDA solver cannot accept a vector of relative tolerances. The relative tolerance must be the same for all the equations.");
			}

			// Set algebraic equations
			for (int i = 0; i<this->n_; i++)
				NV_Ith_S(algebraicSundials_, i) = (this->algebraic_[i] == true) ? 0.0 : 1.0;
			flag = IDASetId(ida_mem_, algebraicSundials_);
			if (check_flag(&flag, std::string("IDASetId"), 1)) exit(-1);

			/* Call Solver */
			if (this->iUseLapack_ == false)
			{
				if (this->mUpper_ == 0 && this->mLower_ == 0)
				{
					std::cout << "IDA Solver: Dense Jacobian (without Lapack)..." << std::endl;

					/* Create dense SUNMatrix for use in linear solves */
					A = SUNDenseMatrix(this->n_, this->n_);
					if (check_flag((void *)A, std::string("SUNDenseMatrix"), 0)) exit(-1);

					/* Create dense SUNLinearSolver object */
					LS = SUNDenseLinearSolver(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNDenseLinearSolver"), 0)) exit(-1);
				}
				else
				{
					std::cout << "IDA Solver: Band Jacobian (without Lapack)..." << std::endl;

					/* Create banded SUNMatrix for use in linear solves -- since this will be factored,
					set the storage bandwidth to be the sum of upper and lower bandwidths */
					/* In Sundials version 4.1.x this is automatically done*/
					// A = SUNBandMatrix(this->n_, this->mUpper_, this->mLower_, (this->mUpper_ + this->mLower_));
					A = SUNBandMatrix(this->n_, this->mUpper_, this->mLower_);
					if (check_flag((void *)A, std::string("SUNBandMatrix"), 0)) exit(-1);

					/* Create banded SUNLinearSolver object for use by CVode */
					LS = SUNBandLinearSolver(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNBandLinearSolver"), 0)) exit(-1);
				}
			}
			else
			{
				if (this->mUpper_ == 0 && this->mLower_ == 0)
				{
					std::cout << "IDA Solver: Dense Jacobian (with Lapack)..." << std::endl;

					/* Create dense SUNMatrix for use in linear solves */
					A = SUNDenseMatrix(this->n_, this->n_);
					if (check_flag((void *)A, std::string("SUNDenseMatrix"), 0)) exit(-1);

					/* Create dense SUNLinearSolver object */
					LS = SUNLapackDense(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNLapackDense"), 0)) exit(-1);
				}
				else
				{
					std::cout << "IDA Solver: Band Jacobian (with Lapack)..." << std::endl;
					
					/* Create banded SUNMatrix for use in linear solves -- since this will be factored,
					set the storage bandwidth to be the sum of upper and lower bandwidths */
					/* In Sundials version 4.1.x this is automatically done*/
					//A = SUNBandMatrix(this->n_, this->mUpper_, this->mLower_, (this->mUpper_ + this->mLower_));
					A = SUNBandMatrix(this->n_, this->mUpper_, this->mLower_);
					if (check_flag((void *)A, std::string("SUNBandMatrix"), 0)) exit(-1);

					/* Create banded SUNLapackBand solver object for use by CVode */
					LS = SUNLapackBand(ySundials_, A);
					if (check_flag((void *)LS, std::string("SUNLapackBand"), 0)) exit(-1);
				}
			}

			/* Attach the matrix and linear solver */
			flag = IDADlsSetLinearSolver(ida_mem_, LS, A);
			if (check_flag(&flag, std::string("IDADlsSetLinearSolver"), 1)) exit(-1);
		}
		else
		{
			flag = IDAReInit(ida_mem_, this->x0_, y0Sundials_, yp0Sundials_);
			if (check_flag(&flag, std::string("IDAReInit"), 1)) exit(-1);
		}

		// Initialize (consistent algebraic equations are automatically found)
		if (calculateConsistentInitialConditions_ == true)
		{
			std::cout << "Calculating consistent initial conditions... " << std::endl;

			flag = IDACalcIC(ida_mem_, IDA_YA_YDP_INIT, (this->xend_ - this->x0_) / 1.e12);
			flag = IDAGetConsistentIC(ida_mem_, y0Sundials_, yp0Sundials_);

			for (int i = 0; i < this->n_; i++)
				this->y0_[i] = NV_Ith_S(y0Sundials_, i);
			for (int i = 0; i < this->n_; i++)
				this->yp0_[i] = NV_Ith_S(yp0Sundials_, i);
		}

		AnalyzeUserOptions();

		/* Solving */
		std::cout << "Solving... " << std::endl;
		this->tStart_ =  this->GetClockTime();
		flag = IDASolve(ida_mem_, this->xend_, &this->x_, ySundials_, ypSundials_, IDA_NORMAL);
		this->tEnd_ =  this->GetClockTime();

		this->x0_ = this->x_;
		for(int i=0;i<this->n_;i++)
			NV_Ith_S(y0Sundials_,i) = NV_Ith_S(ySundials_,i);
		for (int i = 0; i<this->n_; i++)
			NV_Ith_S(yp0Sundials_, i) = NV_Ith_S(ypSundials_, i);

		for(int i=0;i<this->n_;i++)
			this->y_[i] = NV_Ith_S(ySundials_,i);
		for (int i = 0; i<this->n_; i++)
			this->yp_[i] = NV_Ith_S(ypSundials_, i);
	}

	template<typename T>
	void OpenSMOKE_IDA_Sundials<T>::Status() const
	{
		int flag;
		long int nst, nfe, nsetups, netf, nni, ncfn, nje, nfeLS, nge;
		int qcurrent, qlast;
		double hcurrent, hlast;

		flag = IDAGetNumSteps(ida_mem_, &nst);
		check_flag(&flag, std::string("IDAGetNumSteps"), 1);
		flag = IDADlsGetNumJacEvals(ida_mem_, &nje);
		check_flag(&flag, std::string("IDADlsGetNumJacEvals"), 1);
		flag = IDAGetNumResEvals(ida_mem_, &nfe);
		check_flag(&flag, std::string("IDAGetNumResEvals"), 1);

		flag = IDAGetNumLinSolvSetups(ida_mem_, &nsetups);
		check_flag(&flag, std::string("IDAGetNumLinSolvSetups"), 1);
		flag = IDAGetNumErrTestFails(ida_mem_, &netf);
		check_flag(&flag, std::string("IDAGetNumErrTestFails"), 1);
		flag = IDAGetNumNonlinSolvIters(ida_mem_, &nni);
		check_flag(&flag, std::string("IDAGetNumNonlinSolvIters"), 1);
		flag = IDAGetNumNonlinSolvConvFails(ida_mem_, &ncfn);
		check_flag(&flag, std::string("IDAGetNumNonlinSolvConvFails"), 1);
		flag = IDAGetNumGEvals(ida_mem_, &nge);
		check_flag(&flag, std::string("IDAGetNumGEvals"), 1);

		flag = IDADlsGetNumResEvals(ida_mem_, &nfeLS);
		check_flag(&flag, std::string("IDADlsGetNumResEvals"), 1);

		flag = IDAGetLastOrder(ida_mem_, &qlast);
		check_flag(&flag, std::string("IDAGetLastOrder"), 1);
		flag = IDAGetCurrentOrder(ida_mem_, &qcurrent);
		check_flag(&flag, std::string("IDAGetCurrentOrder"), 1);
		flag = IDAGetLastStep(ida_mem_, &hlast);
		check_flag(&flag, std::string("IDAGetLastStep"), 1);
		flag = IDAGetCurrentStep(ida_mem_, &hcurrent);
		check_flag(&flag, std::string("IDAGetCurrentStep"), 1);


		std::cout << "IDA Sundials Status" << std::endl;
		std::cout << " * Absolute tolerance:              " << this->absToleranceScalar_   << std::endl;	// Absolute tolerance
		std::cout << " * Relative tolerance:              " << this->relToleranceScalar_   << std::endl;	// Relative tolerance
		std::cout << " * Number of steps:                 " << nst << std::endl;	// Number of steps taken for the problem so far 
		std::cout << " * Number of function evaluations:  " << nfe << std::endl;	// Number of f evaluations for the problem so far.
		std::cout << " * Number of Jacobians:             " << nje << std::endl;	// Number of Jacobian evaluations (and of matrix LU decompositions) for the problem so far.
		std::cout << " * Last step:                       " << hlast << std::endl;	
		std::cout << " * Next  step:                      " << hcurrent << std::endl;	
		std::cout << " * Last order:                      " << qlast << std::endl;	
		std::cout << " * Next order:                      " << qcurrent << std::endl;
	}

	template<typename T>
	std::string OpenSMOKE_IDA_Sundials<T>::Tag() const
	{
		return "IDA";
	}

	template<typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfSteps() const
	{
		long int nst;
		int flag = IDAGetNumSteps(ida_mem_, &nst);
		check_flag(&flag, std::string("IDAGetNumSteps"), 1);

		return int(nst);
	}

	template<typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfFunctionEvaluations() const
	{
		long int nfe;
		int flag = IDAGetNumResEvals(ida_mem_, &nfe);
		check_flag(&flag, std::string("IDAGetNumResEvals"), 1);

		return int(nfe);
	}

	template<typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfJacobianEvaluations() const
	{
		long int nje;
		int flag = IDADlsGetNumJacEvals(ida_mem_, &nje);
		check_flag(&flag, std::string("IDADlsGetNumJacEvals"), 1);

		return int(nje);
	}
	
	template<typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfFunctionForJacobianEvaluations() const
	{
		long int nfje;
		int flag = IDADlsGetNumResEvals(ida_mem_, &nfje);
		check_flag(&flag, std::string("IDADlsGetNumResEvals"), 1);

		return int(nfje);
	}

	template<typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfLUFactorizations() const
	{
		return int(0);
	}

	template <typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfNonLinearIterations() const
	{
		long int nli;
		int flag = IDAGetNumNonlinSolvIters(ida_mem_, &nli);
		check_flag(&flag, std::string("IDAGetNumNonlinSolvIters"), 1);

		return int(nli);
	}

	template <typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetLastOrderUsed() const
	{
		int qlast;
		int flag = IDAGetLastOrder(ida_mem_, &qlast);
		check_flag(&flag, std::string("IDAGetLastOrder"), 1);
    
		return qlast;
	}

	template <typename T>
	double OpenSMOKE_IDA_Sundials<T>::GetLastStepUsed() const
	{
		double hlast;
		int flag = IDAGetLastStep(ida_mem_, &hlast);
		check_flag(&flag, std::string("IDAGetLastStep"), 1);

		return hlast;
	}

	template <typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfConvergenceFailures() const
	{
		long int nlcf;
		int flag = IDAGetNumNonlinSolvConvFails(ida_mem_, &nlcf);
		check_flag(&flag, std::string("IDAGetNumNonlinSolvConvFails"), 1);
		
		return int(nlcf);
	}

	template <typename T>
	int OpenSMOKE_IDA_Sundials<T>::GetNumberOfErrorTestFailures() const
	{
		long int net;
		int flag = IDAGetNumErrTestFails(ida_mem_, &net);
		check_flag(&flag, std::string("IDAGetNumErrTestFails"), 1);
		
		return int(net);
	}

	template<typename T>
	OpenSMOKE_IDA_Sundials<T>::~OpenSMOKE_IDA_Sundials(void)
	{
		/* Free vectors */
		N_VDestroy_Serial(y0Sundials_);
		N_VDestroy_Serial(yp0Sundials_);
		N_VDestroy_Serial(ySundials_);
		N_VDestroy_Serial(ypSundials_);
		N_VDestroy_Serial(absToleranceSundials_);
		N_VDestroy_Serial(algebraicSundials_);

		/* Free integrator memory */
		IDAFree(&ida_mem_);
		SUNLinSolFree(LS);
		SUNMatDestroy(A);

		delete[] this->y0_;
		delete[] this->y_;
		delete[] this->yp0_;
		delete[] this->yp_;
		delete[] this->algebraic_;
		delete[] this->relToleranceVector_;	
		delete[] this->absToleranceVector_;

	}
}
