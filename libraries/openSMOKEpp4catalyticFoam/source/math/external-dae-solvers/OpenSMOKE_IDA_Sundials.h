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

#ifndef  OpenSMOKE_IDA_Sundials_H
#define  OpenSMOKE_IDA_Sundials_H

#include "math/external-dae-solvers/OpenSMOKE_DaeSystemSolver.h"

#include <ida/ida.h>
#include <sunmatrix/sunmatrix_dense.h>  /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>  /* access to dense SUNLinearSolver      */
#include <ida/ida_direct.h>             /* access to IDADls interface           */
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_dense.h>	/* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>

namespace OpenSMOKE
{
	#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
	#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */
	 
	static int check_flag(void *flagvalue, char *funcname, int opt);

	template <typename T>
	class OpenSMOKE_IDA_Sundials : public OpenSMOKE::OpenSMOKE_DaeSystemSolver<T>
	{
		public:

			OpenSMOKE_IDA_Sundials(T* daeSystem);

			void SetDimensions(const int n);
			void SetLapackSolver(const bool flag);
			void SetCalculateConsistentInitialConditions(const bool flag);
	
			void Solve(const double xf, const bool* index);	

			void Status() const;

			std::string Tag() const;
			int GetNumberOfSteps() const;
			int GetNumberOfFunctionEvaluations() const;
			int GetNumberOfJacobianEvaluations() const;
			int GetNumberOfFunctionForJacobianEvaluations() const;
			int GetNumberOfLUFactorizations() const;
			int GetNumberOfNonLinearIterations() const;
			int GetLastOrderUsed() const;
			int GetNumberOfConvergenceFailures() const;
			int GetNumberOfErrorTestFailures() const;
			double GetLastStepUsed() const;

			/**
			* Default destructor
			*/
			~OpenSMOKE_IDA_Sundials(void);

	private:

		N_Vector y0Sundials_;
		N_Vector yp0Sundials_;		
		N_Vector ySundials_;
		N_Vector ypSundials_;	
		N_Vector absToleranceSundials_;
		N_Vector algebraicSundials_;
		void *ida_mem_;

		SUNMatrix A;
		SUNLinearSolver LS;

		bool firstCall_;
		bool iUseLapack_;
		bool calculateConsistentInitialConditions_;

	private:

		void MemoryAllocation(const int n);
		void AnalyzeUserOptions();
	};
}

#include "OpenSMOKE_IDA_Sundials.hpp"

#endif	// OpenSMOKE_IDA_Sundials_H
