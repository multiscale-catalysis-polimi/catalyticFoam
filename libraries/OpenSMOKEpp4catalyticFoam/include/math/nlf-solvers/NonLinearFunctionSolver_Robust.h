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
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

#ifndef OpenSMOKE_NonLinearFunctionSolver_Robust_H
#define OpenSMOKE_NonLinearFunctionSolver_Robust_H

#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	//!  A class for searching the roots of a non linear function
	/*!
			A class for searching the roots of a non linear function
	*/

	class NonLinearFunctionSolver_Robust
	{
		private:
	
		static const char *const OPENSMOKE_ERROR;
		static const int MAX_ITER;
		static const double TOL_REL;
		static const double TOL_ABS;
		static const double TOL_REL_Y;
		static const double TOL_ABS_Y;
		static const double DT_RATIO;
		static const double DT_MIN;
		static const double U_F_RATIO;

		char	control,unimodal;

		double	t0,y0,
				t1,y1,
				t2,y2,
				t3,y3,
				tMin,tMax,
				tBreak,
				tNew,yNew,
				tUnfeasible,
				tSolution,ySolution,
				lambda,
				lambdaMax,
				dt21,dy21,
				dt31,dy31,
				dt32,dy32,
				tTolRel,tTolAbs,yTol,
				unfeasibleFeasibleRatio,oneMunfeasibleFeasibleRatio,
				dtTol;

		OpenSMOKE::OpenSMOKEVectorDouble t,y,tCriticInf,tCriticSup;
		OpenSMOKE::OpenSMOKEVectorDouble tAux,yAux,tAux1;
		OpenSMOKE::OpenSMOKEVectorDouble x0,pi,pi1,x,x1,x2,x3;

		double (*yMono)(double t);
	//	double (*yMulti)(OpenSMOKE::OpenSMOKEVectorDouble &x);
	//	void (*fMulti)(OpenSMOKE::OpenSMOKEVectorDouble &x,OpenSMOKE::OpenSMOKEVectorDouble &f);
	//	double (*FMulti)(OpenSMOKE::OpenSMOKEVectorDouble &x);
		double Y(double t);
		double Phiw(OpenSMOKE::OpenSMOKEVectorDouble &x);

		int	iVar,
			iter,
			iterFeasible,
			iterTotal,
			iterTotalInConstructor,
			iterFeasibleInConstructor,
			iterTotalFeasible,
			nIter,
			maxIter,
			numPoints,
			nCritic,
			iVar1,iVar2,
			simplexIterations;

		char unfeasible_;

		void InitializeMono(void);
		void InitializeScanningMono(void)
		{	
			std::cout << "TODO InitializeScanningMono in BzzFunctionRootRobust (Cuoci)" << std::endl ;
			std::cout << "Press enter to continue..." << std::endl;
			getchar(); 
			exit(-1);
		}

		void Deinitialize(void);

		void FindtMax(double tt1,double ttUnfeasible);
		void FindtMin(double ttUnfeasible,double ttn);
		void FindtCritic(double tt1,double ttUnfeasible,double ttn);
		char MonoSearch(void);
		char UnconstrainedMultiSearch(void);
		void MonoPrint(void);

		int	numVariables,
			numVertices,
			numVerticesMinusOne,
			oldCase,
			countLeftRight;

		OpenSMOKE::OpenSMOKEVectorInt sorted;
		NonLinearFunctionSolver_Robust *mono;

	public:

		// default constructor
		NonLinearFunctionSolver_Robust(void);

		// copy constructor
		NonLinearFunctionSolver_Robust(const NonLinearFunctionSolver_Robust &rval);

		//	MONO_MONO
		NonLinearFunctionSolver_Robust(const double tt0, const double yy0, double (*y)(double t), const double tmin, const double tmax);
		void operator()(const double tt0, const double yy0,double (*y)(double t), const double tmin, const double tmax);

		NonLinearFunctionSolver_Robust(const double tt0, double (*y)(double t), const double tmin, const double tmax);
		void operator()(const double tt0, double (*y)(double t), const double tmin, const double tmax);

		NonLinearFunctionSolver_Robust(double (*y)(double t), const double tmin, const double tmax);
		void operator()(double (*y)(double t), const double tmin, const double tmax);

		//	MONO_SCANNING
		NonLinearFunctionSolver_Robust(const double tt0, const double yy0, const double tmin, const double tmax);
		void operator()(const double tt0, const double yy0, const double tmin, const double tmax);

		~NonLinearFunctionSolver_Robust(void)
		{
			delete mono;
		}

		double GetTSolution() const { return tSolution; }
		double GetYSolution() const { return ySolution; }
		void   GetSolutions(OpenSMOKE::OpenSMOKEVectorDouble *tSolutions, OpenSMOKE::OpenSMOKEVectorDouble *ySolutions);

		int GetIterTotalFeasible() const { return iterTotalFeasible; }
		int GetIterTotal() const  { return iterTotal; }

		char operator() () {return (*this)(maxIter);}
		char operator()(const int ni);
		void SetTolRel(const double tr);
		void SetTolAbs(const double ta);
		void SetTolY(const double yt);
		void SetMaxIter(const int mxit);

		void SetExtraPointsNumberForMultipleRootSearch(const int numP)
		{
			control = -1;
			if(numP <= 1)
				dtTol = 1.;
			else
				dtTol = 1. / double(numP);
		}

		void SetUnfeasibleFeasibleRatio(const double rat)
		{
			control = -1;
			unfeasibleFeasibleRatio = rat;
			if(unfeasibleFeasibleRatio < .001)
				unfeasibleFeasibleRatio = .001;
			else if(unfeasibleFeasibleRatio > .5)
				unfeasibleFeasibleRatio = .5;
			oneMunfeasibleFeasibleRatio = 1. - unfeasibleFeasibleRatio;
		}

		void OnlyOneRoot() 
		{
			unimodal = 1;
		}
	};
}

#include "NonLinearFunctionSolver_Robust.hpp"

#endif // OpenSMOKE_NonLinearFunctionSolver_Robust_H
