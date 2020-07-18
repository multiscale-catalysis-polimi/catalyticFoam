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
|   Copyright(C) 2016  Alberto Cuoci                                      |
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

#ifndef NlsNonLinearSolver_H
#define NlsNonLinearSolver_H

#include <Eigen/Dense>
#include "math/OpenSMOKEUtilities.h"
#include "SolverUtilities.h"

namespace NlsSMOKE
{
	//!  A class to solve non-linear systems with different techniques
	/*!
	The purpose of this class is the solution of a non-linear system of equations with different techniques:
	Newotn's method, Quasi-Newton's method, Gradient method, 1D searches
	*/

	template <typename Kernel>
	class NonLinearSolver : public Kernel
	{
	public:

		/**
		*@brief Default constructor
		*/
		NonLinearSolver();

		/**
		*@brief Set the first guess solution
		*@param x0 the vector containing the first guess solution
		*/
		void SetFirstGuessSolution(const Eigen::VectorXd& x0);

		/**
		*@brief Summary
		*@param out output stream
		*/
		void NlsSummary(std::ostream& out);

		/**
		*@brief Enables the possibility to print results and data on the screen
		*@param flag if true the possibility to print is turned on
		*/
		void SetPrint(const bool flag);

		/**
		*@brief Solves the system of non linear equations
		*/
		NlsStatus operator() (void);

		/**
		*@brief Solves the system of non linear equations
		*/
		NlsStatus operator() (const double t);

		/**
		*@brief Set the absolute tolerances (default 1e-12)
		*@param abs_tolerances the vector of requested absolute tolerances
		*/
		void SetAbsoluteTolerances(const Eigen::VectorXd& abs_tolerances);

		/**
		*@brief Set the absolute tolerances (default 1e-12)
		*@param abs_tolerances the requested absolute tolerances
		*/
		void SetAbsoluteTolerances(const double abs_tolerance);

		/**
		*@brief Set the relative tolerances
		*@param rel_tolerances the vector of requested relative tolerances
		*/
		void SetRelativeTolerances(const Eigen::VectorXd& rel_tolerances);

		/**
		*@brief Set the relative tolerances
		*@param rel_tolerances the requested relative tolerances
		*/
		void SetRelativeTolerances(const double rel_tolerance);

		/**
		*@brief Set the maximum allowed values for the dependent variables (default: none)
		*@param max_values the requested maximum values
		*/
		void SetMaximumValues(const Eigen::VectorXd& max_values);

		/**
		*@brief Set the maximum allowed values for the dependent variables (default: none)
		*@param max_values the requested maximum value
		*/
		void SetMaximumValues(const double max_value);

		/**
		*@brief Set the minimum allowed values for the dependent variables (default: none)
		*@param min_values the requested minimum values
		*/
		void SetMinimumValues(const Eigen::VectorXd& min_values);

		/**
		*@brief Set the minimum allowed values for the dependent variables (default: none)
		*@param min_value the requested minimum value
		*/
		void SetMinimumValues(const double min_value);

		/**
		*@brief Set the maximum number of allowed iterations (dafault: 10000)
		*@param n the number of allowed iterations
		*/
		void SetMaximumNumberOfIterations(const int n);

		/**
		*@brief Stop the integration because of poor convergence histry
		*@param flag true to stop
		*/
		void SetStopBecauseBadConvergence(const bool flag);

		/**
		*@brief Set the weights for calculating the merit function
		*@param w user-defined weights
		*/
		void SetWeights(const Eigen::VectorXd &w);

		/**
		*@brief Returns the current solution and residuals at the end of calculations
		*@param solution the solution
		*@param solution the residuals
		*/
		void Solution(Eigen::VectorXd& solution, Eigen::VectorXd& residuals) const;
		
		/**
		*@brief Reset the solver to reuse it
		*@param xx0 new first guess solution
		*/
		void Reset(const Eigen::VectorXd& xx0);


	private:

		/**
		*@brief Set default values
		*/
		void SetDefaultValues();

		/**
		*@brief reset the counters
		*/
		void Reset();

		/**
		*@brief Performs the initialization of the class (memory allocation)
		*/
		void MemoryAllocationMethod();

		/**
		*@brief Reset of counters and preparation of parameters (default values)
		*/
		void ResetMethod();

		/**
		*@brief Checks if the new calculated solution satisfies the provided constraints
				Violated constraints are corrected: x=xMin || x=xMax
				Returns the number of violated constraints
		*/
		unsigned int CheckConstraints(Eigen::VectorXd& xx);

		/**
		*@brief Checks if the new calculated solution satisfies the provided constraints
				Returns the number of violated constraints
		*@param gamma is the coefficient to apply for correcting the current solution: x = xMin + gamma*(xMin-x) || x = xMax + gamma*(xMax-x)
		*/
		unsigned int CheckConstraints(Eigen::VectorXd& xx, const double gamma);

		/**
		*@brief Calculates the Jacobian matrix
		*/
		void CalculateJacobian();

		/**
		*@brief Performs the Newton's step: factorization of the Jacobian matrix and solution of the linear system
		*/
		void NewtonPrevision();

		/**
		*@brief Performs the Quasi-Newton's step
		*/
		void QuasiNewtonPrevision();

		/**
		*@brief Checks the current prevision to understand if it can be accepted
		*/
		bool CheckPrevision(const double one_minus_gamma);

		/**
		*@brief Checks if convergence is satisfied, i.e. if the the non-linear system has been solved
		*/
		void CheckObjectiveFunctions();

		/**
		*@brief Calculates the weights to be used for constructing the merit function
		*/
		void Weights();

		/**
		*@brief Calculates the error weights to be used for checking the convergence of the non linear system
		*/
		void ErrorWeights();

		/**
		*@brief Calculates the merit function
		*/
		void CalculatesMeritFunction(const Eigen::VectorXd& ff, double& phi);

		/**
		*@brief Calculates the new correction (reusing the same Jacobian) and the corresponding merit function
		*/
		void CalculatesNewCorrectionAndMeritFunction();

		/**
		*@brief Calculates a proper norm of correction
		*/
		double CheckNewton();

		/**
		*@brief Modifies the first-guess solution in case of serious problems in convergence
		*/
		void ReInitialize();

		/**
		*@brief Analyze the status of the non-linear system solution
		*/
		void PrintErrorState();

	private:

		Eigen::VectorXd x0;						//!< first guess solution
		Eigen::VectorXd f0;						//!< residual of first guess solution

		Eigen::VectorXd x_measure_;				//!< order of magnitude of the first-guess solution elements
		Eigen::VectorXd	weights;				//!< weights to be used for constructing the merit function
		Eigen::VectorXd errorWeights;			//!< weights to be used for checking the final convergence of the non-linear system

		Eigen::VectorXd xi;						//!< current solution
		Eigen::VectorXd fi;						//!< current residual
		Eigen::VectorXd xi1;					//!< next solution (to be accepted or refused)
		Eigen::VectorXd fi1;					//!< next residual
		Eigen::VectorXd pi;						//!< current correction
		Eigen::VectorXd pi1;					//!< next correction
		Eigen::VectorXd dxi;					//!< auxiliary vector
		Eigen::VectorXd dfi;					//!< auxiliary vector
		
		Eigen::VectorXd aux;					//!< auxiliary vector
		Eigen::VectorXd gi;						//!< gradient method: direction
		Eigen::VectorXd gi_abs;					//!< gradient method: direction (absolute value)

		// Stop integration
		bool stopNlsIntegration_;					//!< stop integration from external environment

		// Tolerances
		Eigen::VectorXd abs_tolerances_;			//!< absolute tolerances (vector)
		Eigen::VectorXd rel_tolerances_;			//!< relative tolerances (vector)
		double abs_tolerance_;						//!< absolute tolerance (scalar)
		double rel_tolerance_;						//!< relative tolerance (scalar)
		bool abs_tolerances_scalar_;				//!< types of absolute tolerance (scalar or vector)
		bool rel_tolerances_scalar_;				//!< types of relative tolerance (scalar or vector)

		// Constraints on minimum and maximum values
		bool min_constraints_;						//!< constraints on minimum values are enabled
		bool max_constraints_;						//!< constraints on maximum values are enabled
		Eigen::VectorXd min_values_;				//!< allowed minimum values 
		Eigen::VectorXd max_values_;				//!< allowed maximum values 

		// Tolerance for the merit function
		double phi_new_absolute_tolerance_;			//!< absolute tolerance for the merit function
		double phi_new_relative_tolerance_;			//!< relative tolerance for the merit function
		double phiW_absolute_tolerance_;			//!< absolute tolerance for the merit function


		NlsStatus		nls_status_;						//!< status of the current non-linear system
		NlsMethod		nls_method_;						//!< current method (Newton, Quasi-Newton, 1D searches, Gradient)
		NlsWeights		nls_weights_type_;					//!< type of weight calculation


		int	max_number_iterations_;							//!< maximum number of iterations
		int max_total_functions_;							//!< maximum number of evalutions of system equations
		int max_calls_newton_;								//!< maximum number of newton's calls

		int	number_of_equation_system_calls_;				//!< total number of evalutions of system equations
		int number_of_equation_system_calls_cumulative_;	//!< total (cumulative) number of evalutions of system equations

		int number_iterations_;								//!< total number of iterations
		int number_calls_quasi_newtons_;					//!< total number of calls to the quasi newton's step
		int number_calls_newtons_;							//!< total number of calls to the newton's step
		int number_calls_gradients_;						//!< total number of calls to gradient step
		int number_calls_monodimensional_;					//!< total number of calls to the 1d search
		int number_trials_;									//!< total number of changes to the first guess solution

		bool is_jacobian_singular_;							//!< true if the Jacobian matrix is singular
		bool stop_because_bad_convergence_;					//!< stop in case of bad convergence history
		bool jacobian_must_be_calculated_;					//!< true if the Jacobian matrix must be recalculated
		bool printResults_;									//!< flag to enable the call to the PrintResults function

		double	phiW;				//!< merit function
		double	phiNew;				//!< merit function
		double	phi1W;				//!< merit function
		double	phi1New;			//!< merit function

		double	xiNorm2;			//!< norm of current solution
		double	sqrtInvSize;		//!< inverse of square root of number of unknoewns: 1/sqrt(N)
		
		double	u;					//!< variable used for 1D searches
		double	v;					//!< variable used for 1D searches
		double	z;					//!< variable used for 1D searches
		double	fz;					//!< variable used for 1D searches
		double	fu;					//!< variable used for 1D searches


	private:

		static const double	DEFAULT_ABSOLUTE_TOLERANCE;
		static const double	DEFAULT_RELATIVE_TOLERANCE;
		static const double	DEFAULT_PHI_RELATIVE_TOLERANCE;
		static const double	DEFAULT_PHI_ABSOLUTE_TOLERANCE;
		static const double	DEFAULT_MAX_W;

		static const double	ONE_MINUS_GAMMA_NEWTON;
		static const double	ONE_MINUS_GAMMA_QUASINEWTON;
		static const double	ONE_MINUS_GAMMA_1DSEARCH;

		static const int DEFAULT_MAX_NUMBER_ITERATIONS;
		static const int DEFAULT_MAX_NEWTONS_CALLS;

	};
}

#include "NonLinearSolver.hpp"

#endif // NlsNonLinearSolver_H
