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

#ifndef NlsKernelSparse_H
#define NlsKernelSparse_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/src/IterativeSolvers/DGMRES.h>

#include "math/OpenSMOKE_MatrixSparsityPattern.h"

#if OPENSMOKE_USE_MKL == 1
#include <Eigen/PardisoSupport>
#endif
#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
#include <Eigen/SuperLUSupport>
#endif
#if OPENSMOKE_USE_UMFPACK == 1
#include <Eigen/UmfPackSupport>
#endif


namespace NlsSMOKE
{
	//!  A class to manage sparse Jacobian matrices associated to NLS systems
	/*!
	A class to manage sparse Jacobian matrices associated to NLS systems
	*/

	template <typename NLSSystemObject>
	class KernelSparse : public NLSSystemObject
	{
	public:

		/**
		*@brief Default constructor
		*/
		KernelSparse();

		/**
		*@brief Function to set the sparsity pattern
		*@param i vector containing the indices of rows with non-zero elements
		*@param j vector containing the indices of columns with non-zero elements
		*@param index_zero true means that indices are 0-based, false means that indices are 1-based
		*/
		void SetSparsityPattern(const std::vector<unsigned int>& i, const std::vector<unsigned int>& j, const bool index_zero);

		/**
		*@brief Sets the solver for the solution of the dense linear system
		*@param linear_algebra_solver EigenSparseLU (default) | EigenBiCGSTAB | EigenGMRES | EigenDGMRES | Pardiso | SuperLUSerial | UMFPack
		*/
		void SetLinearAlgebraSolver(const std::string linear_algebra_solver);

		/**
		*@brief Sets the solver for the solution of the dense linear system
		*@param linear_algebra_solver the type of sparse solver to use
		*/
		void SetLinearAlgebraSolver(const OpenSMOKE::SparseSolverType linear_algebra_solver);

		/**
		*@brief Function to set the preconditioner for the iterative solvers
		*@param preconditioner diagonal, ILUT
		*/
		void SetPreconditioner(const std::string preconditioner);

		/**
		*@brief Function to set the preconditioner for the iterative solvers
		*@param preconditioner diagonal, ILUT
		*/
		void SetPreconditioner(const OpenSMOKE::SparsePreconditionerType preconditioner);

		/**
		*@brief Function to set the preconditioner drop tolerance (ILUT)
		*@param drop the drop tolerance (default 1e-6)
		*/
		void SetPreconditionerDropTolerance(const double drop);

		/**
		*@brief Function to set the preconditioner fill factor (ILUT)
		*@param fillfactor the fille factor (default 10)
		*/
		void SetPreconditionerFillFactor(const double fillfactor);


		/**
		*@brief Returns the number of calls to the system of equations for assembling the Jacobian matrix
		Since the Jacobian matrix is full, the number of calls is equal to the number of Jacobian calls times the number
		of equations
		*/
		unsigned int numberOfFunctionCallsForJacobian() const { return numberOfSystemCallsForJacobian_; }

		/**
		*@brief Returns the number of function calls every time the Jacobian matrix is evaluated
		*/
		unsigned int numberFunctionsCallsPerJacobian() const { return numberOfSystemCallsPerJacobian_; }

		/**
		*@brief Returns the number of full constructions (from scratch) of Jacobian matrix
		*/
		unsigned int nnumberOfJacobianFullAssembling() const { return numberOfJacobianFullAssembling_; }

		/**
		*@brief Returns the number approximated constructions (Broyden's formula) of Jacobian matrix
		*/
		unsigned int numberOfJacobianQuasiAssembling() const { return numberOfJacobianQuasiAssembling_; }

		/**
		*@brief Returns the number of Jacobian factorizations
		*/
		unsigned int numberOfJacobianFactorizations() const { return numberOfJacobianFactorizations_; }

		/**
		*@brief Returns the number of linear system solutions
		*/
		unsigned int numberOfLinearSystemSolutions() const { return numberOfLinearSystemSolutions_; }


		/**
		*@brief Returns the cumulative CPU time (in s) for full construction (from scratch) the Jacobian matrix
		*/
		double cpuTimeJacobianFullAssembling() const { return cpuTimeJacobianFullAssembling_; }

		/**
		*@brief Returns the cumulative CPU time (in s) for approximate construction (Broyden's formula) the Jacobian matrix
		*/
		double cpuTimeJacobianQuasiAssembling() const { return cpuTimeJacobianQuasiAssembling_; }

		/**
		*@brief Returns the cumulative CPU time (in s) for factorizing the Jacobian matrix
		*/
		double cpuTimeJacobianFactorization() const { return cpuTimeJacobianFactorization_; }

		/**
		*@brief Returns the cumulative CPU time (in s) for solving the linear system
		*/
		double cpuTimeLinearSystemSolution() const { return cpuTimeLinearSystemSolution_; }


		/**
		*@brief Returns the CPU time for full construction (from scratch) the Jacobian matrix
		*/
		double cpuTimeSingleJacobianFullAssembling() const { return cpuTimeSingleJacobianFullAssembling_; }

		/**
		*@brief Returns the CPU time for approximate construction (Broyden's formula) the Jacobian matrix
		*/
		double cpuTimeSingleJacobianQuasiAssembling() const { return cpuTimeSingleJacobianQuasiAssembling_; }

		/**
		*@brief Returns the CPU time for factorizing a single Jacobian matrix
		*/
		double cpuTimeSingleJacobianFactorization() const { return cpuTimeSingleJacobianFactorization_; }

		/**
		*@brief Returns the CPU time for solving a single linear system
		*/
		double cpuTimeSingleLinearSystemSolution() const { return cpuTimeSingleLinearSystemSolution_; }


		/**
		*@brief Summary
		*/
		void NlsSolverKernelSummary(std::ostream& out);


	protected:

		/**
		*@brief Prepares the object (default option, memory allocation, etc.)
		*/
		void MemoryAllocationKernel();

		/**
		*@brief Allocates the memory for the sparsity pattern
		*/
		void MemoryAllocationSparsityPattern();

		/**
		*@brief Renumbering the Jacobian vector for performing basic operations
		*/
		void RenumberingJacobianVector();

		/**
		*@brief Prepares the object (default option, memory allocation, etc.)
		*/
		void ResetKernel();

		/**
		*@brief Checks if the selected linear solver is really available
		*/
		void CheckLinearAlgebraSolver();

		/**
		*@brief Calculates the Jacobian numerically, using the usual differentiation approach
		*@param x the current vector of dependent variables
		*@param f the current vector of rhs
		*@param x_dimensions TODO
		*@param max_constraints true if constraints on the maximum values of unknowns are imposed
		*@param yMax maximum values (if any) imposed on the unknowns
		*/
		void NumericalJacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& f, Eigen::VectorXd& x_dimensions, const bool max_constraints, const Eigen::VectorXd& xMax);

		/**
		*@brief Calculates the Jacobian analytically
		*@param y the current vector of dependent variables
		*@param t the current value of independent variable
		*/
		void UserDefinedJacobian(const Eigen::VectorXd& x, const double t);

		/**
		*@brief Factorizes the Jacobian matrix
		*/
		void Factorize();

		/**
		*@brief Solves the linear system using the factorization of the Jacobian matrix (see the Factorization() function)
		*@param pi as input variable, is the rhs of the linear system; as output variable is the solution of the linear system
		*/
		void Solve(Eigen::VectorXd& pi);

		/**
		*@brief Calculates an approximation of the Jacobian matrix using the Broyden's formula (Quasi-Newton methods)
		*@param dxi TODO
		*@param dfi TODO
		*/
		void QuasiNewton(const Eigen::VectorXd& dxi, const Eigen::VectorXd& dfi);

		/**
		*@brief Calculates the norm of each row of the Jacobian matrix
		*@param row_norms norms of rows of the Jacobian matrix
		*/
		void CalculatesNormOfJacobianRows(Eigen::VectorXd& row_norms);

		/**
		*@brief Calculates a double product to be used in the context of the gradient method
		*@param gi TODO
		*@param aux TODO
		*/
		void DoubleProduct(Eigen::VectorXd& gi, Eigen::VectorXd& aux);

		/**
		*@brief Default destructor
		*/
		~KernelSparse();

	private:

		std::vector<unsigned int> i_;
		std::vector<unsigned int> j_;

		unsigned int non_zero_elements_;
		Eigen::VectorXi *variables_in_equations_;
		Eigen::VectorXi number_variables_per_equations_;

		Eigen::SparseMatrix<double> J_;					//!< Jacobian matrix

		Eigen::VectorXd				Jaux_;				//!< Jacobian in vector form
		Eigen::VectorXi				indices_;			//!< indices to renumber the Jacobian vector
		Eigen::VectorXi				indices_broyden_;	//!< indices to renumber the Jacobian vector (when the Broyden's formula is adopted)

		Eigen::VectorXd				aux_;		//!< auxiliary vector (dimension equal to the number of equations)
		Eigen::VectorXd				x_plus_;	//!< auxiliary vector (dimension equal to the number of equations)
		Eigen::VectorXd				f_plus_;	//!< auxiliary vector (dimension equal to the number of equations)


		Eigen::SparseLU<Eigen::SparseMatrix<double> > sparse_LU_;													//!< LU solver

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_bicgstab_diagonal_;	//!< BiCGSTAB solver (diagonal)
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_gmres_diagonal_;		//!< GMRES solver (diagonal) 
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::DiagonalPreconditioner<double> > sparse_dgmres_diagonal_;		//!< DGMRES solver (diagonal) 

		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_bicgstab_ilut_;			//!< BiCGSTAB solver (ILUT)
		Eigen::GMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_gmres_ilut_;				//!< GMRES solver (ILUT) 
		Eigen::DGMRES<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > sparse_dgmres_ilut_;				//!< DGMRES solver (ILUT) 

		#if OPENSMOKE_USE_MKL == 1
				Eigen::PardisoLU<Eigen::SparseMatrix<double> > sparse_pardiso_;							//!< MKL PARDISO solver 
		#endif

		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
				Eigen::SuperLU<Eigen::SparseMatrix<double> > sparse_superlu_serial_;					//!< SuperLUSerial solver 
		#endif

		#if OPENSMOKE_USE_UMFPACK == 1
				Eigen::UmfPackLU<Eigen::SparseMatrix<double> > sparse_umfpack_;							//!< UMFPack solver 
		#endif

		#if OPENSMOKE_USE_LIS == 1
				LIS_INT *lis_ptr_;
				LIS_INT *lis_index_;
				LIS_SCALAR *lis_value_;
				LIS_MATRIX lis_G_;
				LIS_VECTOR lis_x_;
				LIS_VECTOR lis_b_;
				LIS_SOLVER lis_solver_;
		#endif

		OpenSMOKE::SparseSolverType			solverType_;			//!< sparse solver type (linear algebra)
		OpenSMOKE::SparsePreconditionerType	preconditionerType_;	//!< preconditioner type (iterative solvers)

		double preconditioner_droptol_;
		int preconditioner_fillfactor_;

		unsigned int numberOfSystemCallsForJacobian_;		//!< number of calls to the system of equation for assembling the Jacobian matrix
		unsigned int numberOfSystemCallsPerJacobian_;		//!< number of calls to the system of equation every time the Jacobian matrix is assembled
		unsigned int numberOfJacobianFullAssembling_;		//!< number of full constructions (from scratch) of Jacobian matrix
		unsigned int numberOfJacobianQuasiAssembling_;		//!< number of approximated constructions (Broyden's formula) of Jacobian matrix
		unsigned int numberOfJacobianFactorizations_;		//!< number of Jacobian factorizations
		unsigned int numberOfLinearSystemSolutions_;		//!< number of linear system solutions

		double cpuTimeJacobianFullAssembling_;			//!< cumulative CPU time for full construction (from scratch) the Jacobian matrix
		double cpuTimeJacobianQuasiAssembling_;			//!< cumulative CPU time for approximate construction (Broyden's formula) the Jacobian matrix
		double cpuTimeJacobianFactorization_;			//!< cumulative CPU time for factorizing the Jacobian matrix
		double cpuTimeLinearSystemSolution_;			//!< cumulative CPU time for solving the linear system

		double cpuTimeSingleJacobianFullAssembling_;	//!< CPU time for full construction (from scratch) the Jacobian matrix
		double cpuTimeSingleJacobianQuasiAssembling_;	//!< CPU time for approximate construction (Broyden's formula) the Jacobian matrix
		double cpuTimeSingleJacobianFactorization_;		//!< CPU time for factorizing a single Jacobian matrix
		double cpuTimeSingleLinearSystemSolution_;		//!< CPU time for solving a single linear system

		OpenSMOKE::OpenSMOKE_MatrixSparsityPattern pattern_;
	};
}

#include "KernelSparse.hpp"

#endif