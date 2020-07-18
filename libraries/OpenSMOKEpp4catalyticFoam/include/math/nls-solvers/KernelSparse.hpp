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

namespace NlsSMOKE
{
	template <typename NLSSystemObject>
	KernelSparse<NLSSystemObject>::KernelSparse()
	{
		solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
		preconditioner_droptol_ = 1e-6;
		preconditioner_fillfactor_ = 10;
		non_zero_elements_ = 0;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::ResetKernel()
	{
		numberOfSystemCallsForJacobian_ = 0;
		numberOfJacobianFullAssembling_ = 0;
		numberOfJacobianQuasiAssembling_ = 0;
		numberOfJacobianFactorizations_ = 0;
		numberOfLinearSystemSolutions_ = 0;

		cpuTimeJacobianFullAssembling_ = 0.;
		cpuTimeJacobianQuasiAssembling_ = 0.;
		cpuTimeJacobianFactorization_ = 0.;
		cpuTimeLinearSystemSolution_ = 0.;

		cpuTimeSingleJacobianFullAssembling_ = 0.;
		cpuTimeSingleJacobianQuasiAssembling_ = 0.;
		cpuTimeSingleJacobianFactorization_ = 0.;
		cpuTimeSingleLinearSystemSolution_ = 0.;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the NlsSystem object
		this->MemoryAllocation();
		
		// Sparsity pattern
		pattern_.resize(this->ne_, this->ne_);
		number_variables_per_equations_.resize(this->ne_);

		// Internal variables
		J_.resize(this->ne_, this->ne_);
		aux_.resize(this->ne_);
		x_plus_.resize(this->ne_);
		f_plus_.resize(this->ne_);

		// Set zero
		number_variables_per_equations_.setZero();
		J_.setZero();
		aux_.setZero();
		x_plus_.setZero();
		f_plus_.setZero();
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::MemoryAllocationSparsityPattern()
	{
		// Set sparsity pattern Jacobian
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(i_.size());
			for (unsigned int k = 0; k < i_.size(); k++)
			{
				tripletList.push_back(T(i_[k], j_[k], 1.));
			}

			J_.setFromTriplets(tripletList.begin(), tripletList.end());
			J_.makeCompressed();
		}

		// Analyze sparsity 
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			sparse_LU_.analyzePattern(J_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
			{
				sparse_bicgstab_diagonal_.analyzePattern(J_);
			}
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_bicgstab_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_bicgstab_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);

				sparse_bicgstab_ilut_.analyzePattern(J_);
			}
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
			{
				sparse_gmres_diagonal_.analyzePattern(J_);
			}
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_gmres_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_gmres_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);

				sparse_gmres_ilut_.analyzePattern(J_);
			}
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
			{
				sparse_dgmres_diagonal_.analyzePattern(J_);
			}
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_dgmres_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_dgmres_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);
				sparse_dgmres_ilut_.analyzePattern(J_);
			}
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			sparse_pardiso_.analyzePattern(J_);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.analyzePattern(J_);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.analyzePattern(J_);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			OpenSMOKE::FatalErrorMessage("The LIS sparse solver for non linear systems is not yet available. Please choose a different solver.");
		}
		#endif

		RenumberingJacobianVector();
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::RenumberingJacobianVector()
	{
		// Count nonzero elemets
		non_zero_elements_ = 0;
		for (int k = 0; k < J_.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			non_zero_elements_++;

		std::cout	<< "* Total number of nonzero elements: " << non_zero_elements_
					<< " (" << double(non_zero_elements_) / double(this->ne_*this->ne_)*100. << "%)" << std::endl;

		// Resize
		Jaux_.resize(non_zero_elements_);
		indices_.resize(non_zero_elements_);
		indices_broyden_.resize(non_zero_elements_);
		Jaux_.setZero();
		indices_.setZero();
		indices_broyden_.setZero();

		// Creates provisional indices
		Eigen::VectorXi prov_rows_;
		Eigen::VectorXi prov_cols_;
		prov_rows_.resize(non_zero_elements_);
		prov_cols_.resize(non_zero_elements_);
		prov_rows_.setZero();
		prov_cols_.setZero();

		// Build list of elements
		{
			unsigned int count = 0;
			for (int group = 1; group <= pattern_.number_groups(); group++)
			{
				for (int i = 1; i <= pattern_.number_variables_in_group(group); i++)
				{
					const int var = pattern_.variable_in_group(group, i) - 1;

					for (int k = 1; k <= pattern_.number_equations_per_variable(var + 1); k++)
					{
						const int eq = pattern_.dependencies(var + 1, k) - 1;

						prov_rows_(count) = eq;
						prov_cols_(count) = var;
						count++;
					}
				}
			}
		}

		const double tStartOrdering = OpenSMOKE::OpenSMOKEGetCpuTime();

		std::vector<unsigned int> prov_cols_ordered_(non_zero_elements_);
		for (unsigned int i = 0; i < non_zero_elements_; i++)
			prov_cols_ordered_[i] = prov_cols_(i);

		std::vector<size_t> prov_cols_ordered_indices_(non_zero_elements_);
		prov_cols_ordered_indices_ = OpenSMOKE::sort_and_track_indices_increasing(prov_cols_ordered_);


		std::vector<unsigned int> prov_rows_ordered_(non_zero_elements_);
		for (unsigned int i = 0; i < non_zero_elements_; i++)
		{
			prov_cols_ordered_[i] = prov_cols_(prov_cols_ordered_indices_[i]);
			prov_rows_ordered_[i] = prov_rows_(prov_cols_ordered_indices_[i]);
		}

		Eigen::VectorXi starting_col(this->ne_ + 1);
		starting_col.setConstant(-1);
		unsigned int flag = 0;
		for (unsigned int k = 0; k < this->ne_; k++)
		{
			for (unsigned int i = flag; i < non_zero_elements_; i++)
			{
				if (prov_cols_ordered_[i] == k)
				{
					starting_col(k) = i;
					flag = i + 1;
					break;
				}
			}
		}
		starting_col(this->ne_) = non_zero_elements_;

		const double tEndOrdering = OpenSMOKE::OpenSMOKEGetCpuTime();

		// Create indices
		std::cout << "* Renumbering Jacobian vector indices..." << std::endl;
		{
			unsigned int count = 0;
			for (int k = 0; k < J_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			{
				for (unsigned int i = starting_col(it.col()); i < starting_col(it.col() + 1); i++)
				{
					if (prov_rows_ordered_[i] == it.row())
					{
						indices_(count) = prov_cols_ordered_indices_[i];
						count++;
						break;
					}
				}
			}
		}

		// Create Broyden indices
		std::cout << "* Renumbering Jacobian vector indices (for Broyden's formula)..." << std::endl;
		{
			unsigned int count = 0;
			for (int j = 1; j <= this->ne_; j++)
			{
				for (int k = 0; k < number_variables_per_equations_(j - 1); k++)
				{
					const int i = variables_in_equations_[j](k);

					for (unsigned int ii = starting_col(i - 1); ii < non_zero_elements_; ii++)
					{
						if (prov_rows_ordered_[ii] == (j - 1))
						{
							indices_broyden_(count) = prov_cols_ordered_indices_[ii];
							count++;
							break;
						}
					}
				}
			}
		}

		const double tEndBuilding = OpenSMOKE::OpenSMOKEGetCpuTime();
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::SetSparsityPattern(const std::vector<unsigned int>& i, const std::vector<unsigned int>& j, const bool index_zero)
	{
		// Assignes the rows and the columns representing the sparsity pattern
		i_ = i;
		j_ = j;

		// In case the rows and columns are provided by using indices starting from 1
		if (index_zero == false)
		{
			for (int k = 0; k < i_.size(); k++)
				i_[k] -= 1;

			for (int k = 0; k < i_.size(); k++)
				j_[k] -= 1;
		}

		// The pattern matrix accepts only indices starting from 1
		for (int k = 0; k < i_.size(); k++)
				pattern_(i_[k] + 1, j_[k] + 1);

		// Analyzes and finds the dependencies
		pattern_.FindDependence();

		// Create additional useful variables
		{
			int ii, jj;

			pattern_.ResetScanning();
			while (pattern_.Scanning(&ii, &jj) != 0)
				number_variables_per_equations_(ii - 1)++;

			variables_in_equations_ = new Eigen::VectorXi[this->ne_ + 1];
			for (int eq = 1; eq <= this->ne_; eq++)
			{
				variables_in_equations_[eq].resize(number_variables_per_equations_(eq - 1));
				variables_in_equations_[eq].setZero();
			}

			int ll = 1;
			int kk = 0;
			pattern_.ResetScanning();
			while (pattern_.Scanning(&ii, &jj) != 0)
			{
				if (ii != ll)
				{
					ll = ii;
					kk = 0;
				}
				variables_in_equations_[ii](kk++) = jj;
			}
		}

		numberOfSystemCallsPerJacobian_ = pattern_.number_groups();

		MemoryAllocationSparsityPattern();
	}

	template <typename NLSSystemObject>
	KernelSparse<NLSSystemObject>::~KernelSparse()
	{
		delete[] variables_in_equations_;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetLinearAlgebraSolver(const OpenSMOKE::SparseSolverType linear_algebra_solver)
	{
			solverType_ = linear_algebra_solver;
			CheckLinearAlgebraSolver();
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "EigenSparseLU")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		}
		else if (linear_algebra_solver == "EigenBiCGSTAB")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB;
		}
		else if (linear_algebra_solver == "EigenGMRES")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES;
		}
		else if (linear_algebra_solver == "EigenDGMRES")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES;
		}
		else if (linear_algebra_solver == "Pardiso")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO;
		}
		else if (linear_algebra_solver == "SuperLUSerial")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL;
		}
		else if (linear_algebra_solver == "UMFPack")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK;
		}
		else if (linear_algebra_solver == "LIS")
		{
			solverType_ = OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS;
		}
		else
		{
			OpenSMOKE::ErrorMessage("KernelSparse<NLSSystemObject>", "Requested linear algebra is not supported! Available options are: EigenSparseLU || EigenBiCGSTAB || EigenGMRES || EigenDGMRES || Pardiso || SuperLUSerial || UMFPack || LIS");
		}

		CheckLinearAlgebraSolver();
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::CheckLinearAlgebraSolver()
	{
		#if OPENSMOKE_USE_MKL == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<NLSSystemObject>", "Requested Pardiso linear algebra solver is not supported because the code was compiled without MKL support!");
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<NLSSystemObject>", "Requested SuperLUSerial linear algebra solver is not supported because the code was compiled without SuperLU support!");
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<NLSSystemObject>", "Requested UMFPack linear algebra solver is not supported because the code was compiled without UMFPACK support!");
		}
		#endif
		#if OPENSMOKE_USE_LIS == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<NLSSystemObject>", "Requested LIS linear algebra solver is not supported because the code was compiled without LIS support!");
		}
		#endif
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::SetPreconditioner(const OpenSMOKE::SparsePreconditionerType preconditioner)
	{
		preconditionerType_ = preconditioner;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::SetPreconditioner(const std::string preconditioner)
	{
		if (preconditioner == "diagonal")
		{
			preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
		}
		else if (preconditioner == "ILUT")
		{
			preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_ILUT;
		}
		else
		{
			OpenSMOKE::ErrorMessage("KernelSparse<NLSSystemObject>", "Requested preconditioner is not supported! Available options are: diagonal || ILUT");
		}
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::SetPreconditionerDropTolerance(const double drop)
	{
		preconditioner_droptol_ = drop;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::SetPreconditionerFillFactor(const double fillfactor)
	{
		preconditioner_fillfactor_ = fillfactor;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		this->Jacobian(y, t, J_);

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += pattern_.number_groups();
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::NumericalJacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& f, Eigen::VectorXd& x_dimensions, const bool max_constraints, const Eigen::VectorXd& xMax)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = 1.e-8;
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);

		unsigned int count = 0;
		for (int group = 1; group <= pattern_.number_groups(); group++)
		{
			x_plus_ = x;

			if (max_constraints == false)
			{
				for (int i = 1; i <= pattern_.number_variables_in_group(group); i++)
				{
					const int var = pattern_.variable_in_group(group, i) - 1;
					const double xh = std::fabs(x(var));
					const double xdh = std::fabs(x_dimensions(var));

					double hj = ETA2*std::max(xh, xdh);
					hj = std::max(hj, ZERO_DER);
					hj = std::min(hj, 0.001 + 0.001*std::fabs(xh));
					x_plus_(var) += hj;
					aux_(var) = hj;
				}
			}
			else
			{
				for (int i = 1; i <= pattern_.number_variables_in_group(group); i++)
				{
					const int var = pattern_.variable_in_group(group, i) - 1;
					const double xh = std::fabs(x(var));
					const double xdh = std::fabs(x_dimensions(var));

					double hj = ETA2*std::max(xh, xdh);
					hj = std::max(hj, ZERO_DER);
					hj = std::min(hj, 0.001 + 0.001*std::fabs(xh));

					if (x_plus_(var) + hj > xMax(var))
						hj = -hj;

					x_plus_(var) += hj;
					aux_(var) = hj;
				}
			}

			this->Equations(x_plus_, f_plus_);

			for (int i = 1; i <= pattern_.number_variables_in_group(group); i++)
			{
				const int var = pattern_.variable_in_group(group, i) - 1;
				const double hInv = 1. / aux_(var);

				for (int k = 1; k <= pattern_.number_equations_per_variable(var + 1); k++)
				{
					const int eq = pattern_.dependencies(var + 1, k) - 1;

					Jaux_(count++) = (f_plus_(eq) - f(eq))*hInv;
				}
			}
		}

		count = 0;
		for (int k = 0; k<J_.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
		{
			it.valueRef() = Jaux_(indices_(count++));
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += pattern_.number_groups();
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::QuasiNewton(const Eigen::VectorXd& dxi, const Eigen::VectorXd& dfi)
	{
		// Updating the Jacobian matrix using the Broyden's formula
		// J(n+1) = J(n) + ( f(n+1)-f(n) - J(n)*Dx) * Dx / (Dx*Dx)

		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		{
			const double eps = 1.e-10;

			unsigned count_global = 0;
			for (int j = 1; j <= this->ne_; j++)
			{
				double dxiNorm2 = 0.;
				for (int k = 0; k < number_variables_per_equations_(j - 1); k++)
				{
					const double coeff = dxi(variables_in_equations_[j](k) - 1);
					dxiNorm2 += coeff*coeff;
				}

				if (dxiNorm2 < eps)
					continue;

				unsigned count_local = count_global;
				double sum = dfi(j - 1);
				for (int k = 0; k < number_variables_per_equations_(j - 1); k++)
				{
					const int i = variables_in_equations_[j](k);
					sum -= Jaux_(indices_broyden_(count_local++))*dxi(i - 1);
				}
				sum /= dxiNorm2;

				count_local = count_global;
				for (int k = 0; k < number_variables_per_equations_(j - 1); k++)
				{
					const int i = variables_in_equations_[j](k);
					Jaux_(indices_broyden_(count_local++)) += sum*dxi(i - 1);
				}
				count_global = count_local;
			}

			count_global = 0;
			for (int k = 0; k < J_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			{
				it.valueRef() = Jaux_(indices_(count_global++));
			}
		}
		

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianQuasiAssembling_++;
		cpuTimeSingleJacobianQuasiAssembling_ = tend - tstart;
		cpuTimeJacobianQuasiAssembling_ += cpuTimeSingleJacobianQuasiAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::DoubleProduct(Eigen::VectorXd& gi, Eigen::VectorXd& aux)
	{
		gi = J_.transpose()*aux;
		aux = J_*gi;
	}
	

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::CalculatesNormOfJacobianRows(Eigen::VectorXd& row_norms)
	{
		// There are two versions, perfectly equivalent to calculate the norms
		// of the rows of the Jacobian matrix
		// TODO: We have to investigate which is the faster between the two

		// Version 1: based on the Eigen++ functions
		/*
		{
			for (int i = 0; i < this->ne_; i++)
				weights(i) = J_.row(i).norm();
		}
		*/

		// Version 2: based on the explicit iteraton over all the nonzero elements
		{
			row_norms.setZero();
			for (int k = 0; k<J_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			{
				const double val = it.valueRef();
				row_norms(it.row()) += val*val;
			}
			for (int i = 0; i < this->ne_; i++)
				row_norms(i) = std::sqrt(row_norms(i));
		}
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::Factorize()
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			sparse_LU_.compute(J_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_bicgstab_diagonal_.compute(J_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_bicgstab_ilut_.compute(J_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_gmres_diagonal_.compute(J_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_gmres_ilut_.compute(J_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_dgmres_diagonal_.compute(J_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_dgmres_ilut_.compute(J_);
		}
#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			sparse_pardiso_.compute(J_);
		}
#endif
#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.compute(J_);
		}
#endif
#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.compute(J_);
		}
#endif
#if OPENSMOKE_USE_LIS == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			OpenSMOKE::FatalErrorMessage("The LIS sparse solver for non linear systems is not yet available. Please choose a different solver.");
		}
#endif
		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		numberOfJacobianFactorizations_++;
		cpuTimeSingleJacobianFactorization_ = tend - tstart;
		cpuTimeJacobianFactorization_ += cpuTimeSingleJacobianFactorization_;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::Solve(Eigen::VectorXd& pi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			pi = sparse_LU_.solve(pi);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				pi = sparse_bicgstab_diagonal_.solve(pi);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				pi = sparse_bicgstab_ilut_.solve(pi);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				pi = sparse_gmres_diagonal_.solve(pi);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				pi = sparse_gmres_ilut_.solve(pi);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				pi = sparse_dgmres_diagonal_.solve(pi);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				pi = sparse_dgmres_ilut_.solve(pi);
		}
#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			pi = sparse_pardiso_.solve(pi);
		}
#endif
#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			pi = sparse_superlu_serial_.solve(pi);
		}
#endif
#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			pi = sparse_umfpack_.solve(pi);
		}
#endif
#if OPENSMOKE_USE_LIS == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			OpenSMOKE::FatalErrorMessage("The LIS sparse solver for non linear systems is not yet available. Please choose a different solver.");
		}
#endif

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfLinearSystemSolutions_++;
		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeLinearSystemSolution_ += cpuTimeSingleLinearSystemSolution_;
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::NlsSolverKernelSummary(std::ostream& out)
	{
		const double totalCpu = cpuTimeJacobianFullAssembling_ + cpuTimeJacobianQuasiAssembling_ + cpuTimeJacobianFactorization_ + cpuTimeLinearSystemSolution_;
		const double totalSingleCpu = cpuTimeSingleJacobianFullAssembling_ + cpuTimeSingleJacobianQuasiAssembling_ + cpuTimeSingleJacobianFactorization_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the Sparse NLS solver Kernel" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << "Number of system calls (only to assemble Jacobian):     " << numberOfSystemCallsForJacobian_ << " (" << numberOfSystemCallsPerJacobian_ << ")" << std::endl;
		out << "Number of full Jacobian constructions (from scratch):   " << numberOfJacobianFullAssembling_ << std::endl;
		out << "Number of approximated Jacobian constructions:          " << numberOfJacobianQuasiAssembling_ << std::endl;
		out << "Number of Jacobian factorizations:                      " << numberOfJacobianFactorizations_ << std::endl;
		out << "Number of linear system solutions:                      " << numberOfLinearSystemSolutions_ << std::endl;

		out << "Cumulative CPU for constructing the Jacobian (full):    " << cpuTimeJacobianFullAssembling_ << " (" << cpuTimeJacobianFullAssembling_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU for constructing the Jacobian (approx.): " << cpuTimeJacobianQuasiAssembling_ << " (" << cpuTimeJacobianQuasiAssembling_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU for factorizing the Jacobian:            " << cpuTimeJacobianFactorization_ << " (" << cpuTimeJacobianFactorization_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU for solving the linear system:           " << cpuTimeLinearSystemSolution_ << " (" << cpuTimeLinearSystemSolution_ / totalCpu *100. << "%)" << std::endl;

		out << "CPU for constructing the Jacobian (full):               " << cpuTimeSingleJacobianFullAssembling_ << " (" << cpuTimeSingleJacobianFullAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU for constructing the Jacobian (approx.):            " << cpuTimeSingleJacobianQuasiAssembling_ << " (" << cpuTimeSingleJacobianQuasiAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU for factorizing the Jacobian:                       " << cpuTimeSingleJacobianFactorization_ << " (" << cpuTimeSingleJacobianFactorization_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU for solving the linear system:                      " << cpuTimeSingleLinearSystemSolution_ << " (" << cpuTimeSingleLinearSystemSolution_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}
}
