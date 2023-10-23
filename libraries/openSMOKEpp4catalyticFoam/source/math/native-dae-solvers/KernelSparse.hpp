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

namespace DaeSMOKE
{
	template <typename DAESystemObject>
	KernelSparse<DAESystemObject>::KernelSparse()
	{
		solverType_			= OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU;
		preconditionerType_ = OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL;
		preconditioner_droptol_ = 1e-6;
		preconditioner_fillfactor_ = 10;
		non_zero_elements_ = 0;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::ResetKernel()
	{
		numberOfFunctionCallsForJacobian_ = 0;

		cpuTimeToAssembleJacobian_ = 0.;
		cpuTimeToFactorize_ = 0.;
		cpuTimeToSolveLinearSystem_ = 0.;

		cpuTimeSingleJacobianAssembling_ = 0.;
		cpuTimeSingleFactorization_ = 0.;
		cpuTimeSingleLinearSystemSolution_ = 0.;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the OdeSystem object
		this->MemoryAllocation();

		// Internal variables
		equation_type_.resize(this->ne_);
		aux_.resize(this->ne_);
		J_.resize(this->ne_, this->ne_);
		G_.resize(this->ne_, this->ne_);
		ones_.resize(this->ne_, this->ne_);
		y_plus_.resize(this->ne_);
		f_plus_.resize(this->ne_);

		// Set zero
		J_.setZero();
		aux_.setZero();
		y_plus_.setZero();
		f_plus_.setZero();

		// Memory allocation sparsity kernel
		MemoryAllocationSparsityPattern();
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
				tripletList.push_back(T(i_[k], j_[k], 1.));

			J_.setFromTriplets(tripletList.begin(), tripletList.end());
			J_.makeCompressed();
		}

		// Set sparsity pattern G Matrix (we have to be sure that diagonal elements are present)
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(i_.size() + this->ne_);
			for (unsigned int k = 0; k < i_.size(); k++)
				tripletList.push_back(T(i_[k], j_[k], 1.));
			for (unsigned int k = 0; k < this->ne_; k++)
				tripletList.push_back(T(k, k, 1.));

			G_.setFromTriplets(tripletList.begin(), tripletList.end());
			G_.makeCompressed();
		}

		// Set identity matrix
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList;
			tripletList.reserve(this->ne_);
			for (unsigned int k = 0; k < this->ne_; k++)
					tripletList.push_back(T(k, k, 1.));

			ones_.setFromTriplets(tripletList.begin(), tripletList.end());
			ones_.makeCompressed();
		}

		// Analyze sparsity 
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			sparse_LU_.analyzePattern(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_bicgstab_diagonal_.analyzePattern(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_bicgstab_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_bicgstab_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);

				sparse_bicgstab_ilut_.analyzePattern(G_);
			}
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_gmres_diagonal_.analyzePattern(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_gmres_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_gmres_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);

				sparse_gmres_ilut_.analyzePattern(G_);
			}
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_dgmres_diagonal_.analyzePattern(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
			{
				sparse_dgmres_ilut_.preconditioner().setDroptol(preconditioner_droptol_);
				sparse_dgmres_ilut_.preconditioner().setFillfactor(preconditioner_fillfactor_);
				sparse_dgmres_ilut_.analyzePattern(G_);
			}
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			sparse_pardiso_.analyzePattern(G_);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.analyzePattern(G_);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.analyzePattern(G_);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			lis_value_ = (LIS_SCALAR *)malloc(G_.nonZeros()*sizeof(LIS_SCALAR));
			lis_ptr_ = (LIS_INT *)malloc((this->ne_ + 1)*sizeof(LIS_INT));
			lis_index_ = (LIS_INT *)malloc(G_.nonZeros()*sizeof(LIS_INT));

			for (unsigned int i = 0; i <= this->ne_; i++)
				lis_ptr_[i] = 0;

			unsigned int count = 0;
			for (int k = 0; k < G_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(G_, k); it; ++it)
			{
				lis_index_[count] = it.row();
				lis_ptr_[it.col() + 1]++;
				count++;
			}

			for (unsigned int i = 1; i <= this->ne_; i++)
				lis_ptr_[i] += lis_ptr_[i - 1];

			// Memory allocation
			LIS_INT err;
			err = lis_matrix_create(LIS_COMM_WORLD, &lis_G_);
			CHKERR(err);
			err = lis_matrix_set_size(lis_G_, 0, this->ne_);
			CHKERR(err);


			lis_vector_create(LIS_COMM_WORLD, &lis_x_);
			lis_vector_set_size(lis_x_, this->ne_, 0);

			lis_vector_create(LIS_COMM_WORLD, &lis_b_);
			lis_vector_set_size(lis_b_, this->ne_, 0);


			lis_solver_create(&lis_solver_);
			lis_solver_set_option("-i gmres -p ilut -initx_zeros 0", lis_solver_);
			lis_solver_set_option("-ilut_drop 1e-8 -ilut_rate 2", lis_solver_);
		}
		#endif

		RenumberingJacobianVector();
	}

	template <typename NLSSystemObject>
	void KernelSparse<NLSSystemObject>::RenumberingJacobianVector()
	{
		// Count nonzero elemets
		non_zero_elements_ = 0;
		for (int k = 0; k<J_.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			non_zero_elements_++;

		std::cout	<< "* Total number of nonzero elements: " << non_zero_elements_
					<< "(" << double(non_zero_elements_) / double(this->ne_*this->ne_)*100. << "%)" << std::endl;

		// Resize
		Jaux_.resize(non_zero_elements_);
		indices_.resize(non_zero_elements_);
		Jaux_.setZero();
		indices_.setZero();

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
		prov_cols_ordered_indices_ = OpenSMOKE::SortAndTrackIndicesIncreasing(prov_cols_ordered_);


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
					flag = i+1;
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
				for (int i = starting_col(it.col()); i < starting_col(it.col()+1); i++)
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
		const double tEndBuilding = OpenSMOKE::OpenSMOKEGetCpuTime();

		//std::cout << tEndOrdering - tStartOrdering << " " << tEndBuilding - tEndOrdering << std::endl;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetSparsityPattern(const unsigned int neq, const std::vector<unsigned int>& i, const std::vector<unsigned int>& j, const bool index_zero)
	{

		// Sparsity pattern
		pattern_.resize(neq, neq);

		// Set
		i_ = i;
		j_ = j;

		// In case the rows and columns are provided by using indices starting from 1
		if (index_zero == false)
		{
			for (unsigned int k = 0; k < i_.size(); k++)
				i_[k] -= 1;

			for (unsigned int k = 0; k < i_.size(); k++)
				j_[k] -= 1;
		}

		// The pattern matrix accepts only indices starting from 1
		for (unsigned int k = 0; k < i_.size(); k++)
			pattern_(i_[k] + 1, j_[k] + 1);

		// Analyzes and finds the dependencies
		pattern_.FindDependence();
		
		// Number of groups
		numberOfSystemCallsPerJacobian_ = pattern_.number_groups();
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetLinearAlgebraSolver(const OpenSMOKE::SparseSolverType linear_algebra_solver)
	{
			solverType_ = linear_algebra_solver;
			CheckLinearAlgebraSolver();
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
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
			OpenSMOKE::ErrorMessage("KernelSparse<DAESystemObject>", "Requested linear algebra is not supported! Available options are: EigenSparseLU || EigenBiCGSTAB || EigenGMRES || EigenDGMRES || Pardiso || SuperLUSerial || UMFPack || LIS");
		}

		CheckLinearAlgebraSolver();
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::CheckLinearAlgebraSolver()
	{
		#if OPENSMOKE_USE_MKL == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<DAESystemObject>", "Requested Pardiso linear algebra solver is not supported because the code was compiled without MKL support!");
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<DAESystemObject>", "Requested SuperLUSerial linear algebra solver is not supported because the code was compiled without SuperLU support!");
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<DAESystemObject>", "Requested UMFPack linear algebra solver is not supported because the code was compiled without UMFPACK support!");
		}
		#endif
		#if OPENSMOKE_USE_LIS == 0
		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			OpenSMOKE::ErrorMessage("KernelSparse<DAESystemObject>", "Requested LIS linear algebra solver is not supported because the code was compiled without LIS support!");
		}
		#endif
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetPreconditioner(const std::string preconditioner)
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
			OpenSMOKE::ErrorMessage("KernelSparse<DAESystemObject>", "Requested preconditioner is not supported! Available options are: diagonal || ILUT");
		}
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetPreconditioner(const OpenSMOKE::SparsePreconditionerType preconditioner)
	{
		preconditionerType_ = preconditioner;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetPreconditionerDropTolerance(const double drop)
	{
		preconditioner_droptol_ = drop;
	}
	
	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::SetPreconditionerFillFactor(const double fillfactor)
	{
		preconditioner_fillfactor_ = fillfactor;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::JacobianTimesVector(const Eigen::VectorXd& v_in, Eigen::VectorXd* v_out)
	{
		*v_out = J_*v_in;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		OpenSMOKE::FatalErrorMessage("User defined JAcobian for OpenSMOKE++ solver not yet available");
		//this->Jacobian(y, t, J_);

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::NumericalJacobian(Eigen::VectorXd& y, const double t, const Eigen::VectorXd& f, const double h, const Eigen::VectorXd& e,
		const bool max_constraints, const Eigen::VectorXd& yMax)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double BETA = 1.e+3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;


		double hf = BETA * std::fabs(h) * OpenSMOKE::ErrorControl(f, e) * double(this->ne_);
		if (hf < 1.e-10)
			hf = 1.;

		// Jacobian vector
		{
			unsigned int count = 0;
			for (int group = 1; group <= pattern_.number_groups(); group++)
			{
				y_plus_ = y;

				if (max_constraints == false)
				{
					for (int i = 1; i <= pattern_.number_variables_in_group(group); i++)
					{
						const int var = pattern_.variable_in_group(group, i) - 1;
						const double yh = y(var);

						double hJ = ETA2*std::fabs(std::max(yh, 1. / e(var)));
						const double hJf = hf / e(var);
						hJ = std::max(hJ, hJf);
						hJ = std::max(hJ, ZERO_DER);
						hJ = std::min(hJ, 0.001 + 0.001*std::fabs(yh));

						y_plus_(var) += hJ;
						aux_(var) = hJ;
					}
				}
				else
				{
					for (int i = 1; i <= pattern_.number_variables_in_group(group); i++)
					{
						const int var = pattern_.variable_in_group(group, i) - 1;
						const double yh = y(var);

						double hJ = ETA2*std::fabs(std::max(yh, 1. / e(var)));
						const double hJf = hf / e(var);
						hJ = std::max(hJ, hJf);
						hJ = std::max(hJ, ZERO_DER);
						hJ = std::min(hJ, 0.001 + 0.001*std::fabs(yh));

						// Check the maximum value
						if (yh + hJ > yMax(var))
							hJ = -hJ;

						y_plus_(var) += hJ;
						aux_(var) = hJ;
					}
				}

				this->Equations(y_plus_, t, f_plus_);

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
		}

		// Jacobian matrix from Jacobian vector
		{
			unsigned int count = 0;
			for (int k = 0; k < J_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(J_, k); it; ++it)
			{
				it.valueRef() = Jaux_(indices_(count++));
			}
		}
		numberOfFunctionCallsForJacobian_ += this->ne_;

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleJacobianAssembling_ = tend - tstart;
		cpuTimeToAssembleJacobian_ += cpuTimeSingleJacobianAssembling_;
	}

	template <typename DAESystemObject>
	bool KernelSparse<DAESystemObject>::BuildAndFactorizeMatrixG(const double hr0, const double rr0)
	{
		// Assembling matrix G
		G_ = J_;

		// Multiplying elements (TOIMPROVE)
		for (int k = 0; k < G_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(G_, k); it; ++it)
			{
				if (this->equation_type_[it.row()] == OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL)
					it.valueRef() *= -hr0;
				else
					it.valueRef() *= -rr0;
			}

		// Adding diagonal (TOIMPROVE)
		for (int k = 0; k < G_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(G_, k); it; ++it)
			{
				if (it.row() == it.col())
					if (this->equation_type_[it.row()] == OpenSMOKE::EQUATION_TYPE_DIFFERENTIAL)
						it.valueRef() += 1.;
			}

		// Factorize
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			//sparse_LU_.compute(G_);
			sparse_LU_.factorize(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_bicgstab_diagonal_.compute(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_bicgstab_ilut_.compute(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_gmres_diagonal_.compute(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_gmres_ilut_.compute(G_);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				sparse_dgmres_diagonal_.compute(G_);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				sparse_dgmres_ilut_.compute(G_);
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
//			sparse_pardiso_.compute(G_);
			sparse_pardiso_.factorize(G_);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			sparse_superlu_serial_.compute(G_);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			sparse_umfpack_.compute(G_);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
		}
		#endif

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();
		cpuTimeSingleFactorization_ = tend - tstart;
		cpuTimeToFactorize_ += cpuTimeSingleFactorization_;

		// TODO: check if the solution was really calculated
		return true;
	}

	template <typename DAESystemObject>
	bool KernelSparse<DAESystemObject>::SolveLinearSystem(Eigen::VectorXd& db)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		Eigen::VectorXd v = db;

		if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SPARSE_LU)
		{
			db = sparse_LU_.solve(v);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_BICGSTAB)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				db = sparse_bicgstab_diagonal_.solve(v);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				db = sparse_bicgstab_ilut_.solve(v);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_GMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				db = sparse_gmres_diagonal_.solve(v);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				db = sparse_gmres_ilut_.solve(v);
		}
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_DGMRES)
		{
			if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_DIAGONAL)
				db = sparse_dgmres_diagonal_.solve(v);
			else if (preconditionerType_ == OpenSMOKE::PRECONDITIONER_SPARSE_ILUT)
				db = sparse_dgmres_ilut_.solve(v);
		}
		#if OPENSMOKE_USE_MKL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_PARDISO)
		{
			db = sparse_pardiso_.solve(v);
		}
		#endif
		#if OPENSMOKE_USE_SUPERLU_SERIAL == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_SUPERLU_SERIAL)
		{
			db = sparse_superlu_serial_.solve(v);
		}
		#endif
		#if OPENSMOKE_USE_UMFPACK == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_UMFPACK)
		{
			db = sparse_umfpack_.solve(v);
		}
		#endif
		#if OPENSMOKE_USE_LIS == 1
		else if (solverType_ == OpenSMOKE::SOLVER_SPARSE_EIGEN_LIS)
		{
			unsigned int count = 0;
			for (int k = 0; k < G_.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(G_, k); it; ++it)
			{
				lis_value_[count] = it.value();
				count++;
			}

			lis_matrix_set_csc(G_.nonZeros(), lis_ptr_, lis_index_, lis_value_, lis_G_);
			lis_matrix_assemble(lis_G_);

			for (unsigned int i = 0; i < this->ne_; i++)
				lis_vector_set_value(LIS_INS_VALUE, i, v(i), lis_b_);
			for (unsigned int i = 0; i < this->ne_; i++)
				lis_vector_set_value(LIS_INS_VALUE, i, aux_(i), lis_x_);

			lis_solve(lis_G_, lis_b_, lis_x_, lis_solver_);

			for (unsigned int i = 0; i < this->ne_; i++)
			{	
				LIS_SCALAR dummy;
				lis_vector_get_value(lis_x_, i, &dummy);
				db(i) = dummy;
			}
		}
		#endif

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeToSolveLinearSystem_ += cpuTimeSingleLinearSystemSolution_;

		// TODO: check if the solution was really calculated
		return true;
	}

	template <typename DAESystemObject>
	void KernelSparse<DAESystemObject>::DaeSolverKernelSummary(std::ostream& out)
	{
		double totalCpu = cpuTimeToAssembleJacobian_ + cpuTimeToFactorize_ + cpuTimeToSolveLinearSystem_;
		double totalSingleCpu = cpuTimeSingleFactorization_ + cpuTimeSingleJacobianAssembling_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the Dense ODE solver Kernel" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << "Number of function calls (only to assemble Jacobian): " << numberOfFunctionCallsForJacobian_ << " (" << numberOfFunctionCallsForJacobian_ / this->ne_ << ")" << std::endl;
		out << "Cumulative CPU time for assembling Jacobian:          " << cpuTimeToAssembleJacobian_ << " (" << cpuTimeToAssembleJacobian_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU time for LU decomposition:             " << cpuTimeToFactorize_ << " (" << cpuTimeToFactorize_ / totalCpu *100. << "%)" << std::endl;
		out << "Cumulative CPU time for solving the linear system:    " << cpuTimeToSolveLinearSystem_ << " (" << cpuTimeToSolveLinearSystem_ / totalCpu *100. << "%)" << std::endl;
		out << "CPU time for assembling Jacobian:                     " << cpuTimeSingleJacobianAssembling_ << " (" << cpuTimeSingleFactorization_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for LU decomposition:                        " << cpuTimeSingleFactorization_ << " (" << cpuTimeSingleJacobianAssembling_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "CPU time for solving the linear system:               " << cpuTimeSingleLinearSystemSolution_ << " (" << cpuTimeSingleLinearSystemSolution_ / totalSingleCpu *100. << "%)" << std::endl;
		out << "---------------------------------------------------------------------------------------------------------" << std::endl;
		out << std::endl;
	}
}
