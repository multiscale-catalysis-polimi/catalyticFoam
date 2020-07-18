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
	KernelDense<NLSSystemObject>::KernelDense()
	{
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::ResetKernel()
	{
		solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		full_pivoting_ = false;

		numberOfSystemCallsForJacobian_ = 0;
		numberOfJacobianFullAssembling_   = 0;
		numberOfJacobianQuasiAssembling_  = 0;
		numberOfJacobianFactorizations_   = 0;
		numberOfLinearSystemSolutions_    = 0;

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
	void KernelDense<NLSSystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the OdeSystem object
		this->MemoryAllocation();

		// Internal variables
		aux_.resize(this->ne_);
		x_plus_.resize(this->ne_);
		J_.resize(this->ne_, this->ne_);

		numberOfSystemCallsPerJacobian_ = this->ne_;
	}

	template <typename NLSSystemObject>
	KernelDense<NLSSystemObject>::~KernelDense()
	{
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "Eigen")
		{
			solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		}
		else
			OpenSMOKE::ErrorMessage("KernelDense<NLSSystemObject>", "Requested linear algebra is not supported!");
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::SetFullPivoting(const bool flag)
	{
		full_pivoting_ = flag;
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		// TODO
		// this->Jacobian(y, t, J_);
		OpenSMOKE::ErrorMessage("KernelDense<NLSSystemObject>", "User defined Jacobian is still not supported!");

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += this->ne_;
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::NumericalJacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& f, Eigen::VectorXd& x_dimensions, const bool max_constraints, const Eigen::VectorXd& xMax)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = 1.e-8;
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);

		x_plus_ = x;

		if (max_constraints == false)
		{
			for (unsigned int j = 0; j < this->ne_; j++)
			{
				const double xh = std::fabs(x(j));
				const double xdh = std::fabs(x_dimensions(j));

				double hJ = ETA2*std::max(xh, xdh);
				hJ = std::max(hJ, ZERO_DER);
				hJ = std::min(hJ, 0.001 + 0.001*std::fabs(xh));

				x_plus_(j) += hJ;

				this->Equations(x_plus_, aux_);

				x_plus_(j) = x(j);

				aux_ -= f;

				const double hInv = 1. / hJ;
				for (unsigned int i = 0; i < this->ne_; i++)
					J_(i, j) = hInv*aux_(i);
			}
		}
		else
		{
			for (unsigned int j = 0; j < this->ne_; j++)
			{
				const double xh = std::fabs(x(j));
				const double xdh = std::fabs(x_dimensions(j));

				double hJ = ETA2*std::max(xh, xdh);
				hJ = std::max(hJ, ZERO_DER);
				hJ = std::min(hJ, 0.001 + 0.001*std::fabs(xh));

				if (xh + hJ > xMax(j))
					hJ = -hJ;

				x_plus_(j) += hJ;

				this->Equations(x_plus_, aux_);

				x_plus_(j) = x(j);

				aux_ -= f;

				const double hInv = 1. / hJ;
				for (unsigned int i = 0; i < this->ne_; i++)
					J_(i, j) = hInv*aux_(i);
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += this->ne_;
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::DoubleProduct(Eigen::VectorXd& gi, Eigen::VectorXd& aux)
	{
		gi = J_.transpose()*aux;
		aux = J_*gi;
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::QuasiNewton(const Eigen::VectorXd& dxi, const Eigen::VectorXd& dfi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		{
			const double dxiNorm = dxi.norm();
			const double dxiNormSquared = dxiNorm*dxiNorm;

			aux_ = J_*dxi;
			aux_ = dfi - aux_;

			J_ += aux_*dxi.transpose() / dxiNormSquared;
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianQuasiAssembling_++;
		cpuTimeSingleJacobianQuasiAssembling_ = tend - tstart;
		cpuTimeJacobianQuasiAssembling_ += cpuTimeSingleJacobianQuasiAssembling_;		
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::Factorize()
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_DENSE_EIGEN)
		{
			if (full_pivoting_ == true)
			{
				full_LU_.compute(J_);
			}
			else
			{
				partial_LU_.compute(J_);
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianFactorizations_++;
		cpuTimeSingleJacobianFactorization_ = tend - tstart;
		cpuTimeJacobianFactorization_ += cpuTimeSingleJacobianFactorization_;
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::Solve(Eigen::VectorXd& pi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		if (solverType_ == OpenSMOKE::SOLVER_DENSE_EIGEN)
		{
			if (full_pivoting_ == true)
			{
				pi = full_LU_.solve(pi);
			}
			else
			{
				pi = partial_LU_.solve(pi);
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfLinearSystemSolutions_++;
		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeLinearSystemSolution_ += cpuTimeSingleLinearSystemSolution_;
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::CalculatesNormOfJacobianRows(Eigen::VectorXd& row_norms)
	{
		for (int i = 0; i < this->ne_; i++)
			row_norms(i) = J_.row(i).norm();
	}


	template <typename NLSSystemObject>
	void KernelDense<NLSSystemObject>::NlsSolverKernelSummary(std::ostream& out)
	{
		const double totalCpu = cpuTimeJacobianFullAssembling_ + cpuTimeJacobianQuasiAssembling_ + cpuTimeJacobianFactorization_ + cpuTimeLinearSystemSolution_;
		const double totalSingleCpu = cpuTimeSingleJacobianFullAssembling_ + cpuTimeSingleJacobianQuasiAssembling_ + cpuTimeSingleJacobianFactorization_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the Dense NLS solver Kernel" << std::endl;
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
