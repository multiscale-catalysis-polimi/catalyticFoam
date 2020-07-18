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

#include "math/OpenSMOKEBandMatrix.h"

namespace NlsSMOKE
{
	#define BAND_COL(A,j) (((A->cols)[j])+(A->s_mu))
	#define BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
	#define BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->s_mu)])
	#define ROW(i,j,smu) (i-j+smu)

	template <typename NLSSystemObject>
	KernelBand<NLSSystemObject>::KernelBand()
	{
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::ResetKernel()
	{
		solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		pre_processing_ = false;

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
	void KernelBand<NLSSystemObject>::MemoryAllocationKernel()
	{
		// Allocate memory (if needed) for the DaeSystem object
		this->MemoryAllocation();

		// Internal variables
		hJ_.resize(this->ne_);
		x_plus_.resize(this->ne_);
		f_plus_.resize(this->ne_);

		// Set zero
		hJ_.setZero();
		x_plus_.setZero();
		f_plus_.setZero();
	}

	template <typename NLSSystemObject>
	KernelBand<NLSSystemObject>::~KernelBand()
	{
		J_->DestroyMat();
		J_factorized_->DestroyMat();
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::SetBandSizes(const unsigned int nUpper, const unsigned int nLower)
	{
		J_ = new OpenSMOKE::OpenSMOKEBandMatrixDouble(this->ne_, nUpper, nLower);
		if (J_ == NULL)
			OpenSMOKE::FatalErrorMessage("Memory allocation for band matrix");

		J_factorized_ = new OpenSMOKE::OpenSMOKEBandMatrixDouble(this->ne_, nUpper, nLower);
		if (J_factorized_ == NULL)
			OpenSMOKE::FatalErrorMessage("Memory allocation for band matrix");

		J_->SetToZero();
		J_factorized_->SetToZero();

		width_ = J_->nUpper() + J_->nLower() + 1;
		ngroups_ = std::min(width_, int(this->ne_));
		numberOfSystemCallsPerJacobian_ = ngroups_;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::PreProcessing()
	{
		pre_processing_ = true;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::SetLinearAlgebraSolver(const std::string linear_algebra_solver)
	{
		if (linear_algebra_solver == "Eigen")
		{
			solverType_ = OpenSMOKE::SOLVER_DENSE_EIGEN;
		}
		else
			OpenSMOKE::ErrorMessage("KernelBand<NLSSystemObject>", "Requested linear algebra is not supported!");
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::UserDefinedJacobian(const Eigen::VectorXd& y, const double t)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		// TODO
		// this->Jacobian(y, t, J_);
		OpenSMOKE::ErrorMessage("KernelBand<NLSSystemObject>", "User defined Jacobian is still not supported!");

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += this->pattern_.number_groups();
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::NumericalJacobian(const Eigen::VectorXd& x, const Eigen::VectorXd& f, Eigen::VectorXd& x_dimensions, const bool max_constraints, const Eigen::VectorXd& xMax)
	{
		if (pre_processing_ == false)
			PreProcessing();

		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const double ZERO_DER = 1.e-8;
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);

		// Save the original vector
		x_plus_ = x;

		// Loop
		for (int group = 1; group <= ngroups_; group++)
		{
			if (max_constraints == false)
			{
				for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
				{
					const double xh = std::fabs(x(j));
					const double xdh = std::fabs(x_dimensions(j));
					hJ_(j) = ETA2*std::max(xh, xdh);
					hJ_(j) = std::max(hJ_(j), ZERO_DER);
					hJ_(j) = std::min(hJ_(j), 0.001 + 0.001*std::fabs(xh));

					x_plus_(j) += hJ_(j);
				}
			}
			else
			{
				for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
				{
					const double xh = std::fabs(x(j));
					const double xdh = std::fabs(x_dimensions(j));
					hJ_(j) = ETA2*std::max(xh, xdh);
					hJ_(j) = std::max(hJ_(j), ZERO_DER);
					hJ_(j) = std::min(hJ_(j), 0.001 + 0.001*std::fabs(xh));

					if (xh + hJ_(j) > xMax(j))
						hJ_(j) = -hJ_(j);

					x_plus_(j) += hJ_(j);
				}
			}

			this->Equations(x_plus_, f_plus_);

			for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
			{
				x_plus_(j) = x(j);

				double* col_j = BAND_COL(J_, j);
				int i1 = std::max(0, j - J_->nUpper());
				int i2 = std::min(j + J_->nLower(), static_cast<int>(this->ne_ - 1));
				for (int i = i1; i <= i2; i++)
					BAND_COL_ELEM(col_j, i, j) = (f_plus_(i) - f(i)) / hJ_(j);
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfSystemCallsForJacobian_ += this->ngroups_;
		numberOfJacobianFullAssembling_++;
		cpuTimeSingleJacobianFullAssembling_ = tend - tstart;
		cpuTimeJacobianFullAssembling_ += cpuTimeSingleJacobianFullAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::DoubleProduct(Eigen::VectorXd& gi, Eigen::VectorXd& aux)
	{
		J_->TProduct(aux.data(), gi.data());
		J_->Product(gi.data(), aux.data());
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::QuasiNewton(const Eigen::VectorXd& dxi, const Eigen::VectorXd& dfi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		{
			// The auxiliary vector named x_plus is used here
			Eigen::VectorXd* normSquared = &x_plus_;
			
			normSquared->setZero();
			for (int group = 1; group <= ngroups_; group++)
			{
				for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
				{
					int i1 = std::max(0, j - J_->nUpper());
					int i2 = std::min(j + J_->nLower(), static_cast<int>(this->ne_ - 1));

					for (int i = i1; i <= i2; i++)
						(*normSquared)(i) += dxi(j)*dxi(j);
				}
			}

			// The auxiliary vector named x_plus is used here
			Eigen::VectorXd* sum_vector = &f_plus_;

			(*sum_vector) = dfi;

			for (int group = 1; group <= ngroups_; group++)
			{
				for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
				{
					double* col_j = BAND_COL(J_, j);
					int i1 = std::max(0, j - J_->nUpper());
					int i2 = std::min(j + J_->nLower(), static_cast<int>(this->ne_ - 1));

					for (int i = i1; i <= i2; i++)
						(*sum_vector)(i) -= BAND_COL_ELEM(col_j, i, j)*dxi(j);
				}
			}

			for (int group = 1; group <= ngroups_; group++)
			{
				for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
				{
					double* col_j = BAND_COL(J_, j);
					int i1 = std::max(0, j - J_->nUpper());
					int i2 = std::min(j + J_->nLower(), static_cast<int>(this->ne_ - 1));
				}
			}

			const double eps = 1.e-10;
			for (int j = 0; j < static_cast<int>(this->ne_); j++)
				(*sum_vector)(j) /= ((*normSquared)(j) + eps);

			for (int group = 1; group <= ngroups_; group++)
			{
				for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
				{
					double* col_j = BAND_COL(J_, j);
					int i1 = std::max(0, j - J_->nUpper());
					int i2 = std::min(j + J_->nLower(), static_cast<int>(this->ne_ - 1));

					for (int i = i1; i <= i2; i++)
						BAND_COL_ELEM(col_j, i, j) += (*sum_vector)(i)*dxi(j);
				}
			}
		}

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianQuasiAssembling_++;
		cpuTimeSingleJacobianQuasiAssembling_ = tend - tstart;
		cpuTimeJacobianQuasiAssembling_ += cpuTimeSingleJacobianQuasiAssembling_;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::Factorize()
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		J_->CopyTo(J_factorized_);
		J_factorized_->Factorize();

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfJacobianFactorizations_++;
		cpuTimeSingleJacobianFactorization_ = tend - tstart;
		cpuTimeJacobianFactorization_ += cpuTimeSingleJacobianFactorization_;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::Solve(Eigen::VectorXd& pi)
	{
		const double tstart = OpenSMOKE::OpenSMOKEGetCpuTime();

		J_factorized_->Solve(pi.data());

		const double tend = OpenSMOKE::OpenSMOKEGetCpuTime();

		numberOfLinearSystemSolutions_++;
		cpuTimeSingleLinearSystemSolution_ = tend - tstart;
		cpuTimeLinearSystemSolution_ += cpuTimeSingleLinearSystemSolution_;
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::CalculatesNormOfJacobianRows(Eigen::VectorXd& row_norms)
	{
		row_norms.setZero();

		for (int group = 1; group <= ngroups_; group++)
		{
			for (int j = group - 1; j < static_cast<int>(this->ne_); j += width_)
			{
				double* col_j = BAND_COL(J_, j);
				int i1 = std::max(0, j - J_->nUpper());
				int i2 = std::min(j + J_->nLower(), static_cast<int>(this->ne_ - 1));

				for (int i = i1; i <= i2; i++)
				{
					const double J = BAND_COL_ELEM(col_j, i, j);
					row_norms(i) += J*J;
				}
			}
		}

		for (int i = 0; i < static_cast<int>(this->ne_); i++)
			row_norms(i) = std::sqrt(row_norms(i));
	}

	template <typename NLSSystemObject>
	void KernelBand<NLSSystemObject>::NlsSolverKernelSummary(std::ostream& out)
	{
		const double totalCpu = cpuTimeJacobianFullAssembling_ + cpuTimeJacobianQuasiAssembling_ + cpuTimeJacobianFactorization_ + cpuTimeLinearSystemSolution_;
		const double totalSingleCpu = cpuTimeSingleJacobianFullAssembling_ + cpuTimeSingleJacobianQuasiAssembling_ + cpuTimeSingleJacobianFactorization_ + cpuTimeSingleLinearSystemSolution_;

		out << std::endl;
		out << "Data for the Banded NLS solver Kernel" << std::endl;
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
