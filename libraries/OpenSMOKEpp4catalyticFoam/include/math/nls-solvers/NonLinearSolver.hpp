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
|   Copyright(C) 2016, 2015, 2014, 2013, 2012  Alberto Cuoci              |
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
	template <typename Kernel>
	const double NonLinearSolver<Kernel>::DEFAULT_ABSOLUTE_TOLERANCE = 1.e-10;
	
	template <typename Kernel>
	const double NonLinearSolver<Kernel>::DEFAULT_RELATIVE_TOLERANCE = std::sqrt(OpenSMOKE::MachEpsFloat());

	template <typename Kernel>
	const double NonLinearSolver<Kernel>::DEFAULT_PHI_RELATIVE_TOLERANCE = OpenSMOKE::MachEpsFloat();

	template <typename Kernel>
	const double NonLinearSolver<Kernel>::DEFAULT_PHI_ABSOLUTE_TOLERANCE = 1.e-15;

	template <typename Kernel>
	const double NonLinearSolver<Kernel>::DEFAULT_MAX_W = 1.e+30;

	template <typename Kernel>
	const double NonLinearSolver<Kernel>::ONE_MINUS_GAMMA_QUASINEWTON = 0.999;

	template <typename Kernel>
	const double NonLinearSolver<Kernel>::ONE_MINUS_GAMMA_NEWTON = 0.9999;

	template <typename Kernel>
	const double NonLinearSolver<Kernel>::ONE_MINUS_GAMMA_1DSEARCH = 0.9999;

	template <typename Kernel>
	const int NonLinearSolver<Kernel>::DEFAULT_MAX_NUMBER_ITERATIONS = 100000;

	template <typename Kernel>
	const int NonLinearSolver<Kernel>::DEFAULT_MAX_NEWTONS_CALLS = 10000;

	template <typename Kernel>
	NonLinearSolver<Kernel>::NonLinearSolver()
	{
		this->ne_ = 0;

		SetDefaultValues();
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetDefaultValues()
	{
		stopNlsIntegration_ = false;

		number_of_equation_system_calls_ = 1;
		number_of_equation_system_calls_cumulative_ = 1;

		phi_new_relative_tolerance_ = DEFAULT_PHI_RELATIVE_TOLERANCE;

		max_number_iterations_ = DEFAULT_MAX_NUMBER_ITERATIONS;

		number_iterations_ = 0;
		number_calls_quasi_newtons_ = 0;
		number_calls_newtons_ = 0;
		number_calls_gradients_ = 0;
		number_calls_monodimensional_ = 0;
		number_trials_ = 0;
		is_jacobian_singular_ = false;

		nls_status_ = NLS_INITIALIZATION_STATE;
		nls_method_ = NEWTON;
		
		jacobian_must_be_calculated_ = true;

		nls_weights_type_ = AUTOMATIC_WEIGHTS;
		max_calls_newton_ = NonLinearSolver<Kernel>::DEFAULT_MAX_NEWTONS_CALLS;
		stop_because_bad_convergence_ = false;

		// Minimum and maximum constraints
		max_constraints_ = false;
		min_constraints_ = false;

		// Tolerances
		abs_tolerances_scalar_ = true;
		rel_tolerances_scalar_ = true;
		abs_tolerance_ = DEFAULT_ABSOLUTE_TOLERANCE;
		rel_tolerance_ = DEFAULT_RELATIVE_TOLERANCE;

		// 
		printResults_ = true;
	}
	
	template <typename Kernel>
	void NonLinearSolver<Kernel>::Reset(const Eigen::VectorXd& xx0)
	{
		number_of_equation_system_calls_ = 1;
		number_of_equation_system_calls_cumulative_ = 1;

		phi_new_relative_tolerance_ = DEFAULT_PHI_RELATIVE_TOLERANCE;
		max_number_iterations_ = DEFAULT_MAX_NUMBER_ITERATIONS;

		number_iterations_ = 0;
		number_calls_quasi_newtons_ = 0;
		number_calls_newtons_ = 0;
		number_calls_gradients_ = 0;
		number_calls_monodimensional_ = 0;
		number_trials_ = 0;
		is_jacobian_singular_ = false;

		nls_status_ = NLS_INITIALIZATION_STATE;
		nls_method_ = NEWTON;
		
		jacobian_must_be_calculated_ = true;

		stop_because_bad_convergence_ = false;
		
		// Check if memory allocation is needed (first time or the number of equations is changed)
		if (this->ne_ != xx0.size())
			OpenSMOKE::FatalErrorMessage("Restart: the new first guess solution has a different size");

		// Initial value of unknowns
		x0 = xx0;
		
		// Memory allocation
		this->MemoryAllocationMethod();

		// Reset kernel
		this->ResetKernel();
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::Reset()
	{
		// TODO
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::ResetMethod()
	{
		// Kernel default values
		this->ResetKernel();
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::MemoryAllocationMethod()
	{
		// Kernel memory allocation
		this->MemoryAllocationKernel();

		// Memory allocation
		this->ne_ = x0.size();

		x_measure_ = x0;

		xi.resize(this->ne_);  xi.setZero();
		fi.resize(this->ne_);  fi.setZero();
		xi1.resize(this->ne_); xi1.setZero();
		fi1.resize(this->ne_); fi1.setZero();

		pi.resize(this->ne_);  pi.setZero();
		pi1.resize(this->ne_); pi1.setZero();
		dxi.resize(this->ne_); dxi.setZero();
		dfi.resize(this->ne_); dfi.setZero();
		f0.resize(this->ne_);  f0.setZero();
		aux.resize(this->ne_);  aux.setZero();
		gi.resize(this->ne_);  gi.setZero();
		gi_abs.resize(this->ne_);  gi_abs.setZero();

		weights.resize(this->ne_); weights.setZero();
		errorWeights.resize(this->ne_);  errorWeights.setZero();

		// Sets the current solution equal to the first-guess solution
		xi = x0;

		// Checks if the provided first guess satisfies the provided (if any) contraints
		CheckConstraints(xi);

		// Calculates the norm of the current solution
		xiNorm2 = xi.norm();

		for (unsigned int i = 0; i < this->ne_; i++)
		{
			if (x_measure_(i) == 0.)
				x_measure_(i) = 1.;
			else
				x_measure_(i) = std::fabs(x_measure_(i)*DEFAULT_RELATIVE_TOLERANCE);
		}

		double sqrtvar = std::sqrt(double(this->ne_));
		sqrtInvSize = 1. / std::sqrt(double(this->ne_));

		phi_new_absolute_tolerance_ = DEFAULT_PHI_ABSOLUTE_TOLERANCE*sqrtvar;
		phiW_absolute_tolerance_ = DEFAULT_PHI_ABSOLUTE_TOLERANCE*sqrtvar;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetFirstGuessSolution(const Eigen::VectorXd& xx0)
	{
		// Check if memory allocation is needed (first time or the number of equations is changed)
		bool memory_allocation_is_needed = false;
		if (this->ne_ != xx0.size())
			memory_allocation_is_needed = true;

		// Initial value of unknowns
		x0 = xx0;

		// Total number of equations
		this->ne_ = boost::lexical_cast<unsigned int>(x0.size());
		if (this->ne_ <= 1)
			OpenSMOKE::FatalErrorMessage("Wrong number of equations. The number of equations must be > 1");

		// Allocating local memory
		if (memory_allocation_is_needed == true)
		{
			this->MemoryAllocationMethod();
		}

		// Reset the method class
		this->ResetMethod();

		// Reset counters
		Reset();
	}
		
	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetMaximumNumberOfIterations(const int n)
	{
		max_number_iterations_ = n;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetStopBecauseBadConvergence(const bool flag)
	{
		stop_because_bad_convergence_ = flag;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::ReInitialize(void)
	{
		number_trials_++;
		
		bool j = false;
		if(number_trials_ == 1)
		{
			// The first guess solution is completely changed:
			// What was 0. becomes 1.
			for (unsigned int i = 0; i < this->ne_; i++)
			{
				if (x0(i) == 0.)
				{
					x0(i) = 1.;
					j = true;
				}
			}
		}
	
		if(j == false)
		{
			// The first guess solution is completely changed:
			// The signes are inverted
			for (unsigned int i = 0; i < this->ne_; i++)
			{
				if(number_trials_ & 1)
					x0(i) = -x0(i);
				else if(i & 1)
					x0(i) = -x0(i);
			}
		}

		// Sets the current solution equal to the first-guess solution
		xi = x0;

		// Checks if the provided first guess satisfies the provided (if any) contraints
		CheckConstraints(xi);

		// Calculates the norm of the current solution
		xiNorm2 = xi.norm();

		// Calculates the residual in the current solution
		this->Equations(xi, f0);
		number_of_equation_system_calls_++;
		number_of_equation_system_calls_cumulative_++;	
		number_calls_gradients_ = 0;

		// Update the current residuals
		fi = f0;
		
		// Update additional variables
		nls_method_ = NEWTON;
		jacobian_must_be_calculated_ = true;
		nls_weights_type_ = AUTOMATIC_WEIGHTS;
		max_calls_newton_ = NonLinearSolver<Kernel>::DEFAULT_MAX_NEWTONS_CALLS;
		stop_because_bad_convergence_ = false;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::PrintErrorState()
	{
		switch(nls_status_)
		{
			case NLS_EXCESSIVE_WORK_STATE:
				std::cout << "Impossible to reach the solution because of the excessive number of iterations" << std::endl;
				break;	

			case NLS_STOP_FOUND:
				std::cout << "The solution of the non linear system was stopped by the user" << std::endl;
				break;	

			case NLS_NOT_INITIALIZED:
				std::cout << "The non-linear system was not initialized" << std::endl;
				break;

				break;
			case NLS_STOP_FOR_BAD_CONVERGENCE:
				std::cout << "Impossible to reach the olution because of the bad convergence history" << std::endl;
				break;
		}
	}

	template <typename Kernel>
	unsigned int NonLinearSolver<Kernel>::CheckConstraints(Eigen::VectorXd& xx, const double gamma)
	{
		unsigned int violated_constraints = 0;
	
		if (min_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			{
				if (xx(i) < min_values_(i))
				{
					violated_constraints++;
					xx(i) = min_values_(i) + gamma*(min_values_(i) - xx(i));
				}
			}
		}

		if (max_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			{
				if (xx(i) > max_values_(i))
				{
					violated_constraints++;
					xx(i) = max_values_(i) + gamma*(max_values_(i) - xx(i));
				}
			}
		}

		return violated_constraints;
	}


	template <typename Kernel>
	unsigned int NonLinearSolver<Kernel>::CheckConstraints(Eigen::VectorXd& xx)
	{
		unsigned int violated_constraints = 0;

		if (min_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			{
				if (xx(i) < min_values_(i))
				{
					violated_constraints++;
					xx(i) = min_values_(i);
				}
			}
		}

		if (max_constraints_ == true)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
			{
				if (xx(i) > max_values_(i))
				{
					violated_constraints++;
					xx(i) = max_values_(i);
				}
			}
		}

		return violated_constraints;
	}


	template <typename Kernel>
	void NonLinearSolver<Kernel>::ErrorWeights()
	{
		if (abs_tolerances_scalar_ == true && rel_tolerances_scalar_ == true)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				errorWeights(i) = abs_tolerance_ + rel_tolerance_*xi(i);
		}
		else if (abs_tolerances_scalar_ == false && rel_tolerances_scalar_ == true)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				errorWeights(i) = abs_tolerances_(i) + rel_tolerance_*xi(i);
		}
		else if (abs_tolerances_scalar_ == true && rel_tolerances_scalar_ == false)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				errorWeights(i) = abs_tolerance_ + rel_tolerances_(i)*xi(i);
		}
		else if (abs_tolerances_scalar_ == false && rel_tolerances_scalar_ == false)
		{
			for (unsigned int i = 0; i<this->ne_; i++)
				errorWeights(i) = abs_tolerances_(i) + rel_tolerances_(i)*xi(i);
		}
	
		for (unsigned int i = 0; i<this->ne_; i++)
		{
			const double den = errorWeights(i);

			if(den == 0.)
				OpenSMOKE::FatalErrorMessage("Numerical Problems.""\nTry with Improve Float Consistency");
			errorWeights(i) = 1. / den;
		}
	}


	template <typename Kernel>
	void NonLinearSolver<Kernel>::CheckObjectiveFunctions()
	{						 
		if ((phiNew < xiNorm2*phi_new_relative_tolerance_ || phiNew < phi_new_absolute_tolerance_) && is_jacobian_singular_ == false)
			nls_status_ = NLS_PHINEW_OK_STATE;

		if(phiW < phiW_absolute_tolerance_)	
		{
			nls_status_ = NLS_PHIW_OK_STATE;

			if(weights.minCoeff() == 1.e-8)
				nls_status_ = NLS_DUBIOUS_STATE;
		}
	}


	template <typename Kernel>
	double NonLinearSolver<Kernel>::CheckNewton()
	{
		for(unsigned int i = 0;i < this->ne_; i++)
			aux(i) = pi(i)*errorWeights(i);
		return sqrtInvSize*aux.norm();
	}


	template <typename Kernel>
	bool NonLinearSolver<Kernel>::CheckPrevision(const double one_minus_gamma)
	{
		// Check if the Newton's prevision can be accepted or not
		if ((phi1W <= phiW*one_minus_gamma) || phi1New <= phiNew*one_minus_gamma && is_jacobian_singular_ == false)
		{
			xi = xi1;
			fi = fi1;
			pi = pi1;
			phiW = phi1W;
			phiNew = phi1New;
			xiNorm2 = xi.norm();

			return true;
		}

		return false;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::CalculateJacobian()
	{
		this->NumericalJacobian(xi, fi, x_measure_, false, max_values_);

		number_of_equation_system_calls_ += this->numberFunctionsCallsPerJacobian();
		number_of_equation_system_calls_cumulative_ += this->numberFunctionsCallsPerJacobian();
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::NewtonPrevision()
	{
		// Newton's step: J*deltax =-f (25.31, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 550)
		// Here pi is both the RHS and the solution, i.e. the correction deltax is pi
		pi = -fi;
		this->Factorize();
		this->Solve(pi);

		// The Jacobian matrix is not singular
		is_jacobian_singular_ = false;

		// Scaling of the correction step
		// The correction step is scaled if ||pi|| > (1+||x||)
		phiNew = pi.norm();			 
		if(phiNew > 1. + xiNorm2)
		{
			pi /= (phiNew / (1. + xiNorm2));
			phiNew = pi.norm();
		}
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::QuasiNewtonPrevision()
	{
		pi = -fi;

	
		dxi = xi - dxi;
		dfi = fi - dfi;
		this->QuasiNewton(dxi, dfi);
		this->Factorize();
		this->Solve(pi);

		phiNew = pi.norm();
		if(phiNew > 1. + xiNorm2)
		{
			pi /= (phiNew / (1. + xiNorm2));
			phiNew = pi.norm();
		}
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::Weights()
	{
		// Calculation of norms of Jacobian rows
		this->CalculatesNormOfJacobianRows(weights);

		// Weights are defined as the inverse of the Jacobian rows
		// See 25.11, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 545
		{
			double wmax = weights.maxCoeff();

			if (wmax > DEFAULT_MAX_W)
			{
				for (unsigned int i = 0; i < this->ne_; i++)
					if (weights(i) > DEFAULT_MAX_W)
						weights(i) = DEFAULT_MAX_W;
				wmax = DEFAULT_MAX_W;
			}

			if (wmax == 0.)
			{
				for (unsigned int i = 0; i < this->ne_; i++)
					weights(i) = 1.;
			}
			else
			{
				double wmin = weights.minCoeff();

				if (wmin == 0.)
					wmin = 1.;

				for (unsigned int i = 0; i < this->ne_; i++)
				{
					if (weights(i) == 0.)
						weights(i) = 1.;
					else
						weights(i) = wmin / weights(i);

					if (weights(i) < 1.e-8)
						weights(i) = 1.e-8;
				}
			}
		}
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::CalculatesNewCorrectionAndMeritFunction()
	{
		pi1 = -fi1;
		this->Solve(pi1);
		phi1New = pi1.norm();
	}								
	
	template <typename Kernel>
	void NonLinearSolver<Kernel>::CalculatesMeritFunction(const Eigen::VectorXd& ff, double& phi)
	{
		if (nls_weights_type_ != UNITARY_WEIGHTS)
		{
			for (unsigned int i = 0; i < this->ne_; i++)
				aux(i) = ff(i)*weights(i);
			phi = aux.norm();
		}
		else
		{
			phi = ff.norm();
		}

		phi *= 0.5*phi;
	}


	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetMaximumValues(const Eigen::VectorXd& max_values)
	{
		if (max_values.size() != this->ne_)
			OpenSMOKE::FatalErrorMessage("The size of the maximum value vector is wrong");

		this->max_values_.resize(this->ne_);
		this->max_values_ = max_values;
		this->max_constraints_ = true;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetMaximumValues(const double max_value)
	{
		this->max_values_.resize(this->ne_);
		this->max_values_.setConstant(max_value);
		this->max_constraints_ = true;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetMinimumValues(const Eigen::VectorXd& min_values)
	{
		if (min_values.size() != this->ne_)
			OpenSMOKE::FatalErrorMessage("The size of the maximum value vector is wrong");

		this->min_values_.resize(this->ne_);
		this->min_values_ = min_values;
		this->min_constraints_ = true;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetMinimumValues(const double min_value)
	{
		this->min_values_.resize(this->ne_);
		this->min_values_.setConstant(min_value);
		this->min_constraints_ = true;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::Solution(Eigen::VectorXd& solution, Eigen::VectorXd& residuals) const
	{
		solution = xi;
		residuals = fi;
	}

	template <typename Kernel>
	NlsStatus NonLinearSolver<Kernel>::operator()(void)
	{
		return (*this)(0.);
	}

	template <typename Kernel>
	NlsStatus NonLinearSolver<Kernel>::operator()(const double t)
	{
		// Set auxialiary vectors
		Eigen::VectorXd aux1(this->ne_);

		// The functions are called only the very first time
		if(number_of_equation_system_calls_cumulative_ == 1)	
		{
			// Calculate the residual and calculate the CPU time for calling a single system of equations
			// Maybe it can be useful in the fuure
			double cpuTimeEquationSystem_ = OpenSMOKE::OpenSMOKEGetCpuTime();
			this->Equations(xi,f0);
			cpuTimeEquationSystem_ = OpenSMOKE::OpenSMOKEGetCpuTime() - cpuTimeEquationSystem_;
			
			// Update the current residuals
			fi = f0;
		}

		bool control;
	
		// Stop integration if requested by the user through the extern variable
		if (stopNlsIntegration_ == true)
			OpenSMOKE::FatalErrorMessage("The solution of the non linear system is not performed because the stopNlsIntegration variable is set equal to true");

		// Stop the calculations if previous calls failed
		if(nls_status_ < 0)
		{
			PrintErrorState();
			return nls_status_;
		}		

		// Check if previous calls (if any) of the solver were completed successfully
		if (nls_status_ == NLS_INITIALIZATION_STATE)
		{
			nls_status_ = NLS_CONTINUATION_STATE;
		}
		else if(nls_status_ != NLS_CONTINUATION_STATE)
		{
			return nls_status_;
		}
		
		// Reset the number of iteration calls
		number_of_equation_system_calls_ = 0;

		for (number_iterations_ = 1; number_iterations_ <= max_number_iterations_; number_iterations_++)
		{
			// Call the print function (defined in the NLS object)
			if (printResults_ == true)
			{
				this->Print(number_iterations_, t, phiW, xi, fi);
			}

			// Stop the integration if requested by the user
			if(stopNlsIntegration_ == true)
			{
				nls_status_ = NLS_STOP_FOUND;
				return nls_status_;
			}
		
			// Calculates the error weights to check if the convergence of the non-linear system is reached
			ErrorWeights();

			// Select the proper method
			switch(nls_method_)
			{
				case NEWTON:		

					// Returns an error message if the maximum number of Newton's calls has been reached
					if(number_calls_newtons_ >= max_calls_newton_)
					{
						nls_status_ = NLS_MAX_NEWTON_CALLS;
						return nls_status_;
					}

					// Update the counters of Newton's calls
					number_calls_newtons_++;

					// Store the current solution and the residual function
					dxi = xi;
					dfi = fi;

					// Updating the Jacobian matrix
					if (jacobian_must_be_calculated_ == true)
						CalculateJacobian();
					jacobian_must_be_calculated_ = true;

					// Factorize and solve the linear system corresponding to the Newton's step 
					NewtonPrevision();

					// Calculation of weights for calculation of merit function phi
					if(nls_weights_type_ == AUTOMATIC_WEIGHTS)
						Weights();

					// Calculation of merit function: phi = f*W*W*f' (25.6, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 544)
					CalculatesMeritFunction(fi, phiW);
	 		
					// Applies the Newton's correction
					xi1 = xi + pi;

					// Check the Newton's step
					// Estimation of a proper norm of the correction step
					double controlNewton; controlNewton = CheckNewton();

					// Check if the new prevision satisfies the provided constraints (if any)
					{
						const unsigned int violated_constraints = CheckConstraints(xi1);;
						if (violated_constraints != 0 && printResults_ == true)
							std::cout << "NonLinearSolver<Kernel>: number of constraints violated in Newton's step = " << violated_constraints << std::endl;
					}

					// Calculates the residuals corresponding to the new prevision
					this->Equations(xi1, fi1);
					number_of_equation_system_calls_++; 
					number_of_equation_system_calls_cumulative_++;

					// 1. Reuse the same Jacobian matrix (already factorized) in order to perform a new Newton's step
					// 2. Calculates the new merit function in the new prevision
					CalculatesNewCorrectionAndMeritFunction();
					CalculatesMeritFunction(fi1, phi1W);

					// If the proper norm of the correction step is sufficiently small, this means that the Newton's step is satisfactory
					if (controlNewton < 1. && is_jacobian_singular_ == false)
					{
						// If the merit function is smaller, accept the second prevision (obtained reusing the same Jacobian)
						if(phi1W < phiW)
						{
							xi = xi1;
							fi = fi1;
							phiW = phi1W;
						}

						nls_status_ = NLS_NEWTON_OK_STATE;

						return nls_status_;
					}

					// Check if the prevision can be accepted (true) or not (false)
					control = CheckPrevision(NonLinearSolver<Kernel>::ONE_MINUS_GAMMA_NEWTON);

					// The test was successfull (i.e. the prevision can be accepted)
					if(control == true)
					{
						// Checks if the residuals are sufficiently small
						// If yes, this means the algorithm reached the solution
						CheckObjectiveFunctions();

						if(nls_status_ != NLS_CONTINUATION_STATE)
							return nls_status_;

						if (is_jacobian_singular_ == false)
							nls_method_ = QUASI_NEWTON;

						break;
					}
							
					// If the Newton's step failed, a 1D searh is performed in order to have a more appropriate first guess solution
					nls_method_ = ONED_SEARCH_FIRSTCALL;

					break;

				case QUASI_NEWTON:

					number_calls_quasi_newtons_++;

					// Quasi-Newton prediction
					QuasiNewtonPrevision();

					// Applicates the calculated correction
					xi1 = xi + pi;

					// Check the Quasi-Newton's step
					// Estimation of a proper norm of the correction step
					double controlQuasiNewton; controlQuasiNewton = CheckNewton();

					// Check if the new prevision satisfies the provided constraints (if any)
					{
						unsigned int violated_constraints = CheckConstraints(xi1);;
					//	if (violated_constraints != 0)
					//		std::cout << "NonLinearSolver<Kernel>: number of contraints violated in Quasi-Newton's step = " << violated_constraints << std::endl;
					}

					// Calculates the residuals corresponding to the new prevision
					this->Equations(xi1, fi1);
					number_of_equation_system_calls_++; 
					number_of_equation_system_calls_cumulative_++;

					// 1. Reuse the same Jacobian matrix (already factorized) in order to perform a new Newton's step
					// 2. Calculates the new merit function in the new prevision
					CalculatesNewCorrectionAndMeritFunction();
					CalculatesMeritFunction(fi1, phi1W);
					
					// If the proper norm of the correction step is sufficiently small, this means that the Newton's step is satisfactory
					if (controlQuasiNewton < 0.1)
					{
						// If the Jacobian matrix is singular, it is more appropriate to perform a Newton's step
						if (is_jacobian_singular_ == true)
						{
							nls_method_ = NEWTON;
							break;
						}
						else
						{
							// If the merit function is smaller, accept the Quais-Newton's step
							if(phi1W < phiW)
							{
								xi = xi1;
								fi = fi1;
								phiW = phi1W;
							}
						}

						nls_status_ = NLS_QUASI_NEWTON_OK_STATE;

						return nls_status_;
					}		
					
					// Store the current solution and residuals
					dxi = xi;
					dfi = fi;

					// Check if the prevision can be accepted (true) or not (false)
					control = CheckPrevision(NonLinearSolver<Kernel>::ONE_MINUS_GAMMA_QUASINEWTON);

					// The test was successfull (i.e. the prevision can be accepted)
					if(control == true)
					{
						// Checks if the residuals are sufficiently small
						// If yes, this means the algorithm reached the solution
						CheckObjectiveFunctions();

						if(nls_status_ != NLS_CONTINUATION_STATE)
							return nls_status_;
						
						nls_method_ = QUASI_NEWTON;
						
						break;
					}
							
					// If the Quasi-Newton's step failed, the Newton's step is performed
					nls_method_ = NEWTON;

					break;

				case ONED_SEARCH_FIRSTCALL:

					// Returns an error message in case of poor convergence history
					if (stop_because_bad_convergence_ == true)
					{
						nls_status_ = NLS_STOP_FOR_BAD_CONVERGENCE;
						return nls_status_;
					}

					// Stores the merit functions
					v = phiW;
					z = phi1New;

					for(int i = 1;i <= 2;i++)
					{		
						// Calculates the minimum point of the parabola
						// See 25.65, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 561
						if (i == 1)
						{
							u = phiNew / (phi1New + phiNew);
						}
						else
						{
							phiW = v;
							phi1New = z;
							u = phiW/(phi1W + phiW);
						}

						// Better to check if the minimum point of the variable is too small
						// See 25.66, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 561
						if(u < 0.1)
							u = 0.1;

						// 
						fz = phi1W;

						// Calculation of correction to be applied to the current solution
						// and calculation of the new prevision
						// See Point 2, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 509
						aux1 = pi*u;
						xi1 = xi + aux1;

						// Check if the new prevision satisfies the provided contraints (if any)
						// If not, proper corrections are automatically applied
						CheckConstraints(xi1, 1.);

						// Calculates the residuals corresponding to the new prevision
						this->Equations(xi1, fi1);
						number_of_equation_system_calls_++; 
						number_of_equation_system_calls_cumulative_++;
						number_calls_monodimensional_++;

						if(i == 1)
						{
							CalculatesNewCorrectionAndMeritFunction();
							fu = phi1New;
						}
						else
						{			 
							CalculatesMeritFunction(fi1, phi1W);
							fu = phi1W;
						}

						// Check if the prevision can be accepted (true) or not (false)
						control = CheckPrevision(NonLinearSolver<Kernel>::ONE_MINUS_GAMMA_1DSEARCH);

						// The test was successfull (i.e. the prevision can be accepted)
						if(control == true)
						{
							// Checks if the residuals are sufficiently small
							// If yes, this means the algorithm reached the solution
							CheckObjectiveFunctions();

							if(nls_status_ != NLS_CONTINUATION_STATE)
								return nls_status_;

							nls_method_ = QUASI_NEWTON;

							break;
						}
					}
				
					if (control == true)
						break;
				
					// If the 1D search step failed, a second 1D search is performed
					nls_method_ = ONED_SEARCH_SECONDCALL;
				
					break;

				case ONED_SEARCH_SECONDCALL:

					// Returns an error message in case of poor convergence history
					if (stop_because_bad_convergence_ == true)
					{
						nls_status_ = NLS_STOP_FOR_BAD_CONVERGENCE;
						return nls_status_;
					}
				
					// Step 2, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 509
					z = 1.;
					phi1New = 2.*phiNew;
				
					int iSecond;
					for(iSecond = 1;iSecond <= 3; iSecond++)
					{	
						{
							// Calculation of a1 and a2 coefficients
							// See 23.104, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 509
							// Please, consider that F'(0) in  Buzzi Ferraris is equal to -phiW/2
							const double b1 = fu - phiW + 2.*u*phiW;
							const double b2 = fz - phiW + 2.*z*phiW;
							const double bu = b1 / (u*u);
							const double bz = b2 / (z*z);
							const double uz = u - z;
							const double a1 = (-z*bu + u*bz) / uz;
							const double a2 = (bu - bz) / uz;
							double a1a2 = a1*a1 + 6.*a2*phiW;

							// Solve the quadratic equation
							// See 23.105, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 510
							{
								if (a1a2 < 0.)
									a1a2 = 0.;

								if (a2 == 0.)
									v = 0.5*u;
								else if (a1 == 0. && a1a2 == 0.)
									v = 0.1*v;
								else if (a1 < 0.)
									v = (-a1 + std::sqrt(a1a2)) / (3.*a2);
								else
									v = 6.*a2*phiW / (3.*a2*(std::sqrt(a1a2) + a1));

								if (v < 0.1*u)
									v = 0.1*u;

								if (v > 0.5*u)
									v = 0.5*u;
							}
						}

						// Step 6, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 510
						z = u;
						fz = fu;
						u = v;

						// Calculates the provisional correction
						aux1 = pi*u;
						xi1 = xi + aux1;
						
						// Check if the new prevision satisfies the provided contraints (if any)
						// If not, proper corrections are automatically applied
						CheckConstraints(xi1, 1.);

						// Calculates the residuals corresponding to the new prevision
						this->Equations(xi1, fi1);
						number_of_equation_system_calls_++; 
						number_of_equation_system_calls_cumulative_++;
						number_calls_monodimensional_++;

						// Update the merit function
						CalculatesMeritFunction(fi1, phi1W);
						fu = phi1W;		

						// Check if the prevision can be accepted (true) or not (false)
						control = CheckPrevision(1.);
						
						// The test was successfull (i.e. the prevision can be accepted)
						if(control == true)
						{
							// Checks if the residuals are sufficiently small
							// If yes, this means the algorithm reached the solution
							CheckObjectiveFunctions();
							
							if(nls_status_ != NLS_CONTINUATION_STATE)
								return nls_status_;
							
							nls_method_ = GRADIENT;

							break;
						}
					}

					// In case three iterations were not enough to reach a satisfactory prevision, better to move to the gradient method
					if(iSecond > 3)
						nls_method_ = GRADIENT;
				
					break;

				case GRADIENT:
				
					// Returns an error message in case of poor convergence history
					if (stop_because_bad_convergence_ == true)
					{
						nls_status_ = NLS_STOP_FOR_BAD_CONVERGENCE;
						return nls_status_;
					}

					number_calls_gradients_++;

					if(number_calls_gradients_ >= 4)
					{
						if (number_trials_ < 4)
						{
							ReInitialize();
						}
						else
						{
							nls_status_ = NLS_EXCESSIVE_WORK_STATE;
							return nls_status_;
						}
						break;
					}
				
					if(nls_weights_type_ != UNITARY_WEIGHTS)
					{
						for(unsigned int j = 0; j < this->ne_; j++)
							aux(j) = -fi(j)*(weights(j)*weights(j));
					}
					else
					{
						aux = -fi;
					}

					// Calculates the direction for looking for the correction
					// 25.24, Buzzi-Ferraris, Metodi Numerici e Software in C++, p. 548
					// g = J'*f and aux = J*J'*f (please note that originally aux = f)
					this->DoubleProduct(gi, aux);
				
					// Norm of gradient
					const double giNorm2 = gi.norm();

					//
					{
						double an = aux.norm();
						if(an == 0.)
							an = 1.;
						u = giNorm2/an;
					}
					u *= u;
				
					int i;
					if(giNorm2 == 0.)
					{
						i = 1;
						break;
					}

					v = 0.001*xiNorm2/giNorm2;

					if(u < v)
						u = v;
				
					for(i = 1;;i++)
					{
						aux1 = gi*u;
						xi1 = xi + aux1;

						// Check if the new prevision satisfies the provided contraints (if any)
						// If not, proper corrections are automatically applied
						CheckConstraints(xi1, 1.);

						// Calculates the residuals corresponding to the new prevision
						this->Equations(xi1, fi1);
						number_of_equation_system_calls_++; 
						number_of_equation_system_calls_cumulative_++;

						// Update the merit function
						CalculatesMeritFunction(fi1, phi1W);
						phi1New = 1.1*phiNew;

						// Check if the prevision can be accepted (true) or not (false)
						control = CheckPrevision(1.);

						// The test was successfull (i.e. the prevision can be accepted)
						if(control == true)
						{
							// Checks if the residuals are sufficiently small
							// If yes, this means the algorithm reached the solution
							CheckObjectiveFunctions();

							if(nls_status_ != NLS_CONTINUATION_STATE)
								return nls_status_;
						}
						else
							break;
					
						u *= 2.;
					} 

					{	
						int jMax;

						for (unsigned int j = 0; j < this->ne_; j++)
							gi_abs(j) = std::fabs(gi(j));

						double maxGi = gi_abs.maxCoeff(&jMax);
						double minGi = gi_abs.minCoeff();

						if(minGi == 0.)
						{						
							if(maxGi != 0.)
							{ 
								minGi = 0.1*maxGi;
								for (unsigned int k = 0; k < this->ne_; k++)
								{
									if (k != jMax)
										gi(k) = minGi;
								}
							}
							else																		 
							{
								for(unsigned int k = 0;k < this->ne_;k++)
									gi(k) = 1.;														
							}
						}
						if(number_calls_gradients_%3 == 0)
						{
							gi.setConstant(1.);
						}

						for(unsigned int k = 1;k <= std::min(this->ne_, static_cast<unsigned int>(20) );k++)
						{
							for (unsigned int j = 0; j < this->ne_; j++)
								gi_abs(j) = std::fabs(gi(j));

							maxGi = gi_abs.maxCoeff(&jMax);

							u = 0.001*OpenSMOKE::MaxAbs(x_measure_(jMax), xi(jMax));

							gi(jMax) = 0.;
							xi1 = xi;

							int j;
							for(j = 1;;j++)
							{
								xi1(jMax) += u;

								// Check if the new prevision satisfies the provided contraints (if any)
								// If not, proper corrections are automatically applied
								CheckConstraints(xi1, 1.);

								// Calculates the residuals corresponding to the new prevision
								this->Equations(xi1, fi1);
								number_of_equation_system_calls_++; 
								number_of_equation_system_calls_cumulative_++;

								// Update the merit function
								CalculatesMeritFunction(fi1, phi1W);
								phi1New = 1.1*phiNew;

								// Check if the prevision can be accepted (true) or not (false)
								control = CheckPrevision(1.);

								// The test was unsuccessfull (i.e. the prevision cannot be accepted)
								if(control == false)
								{
									if(j == 1)
										u = -u;
									else
										break;					 
								}
								// The test was successfull (i.e. the prevision can be accepted)
								else
								{
									u *= 2.;
								}
							}
				 		
							if(j > 2)
							{
								i = 2; 
								if(number_calls_gradients_ < 3)
									break;
							}
						}
					}
			
					if(i == 1)
					{
						nls_status_ = NLS_EXCESSIVE_WORK_STATE;
						PrintErrorState();
						return nls_status_;
					}
			
					nls_method_ = NEWTON;
			
					break;
				}
			}
		
			nls_status_ = NLS_CONTINUATION_STATE;
	
			return nls_status_;
		}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetAbsoluteTolerances(const Eigen::VectorXd& abs_tolerances)
	{
		if (abs_tolerances.size() != this->ne_)
			OpenSMOKE::FatalErrorMessage("The size of the absolute tolerance vector is wrong");

		abs_tolerances_.resize(this->ne_);
		abs_tolerances_ = abs_tolerances;
		abs_tolerances_scalar_ = false;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetRelativeTolerances(const Eigen::VectorXd& rel_tolerances)
	{
		if (rel_tolerances.size() != this->ne_)
			OpenSMOKE::FatalErrorMessage("The size of the relative tolerance vector is wrong");

		rel_tolerances_.resize(this->ne_);
		rel_tolerances_ = rel_tolerances;
		rel_tolerances_scalar_ = false;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetAbsoluteTolerances(const double abs_tolerance)
	{
		if (abs_tolerance < 1.e-32)
			OpenSMOKE::FatalErrorMessage("The user-defined absolute tolerance is too small");

		abs_tolerance_ = abs_tolerance;
		abs_tolerances_scalar_ = true;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetRelativeTolerances(const double rel_tolerance)
	{
		if (rel_tolerance < 1.e-32)
			OpenSMOKE::FatalErrorMessage("The user-defined relative tolerance is too small");

		rel_tolerance_ = rel_tolerance;
		rel_tolerances_scalar_ = true;
	}


	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetWeights(const Eigen::VectorXd &w)
	{
		if(this->ne_ == w.size())
		{
			nls_weights_type_ = FIXED_WEIGHTS;
			weights = w;
		}
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::SetPrint(const bool flag)
	{
		printResults_ = flag;
	}

	template <typename Kernel>
	void NonLinearSolver<Kernel>::NlsSummary(std::ostream& out)
	{
		if (printResults_ == true)
		{
			out << std::endl;
			out << "NLS system solution" << std::endl;
			out << "---------------------------------------------------------------------------------------" << std::endl;
			out << "* Number of iterations:                     " << number_iterations_ << std::endl;
			out << "* Number of system calls:                   " << number_of_equation_system_calls_ << std::endl;
			out << "* Number of Newton's calls:                 " << number_calls_newtons_ << std::endl;
			out << "* Number of Quasi-Newton's calls:           " << number_calls_quasi_newtons_ << std::endl;
			out << "* Number of 1D searches:                    " << number_calls_monodimensional_ << std::endl;
			out << "* Number of gradient calls:                 " << number_calls_gradients_ << std::endl;
			out << "---------------------------------------------------------------------------------------" << std::endl;
			out << std::endl;

			this->NlsSolverKernelSummary(out);
		}
	}
}
