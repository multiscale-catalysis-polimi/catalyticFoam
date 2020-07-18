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

#ifndef DaeSolverUtilities_H
#define DaeSolverUtilities_H

namespace DaeSMOKE
{

	enum DaeHStatus
	{
		H_STATUS_DECREASED,
		H_STATUS_CONST,
		H_STATUS_INCREASED
	};

	enum DaeOrderStatus
	{
		ORDER_STATUS_DECREASED,
		ORDER_STATUS_CONST,
		ORDER_STATUS_INCREASED
	};

	enum JacobianType
	{
		JACOBIAN_TYPE_CONST,
		JACOBIAN_TYPE_USERDEFINED,
		JACOBIAN_TYPE_NUMERICAL
	};

	enum JacobianStatus
	{
		JACOBIAN_STATUS_HAS_TO_BE_CHANGED,
		JACOBIAN_STATUS_MODIFIED,
		JACOBIAN_STATUS_OK
	};

	enum FactorizationStatus
	{
		MATRIX_HAS_TO_BE_FACTORIZED,
		MATRIX_FACTORIZED,
	};

	enum DaeConvergence
	{
		CONVERGENCE_STATUS_FAILURE,
		CONVERGENCE_STATUS_OK
	};

	enum DaeStatus
	{
		DAE_STATUS_TO_BE_INITIALIZED = 1,
		DAE_STATUS_CONTINUATION = 2,
		DAE_STATUS_STOP_INTEGRATION_BEFORE_RECALCULATING_JACOBIAN = 10,
		DAE_STATUS_STOP_INTEGRATION_FOR_SMALL_YPRIME_NORM1 = 11,
		DAE_STATUS_MAX_NUMBER_OF_STEPS_REACHED = 12,

		DAE_STATUS_TOO_STRICT_TOLERANCES = -2,		
		DAE_STATUS_ILLEGAL_MAX_INDEPENDENT_VARIABLE = -3,
		DAE_STATUS_MAX_NUMBER_ERRORTEST_FAILURES = -4,
		DAE_STATUS_MAX_NUMBER_CONVERGENCETEST_FAILURES = -5,
		DAE_STATUS_TOO_SMALL_STEP_SIZE = -6,
		DAE_STATUS_YOU_MUST_USE_TCRITIC_STATE = -7,
		DAE_STATUS_ILLEGAL_CONTINUATION_REQUEST = -8,
		DAE_STATUS_ILLEGAL_CONSTRAINTS = -9,
		DAE_STATUS_EXCEPTION_HANDLING_STOP = -10,
		DAE_STATUS_DEINITIALIZE_STATE = -11,
		DAE_STATUS_YOU_CANNOT_OVERSHOOT_TCRITIC = -12,
	};

	unsigned int Factorial(unsigned int n);
}

#include "SolverUtilities.hpp"

#endif // DAESolverUtilities_H