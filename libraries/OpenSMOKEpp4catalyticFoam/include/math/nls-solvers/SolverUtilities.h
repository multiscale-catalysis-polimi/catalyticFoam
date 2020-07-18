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

#ifndef NlsSolverUtilities_H
#define NlsSolverUtilities_H

namespace NlsSMOKE
{

	enum NlsStatus
	{
		NLS_INITIALIZATION_STATE = 0,
		NLS_CONTINUATION_STATE = 1,
		NLS_NEWTON_OK_STATE = 2,
		NLS_QUASI_NEWTON_OK_STATE = 3,
		NLS_GRADIENT_OK_STATE = 4,
		NLS_PHINEW_OK_STATE = 5,
		NLS_PHIW_OK_STATE = 6,
		NLS_DUBIOUS_STATE = 7,
		NLS_MAX_NEWTON_CALLS = 8,

		NLS_EXCESSIVE_WORK_STATE = -1,
		NLS_STOP_FOUND = -2,
		NLS_NOT_INITIALIZED = -3,
		NLS_STOP_FOR_BAD_CONVERGENCE = -5
	};


	enum NlsMethod
	{
		NEWTON,
		QUASI_NEWTON,
		ONED_SEARCH_FIRSTCALL,
		ONED_SEARCH_SECONDCALL,
		GRADIENT
	};


	enum NlsWeights
	{
		AUTOMATIC_WEIGHTS,
		FIXED_WEIGHTS,
		UNITARY_WEIGHTS
	};

}

#endif // NlsSolverUtilities_H