/*------------------------------------------------------------------------\
|   catalyticFOAM                                                         |
|   http://www.catalyticfoam.polimi.it/                                   |
|                                                                         |
|   Authors:                                                              |
|                                                                         |
|   Alberto Cuoci <alberto.cuoci@polimi.it>                               |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|   Matteo Maestri <matteo.maestri@polimi.it>                             |
|   Department of Energy                                                  |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|   Mauro Bracconi <mauro.bracconi@polimi.it>                             |
|   Department of Energy                                                  |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of catalyticFOAM framework.                         |
|                                                                         |
|   Copyright(C) 2020-2011, A.Cuoci, M.Maestri,                           |
|                2020-2014, M. Bracconi                                   |
|                2015-2013, S.Rebughini                                   |
|                     2013, T.Maffei                                      |
|                     2013, G.Gentile, F.Manelli                          |
|                     2012, M.Calonaci, F.Furnari                         |
|                     2011, S.Goisis, A.Osio                              |
|                                                                         |
|   catalyticFOAM is distributed in the hope that it will be useful,      |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with catalyticFOAM. If not, see <http://www.gnu.org/licenses/>. |
|                                                                         |
\*-----------------------------------------------------------------------*/

	//- Memory allocation: gas-phase properties
	OpenSMOKE::OpenSMOKEVectorDouble massFractions(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble moleFractions(thermodynamicsMapXML->NumberOfSpecies());

	//- Memory allocation: gas-phase chemistry
	OpenSMOKE::OpenSMOKEVectorDouble omega(thermodynamicsMapXML->NumberOfSpecies());
	
	// Batch reactor homogeneous
	BatchReactorHomogeneousConstantPressure batchReactorHomogeneousConstantPressure(*thermodynamicsMapXML, *kineticsMapXML);
	BatchReactorHomogeneousConstantVolume   batchReactorHomogeneousConstantVolume(*thermodynamicsMapXML, *kineticsMapXML);
	#include "createHomogeneousOdeSolvers.H"

	// Batch reactor heterogeneous
	BatchReactorHeterogeneousConstantPressure 	batchReactorHeterogeneousConstantPressure(*thermodynamicsMapXML, *kineticsMapXML, *surfaceThermodynamicsMapXML, *surfaceKineticsMapXML);
	BatchReactorHeterogeneousConstantVolume 	batchReactorHeterogeneousConstantVolume(*thermodynamicsMapXML, *kineticsMapXML, *surfaceThermodynamicsMapXML, *surfaceKineticsMapXML);
	#include "createHeterogeneousOdeSolvers.H"

	#if OPENSMOKE_USE_ISAT == 1
    		#include "memoryAllocation_ISAT.H"
	#endif

