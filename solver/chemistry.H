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

{
	//- Initial conditions
	scalarField& TCells = T.ref();
	const scalarField& pCells = p.ref();
	const scalarField& rhoCells = rho.ref();
	const scalarField& vCells = mesh.V();

	#if OPENSMOKE_USE_ISAT == 1

    		if (isatHomoCheck == true) 
		{
			#include "chemistryHomogeneous_ISAT.H"
		}
		else 
		{	
			#include "chemistryHomogeneous_DI.H"
		}
		

		if (isatHeteroCheck == true)
		{
			#include "chemistryHeterogeneous_ISAT.H"
		}
		else
		{
			#include "chemistryHeterogeneous_DI.H"
		}

	#else

		#include "chemistryHomogeneous_DI.H"
		#include "chemistryHeterogeneous_DI.H"

	#endif

	

	
	// Spark (TOCHECK)
	spark.SetSpark(runTime.value(), T);
	
	// Info
	Info << endl;
    	Info <<" * Chemistry solved with success!"<<endl;
    	Info<< " * T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;
}

				
