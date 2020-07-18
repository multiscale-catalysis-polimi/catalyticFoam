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

	Info << endl;
	Info<< " * Reading chemical species mole fractions or mass fractions..." << endl;

	PtrList<volScalarField> Y(thermodynamicsMapXML->NumberOfSpecies());
	PtrList<volScalarField> Dmix(thermodynamicsMapXML->NumberOfSpecies());
	PtrList<volScalarField> HSpecies(thermodynamicsMapXML->NumberOfSpecies());
	PtrList<volScalarField> Ygrad(thermodynamicsMapXML->NumberOfSpecies());
	

	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	{

		Info << "   Reading species: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " > ";

		IOobject header
		(
			thermodynamicsMapXML->NamesOfSpecies()[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		#if OFVERSION > 40
                if (header.typeHeaderOk<volScalarField>(true))
                #else
                if (header.headerOk())
                #endif
		{
			Info << thermodynamicsMapXML->NamesOfSpecies()[i] << endl;
			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);

			Ygrad.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"grad_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, -1, 0, 0, 0)
				)
			);
		}
		else
		{
			Info << "Ydefault" << endl;

			volScalarField Ydefault
			(
				IOobject
				(
					"Ydefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			Y.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					Ydefault
				)
			);

			Ygrad.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"grad_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, 0, -1, 0, 0)
				)
			);
		}

		// Sets of diffusivities and specific heats of species
		{
			Dmix.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"gas::Dmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, 2, -1, 0, 0)
				)
			);

			HSpecies.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
					"gas::H_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
					),
					mesh,
					dimensionSet(0, 2, -2, 0, 0)
				)
			);
		}

	}

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
	//	Viscosity and mass diffusion coefficients initialization
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
	// The fields of viscosity and mass diffusion coefficients are initialized according to
	// the values suggested by the user

	mu = dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 1.8e-5);
	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		Dmix[i] = dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 1.8e-5);

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//

	// Soret effect
	PtrList<volScalarField> Dsoret(transportMapXML->iThermalDiffusionRatios().size());
	if (soretEffect == true)
	{
		for(unsigned int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
		{ 
			unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
			indexSoret[index-1] = i;
			soretEffectList[index-1] = true;

			std::cout << "Soret effect: species " 	<< i+1 << "/" << transportMapXML->iThermalDiffusionRatios().size() 
								<< " " << thermodynamicsMapXML->NamesOfSpecies()[index-1] << " (" << index << ")" << std::endl;

			Dsoret.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"gas::Dsoret_"+ thermodynamicsMapXML->NamesOfSpecies()[index-1],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					mesh,
					dimensionSet(0, 2, -1, 0, 0)
				)
			);
		}
	}

	// Create fluxes fields
	#include "createFluxes.H"

