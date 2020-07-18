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

	bool catalyst = false;

	forAll(Y[0].boundaryField(), patchi)
    	if (Y[0].boundaryFieldRef()[patchi].type() == "catalyticWall")
        {
			catalyst = true;
			break;
		}

	if (catalyst == false)
	{
		Info << "WARNING: Catalyst not found!" << endl;        
		getchar();
	}

	//- Reading adsorbed species
	Info << endl;	
	Info << " * Reading Surface Species:"<<endl;
	PtrList<volScalarField> siteFraction(surfaceThermodynamicsMapXML->number_of_site_species());

	for(int i=0;i<surfaceThermodynamicsMapXML->number_of_site_species();i++)
	{
		std::string name = surfaceThermodynamicsMapXML->NamesOfSpecies()[i+surfaceThermodynamicsMapXML->number_of_gas_species()];
		Info << "   Reading adsorbed specie: " << name << " > ";
	
		IOobject header
		(
			name,
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
			Info << name << " found at time " << mesh.time().timeName() << endl;
				
			siteFraction.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						name,
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
			volScalarField thetaDefault
			(
				IOobject
				(
					"thetaDefault",
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			siteFraction.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						name,
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					thetaDefault
				)
			);

			if(i == 0)
			{
				forAll(siteFraction[i].boundaryField(), patchi)
				{
					if (Y[i].boundaryFieldRef()[patchi].type() == "catalyticWall")
					{
						forAll(siteFraction[i].boundaryFieldRef()[patchi], facei)
							siteFraction[i].boundaryFieldRef()[patchi][facei] = 1.;
					}
					else
					{
						forAll(siteFraction[i].boundaryFieldRef()[patchi], facei)
							siteFraction[i].boundaryFieldRef()[patchi][facei] = 0.;
					}
				}

				Info << "Catalyst's material " <<endl;
			}
			else
			{
				Info << "specie not adsorbed" << endl;
			}
		}
	}
