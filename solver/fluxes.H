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
	// Mass diffusion
	forAll (Y,i)
	{
		volScalarField& Yi = Y[i];
		volScalarField& Dmixi = Dmix[i];

		J[i] = ( fvc::interpolate ( -rho*Dmixi ) ) * ( fvc::interpolate ( fvc::grad(Yi) ) & mesh.Sf() );
	}

	// Soret effect
	forAll (Y,i)
	{
		if (soretEffect == true)
		{
			if (soretEffectList[i] == true)
			{
				J[i] += ( fvc::interpolate ( -rho*Dsoret[indexSoret[i]]/T ) ) *
					( fvc::interpolate (  fvc::grad(T) ) & mesh.Sf () ) ;
			}
		}
	}

	// Thermophoretic effect (TODO)
	/*
	forAll (Y,i)
	{
		if (thermophoreticEffect == true)
		{
			if (thermophoreticEffectList[i] == true)
			{
				volScalarField& Yi = Y[i];

				J[i] += ( fvc::interpolate ( -0.55*mu/T*Yi ) ) *
					( fvc::interpolate ( fvc::grad(T) ) & mesh.Sf () ) ;
			}
		}
	}
	*/

	// Correction flux
	Jc *= 0.;
	forAll (Y,i)
	{
		Jc -= J[i];
	}

	// Corrected fluxes
	forAll (Y,i)
	{
		J[i] += Jc*fvc::interpolate(Y[i]);
	}
}
