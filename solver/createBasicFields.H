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

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Reading pressure [Pa]
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Reading temperature [K]
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Temperature gradient [K/s]
Info<< "Creating field dTdt\n" << endl;
volScalarField grad_T
(
    IOobject
    (
        "grad_T",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet (0, -1, 0, 1, 0)
);

//- Allocating memory for mean molecular weight [kg/mol]
Info<< "Creating field MWmix\n" << endl;
volScalarField MWmix
(
    IOobject
    (
        "gas::MW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 0, 0, 0, -1)
);

//- Allocating memory for total concentration [kmol/m3]
Info<< "Creating field cTot\n" << endl;
volScalarField cTot
(
    IOobject
    (
        "gas::cTot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, -3, 0, 0, 1)
);


//- Allocating memory for p/rho [s2/m2]
Info<< "Creating field psi\n" << endl;
volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, -2, 2, 0, 0)
);

//- Allocating memory for dynamic viscosity [kg/m/s]
Info<< "Creating field mu\n" << endl;
volScalarField mu
(
    IOobject
    (
        "gas::mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -1, 0, 0)
);

//- Allocating memory for thermal conductivity [W/m/K]
Info<< "Creating field lambda\n" << endl;
volScalarField lambda
(
    IOobject
    (
        "gas::lambda",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(1, 1, -3, -1, 0)
);

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field Cp\n" << endl;
volScalarField cp
(
    IOobject
    (
        "gas::Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, -1, 0)
);

//- Allocating memory for mean specific heat [J/kg/K]
Info<< "Creating field Cv\n" << endl;
volScalarField cv
(
    IOobject
    (
        "gas::Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, -1, 0)
);

//- Allocating memory for total enthalpy [m2/s2]
Info<< "Creating field enthalpy\n" << endl;
volScalarField h
(
    IOobject
    (
        "gas::h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0, 2, -2, 0, 0)
);


//- Allocating memory for radiative heat transfer term [W/m3]
Info<< "Creating field Qrad\n" << endl;
volScalarField Qrad
(
    IOobject
    (
        "rad::Qrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rad::Qrad", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

//- Allocating memory for surface heat flux [W/m2]
Info<< "Creating field Qsurf\n" << endl;
volScalarField Qsurf
(
    IOobject
    (
        "rad::Qsurf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rad::Qrad", dimensionSet(1, 2, -3, 0, 0), 0.0)
);

//- Allocating memory for Planck mean absorption coefficient [1/m]
volScalarField aPlanck
(
    IOobject
    (
        "rad::aPlanck",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rad::aPlanck", dimensionSet(0, 0, -1, 0, 0), 0.0)
);

Info<< "Creating massDiffusionInEnergyEquation\n" << endl;
volScalarField massDiffusionInEnergyEquation
(
    IOobject
    (
		"gas::massDiffusionInEnergyEquation",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::NO_WRITE
    ),
    mesh,
    dimensionSet(1, -1, -3, 0, 0)
);
massDiffusionInEnergyEquation = dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.);
grad_T = dimensionedScalar("zero",dimensionSet (0, -1, 0, 1, 0), 0.);
