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
    if(energyEquation == true)
    {

        massDiffusionInEnergyEquation *= 0.;
        if (iMassDiffusionInEnergyEquation == true)
        {
            for (label i=0; i<Y.size(); i++)
                massDiffusionInEnergyEquation -= fvc::div(J[i],HSpecies[i], "div(Ji,hsi)") - HSpecies[i] * fvc::div(J[i]);
        }
    
        #if OFVERSION > 80
        fvScalarMatrix TEqn
        (
            cp*fvm::ddt(rho,T) + 
            cp*fvm::div(phi, T)
            ==
            fvm::laplacian(lambda,T) + 
            massDiffusionInEnergyEquation +  
            DpDt +
            cp*fvModels.source(rho, T)
        );
        #else
        fvScalarMatrix TEqn
        (
            cp*fvm::ddt(rho,T) + 
            cp*fvm::div(phi, T)
            ==
            fvm::laplacian(lambda,T) + 
            massDiffusionInEnergyEquation +  
            DpDt +
            cp*fvOptions(rho, T)
        );
        #endif

        TEqn.relax();
        #if OFVERSION > 80
            fvConstraints.constrain(TEqn);
        #else
            fvOptions.constrain(TEqn);
        #endif
            
        #if OFVERSION > 60
            TEqn.solve("T");
        #else
            TEqn.solve(mesh.solver("T"));
        #endif

        #if OFVERSION > 80
            fvConstraints.constrain(T);
        #else
            fvOptions.correct(T);
        #endif             

    }

    Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

    // Spark
    spark.SetSpark(runTime.value(), T);
}
