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

// Mass Fluxes
#include "fluxes.H"

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
    
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (i != inertIndex)
        {
            volScalarField& Yi = Y[i];
	    	volScalarField& Dmixi = Dmix[i];

		// Basic equations
		fvScalarMatrix YiEqn
		(
		     fvm::ddt(rho, Yi)
		   + mvConvection->fvmDiv(phi, Yi)
		   - fvm::laplacian(rho*Dmixi, Yi)
		     ==
                   - fvm::div(Jc,Yi, "div(Jc,Yi)") + 						
		     fvOptions(rho, Yi)
		);

		// Add Soret effect
		if (soretEffect == true)
		{ 
			if (soretEffectList[i] == true)
				YiEqn -= fvc::laplacian(rho*Dsoret[indexSoret[i]]/T, T, "laplacian(teta,Yi)");
		}

		// Add thermophoretic effect (TODO)
		// if (thermophoreticEffect == true)
		// {
		//	if (thermophoreticEffectList[i] == true)
		//		YiEqn -= fvc::laplacian(0.55*mu/T*Yi, T, "laplacian(teta,Yi)");
		// }

		// Solve
		YiEqn.relax();
		fvOptions.constrain(YiEqn);
		#if OFVERSION > 60
			YiEqn.solve("Yi");
		#else
			YiEqn.solve(mesh.solver("Yi"));
		#endif
		fvOptions.correct(Yi);	
			
		// Sum of mass fractions
            	Yi.max(0.0);
           	Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
    Info << "   Transport solved in " << tEnd - tStart << " s " << endl;

}
