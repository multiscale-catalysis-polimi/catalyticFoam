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


    #if OFVERSION > 80
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      - fvm::laplacian(mu, U)
      - fvc::div(mu*dev2(fvc::grad(U)().T()))
         ==
        rho*g
        + fvModels.source(rho, U)
    );
    #else
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      - fvm::laplacian(mu, U)
      - fvc::div(mu*dev2(fvc::grad(U)().T()))
     ==
        rho*g
      + fvOptions(rho, U)
    );
    #endif

    fvVectorMatrix& UEqn = tUEqn.ref();
    
    UEqn.relax();

    #if OFVERSION > 80
        fvConstraints.constrain(UEqn);
    #else
        fvOptions.constrain(UEqn);
    #endif

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        #if OFVERSION > 80
            fvConstraints.constrain(U);
        #else
            fvOptions.correct(U);
        #endif
    }

     Info << " Umin: " << min(U).value() << " Umax: " << max(U).value() << endl;

