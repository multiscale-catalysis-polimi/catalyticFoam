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


// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Additional classes
#include "catalyticReactorClass.H"

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#if OFVERSION > 80
    #include "fvModels.H"
    #include "fvConstraints.H"
#else
    #include "fvOptions.H"
#endif

// Additional include files
#include "sparkModel.H"
#include "userDefinedFunctions.H"

// Homogeneous reactors
#include "BatchReactorHomogeneousConstantPressure.H"
#include "BatchReactorHomogeneousConstantPressure_ODE_Interface.H"
#include "BatchReactorHomogeneousConstantVolume.H"
#include "BatchReactorHomogeneousConstantVolume_ODE_Interface.H"

// Heterogeneous reactors
#include "BatchReactorHeterogeneousConstantPressure.H"
#include "BatchReactorHeterogeneousConstantPressure_ODE_Interface.H"
#include "BatchReactorHeterogeneousConstantVolume.H"
#include "BatchReactorHeterogeneousConstantVolume_ODE_Interface.H"

#include "reactions/reactionRates.H"

// ISAT
#if OPENSMOKE_USE_ISAT == 1
    #include "ISAT.h"
    #include "numericalJacobian4ISAT.H"
    #include "mappingGradient4OpenFOAM.h"
#endif
/*
template<typename Solver, typename OdeBatch>
void SolveOpenSourceSolvers(OdeBatch& ode, const double t0, const double tf, const OpenSMOKE::OpenSMOKEVectorDouble& y0, OpenSMOKE::OpenSMOKEVectorDouble& yf, const OpenSMOKE::ODE_Parameters& parameters)
{
    Solver o(ode);
    o.SetDimensions(y0.Size());
    o.SetAbsoluteTolerance(parameters.absolute_tolerance());
    o.SetRelativeTolerance(parameters.relative_tolerance());
    o.SetAnalyticalJacobian(false);
    o.SetInitialValues(t0, y0.GetHandle());
    o.Solve(tf);
    o.Solution(yf.GetHandle());
}    

*/
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define vector Foam::Vector<double>
    #include "setRootCase.H"
    #undef vector

    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "readGravitationalAcceleration.H"
    
    #include "createBasicFields.H"
    #include "readSolverOptions.H"
    #include "createAdditionalFields.H"
    #include "createCatalyticFields.H"
    #if OFVERSION < 90
        #include "createFvOptions.H"
    #else
        #include "createFvModels.H"
        #include "createFvConstraints.H"
    #endif
    #include "catalystTopology.H"
    #include "memoryAllocation.H"
    #include "properties.H"
    #include "finalSetupOperations.H"
    #include "createMRF.H"
    
    #include "disclaimer.H"
    
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
        
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (strangAlgorithm == STRANG_REACTION_TRANSPORT)
        {
            #include "Policy_ReactionTransport.H"
        }
        else if (strangAlgorithm == STRANG_TRANSPORT_REACTION)
        {
            #include "Policy_TransportReaction.H"
        }
        else if (strangAlgorithm == STRANG_TRANSPORT_REACTION_MOMENTUM)
        {
            #include "Policy_TransportReactionMomentum.H"
        }
        else if (strangAlgorithm == STRANG_REACTION_TRANSPORT_REACTION)
        {
            #include "Policy_ReactionTransportReaction.H"
        }
        else if (strangAlgorithm == STRANG_REACTION_TRANSPORT_HYBRID)
        {
            #include "Policy_ReactionTransportHybrid.H"
        }
                
        runTime.write();

        Info    << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
    }
    
    #include "disclaimer.H"
    Info<< "End\n" << endl;
    
    return 0;
}

// ************************************************************************* //
