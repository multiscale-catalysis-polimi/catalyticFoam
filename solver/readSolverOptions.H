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
	
	Info<< "Reading Solver options\n" << endl;
	IOdictionary solverOptions
	(
		IOobject
		(
			"solverOptions",
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);
	
	// Read the kinetic scheme in XML format

	label inertIndex;
	OpenSMOKE::ThermodynamicsMap_CHEMKIN*			thermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_CHEMKIN* 			kineticsMapXML;
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN*		surfaceThermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped* 		surfaceKineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN* 		transportMapXML;
	
	const dictionary& kineticsDictionary = solverOptions.subDict("Kinetics");
	{
		//- Kinetics folder
		Foam::string kinetics_folder= kineticsDictionary.lookup("folder");
		
		// Reading homogeneous data
		Info<< " * loading homogeneous kinetics...\n" << endl;
		{	
			word kinetics("opensmoke");
			boost::filesystem::path path_kinetics = kinetics_folder;        
		    boost::property_tree::ptree ptree;
		    boost::property_tree::read_xml( (path_kinetics / "kinetics.xml").string(), ptree );

		    double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		    
		    thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree); 
		    transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(ptree); 
		    kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree); 
								
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
		}
		
		// Reading heterogeneous data
		Info<< " * loading heterogeneous kinetics...\n" << endl;
		{	
			word kinetics("opensmoke");
			boost::filesystem::path path_kinetics = kinetics_folder;        
            boost::property_tree::ptree ptree;
            boost::property_tree::read_xml( (path_kinetics / "kinetics.surface.xml").string(), ptree );

            double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
            
            surfaceThermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN(ptree); 
            surfaceKineticsMapXML = new OpenSMOKE::KineticsMap_Surface_CHEMKIN_Lumped(*surfaceThermodynamicsMapXML, ptree); 
								
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
		}	
		
		//- Inert species
		word inertSpecies(kineticsDictionary.lookup("inertSpecies"));
		inertIndex = thermodynamicsMapXML->IndexOfSpecies(inertSpecies)-1;
	}
	
	
	//- Mass fractions tolerance
	scalar massFractionsTol = 1e-8;
	scalar siteFractionsTol = 1e-6;

	//- Internal models
	enum { STRANG_REACTION_TRANSPORT_REACTION, STRANG_REACTION_TRANSPORT, STRANG_TRANSPORT_REACTION, STRANG_REACTION_TRANSPORT_HYBRID, STRANG_TRANSPORT_REACTION_MOMENTUM} strangAlgorithm;
	enum { DIFFUSIVITY_MODEL_MULTICOMPONENT, DIFFUSIVITY_MODEL_LEWIS_NUMBERS} diffusivityModel;

	// Physical model
	Switch energyEquation;
	Switch reactionHeatFromHeterogeneousReactions;
	Switch homogeneousReactions;
	Switch heterogeneousReactions;
	Switch iMassDiffusionInEnergyEquation;
	Switch constPressureBatchReactor;
	Switch soretEffect;
	std::vector<bool> soretEffectList;
	std::vector<int> indexSoret;
	Switch includeDpDt;
	scalar AlfaCatalyst;
	Eigen::VectorXd LewisNumbers;
	List<word> catalyticWalls;
	const dictionary& physicalModelDictionary = solverOptions.subDict("PhysicalModel");
	{
		energyEquation = Switch(physicalModelDictionary.lookup(word("energyEquation")));
		reactionHeatFromHeterogeneousReactions = Switch(physicalModelDictionary.lookup(word("reactionHeatFromHeterogeneousReactions")));
		homogeneousReactions = Switch(physicalModelDictionary.lookup(word("homogeneousReactions")));
		heterogeneousReactions = Switch(physicalModelDictionary.lookup(word("heterogeneousReactions")));
		includeDpDt = Switch(physicalModelDictionary.lookup(word("includeDpDt")));
		iMassDiffusionInEnergyEquation = Switch(physicalModelDictionary.lookup(word("massDiffusionInEnergyEquation")));

		// Soret effect
		soretEffect = Switch(physicalModelDictionary.lookup(word("soretEffect")));
		soretEffectList.resize(thermodynamicsMapXML->NumberOfSpecies());
		indexSoret.resize(thermodynamicsMapXML->NumberOfSpecies());

		for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			soretEffectList[i] = false;	
		for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			indexSoret[i] = -1;
		
		// Algorithm
		{
			word strang(physicalModelDictionary.lookup("strangAlgorithm"));
			     if (strang == "ReactionTransportReaction")	strangAlgorithm = STRANG_REACTION_TRANSPORT_REACTION;
			else if (strang == "ReactionTransport")			strangAlgorithm = STRANG_REACTION_TRANSPORT;
			else if (strang == "TransportReaction")			strangAlgorithm = STRANG_TRANSPORT_REACTION;
			else if (strang == "ReactionTransportHybrid")	strangAlgorithm = STRANG_REACTION_TRANSPORT_HYBRID;
			else if (strang == "TransportReactionMomentum")	strangAlgorithm = STRANG_TRANSPORT_REACTION_MOMENTUM;
			
			else
			{
				Info << "Wrong strangAlgorithm option: ReactionTransportReaction || TransportReaction || ReactionTransport || ReactionTransportHybrid || TransportReactionMomentum" << endl;
				abort();
			}
		}
		
		// Diffusivity of species
		{
			word diffusivity(physicalModelDictionary.lookup("diffusivityModel"));
			if (diffusivity == "multi-component")		diffusivityModel = DIFFUSIVITY_MODEL_MULTICOMPONENT;
			else if (diffusivity == "lewis-numbers")	diffusivityModel = DIFFUSIVITY_MODEL_LEWIS_NUMBERS;
			else
			{
				Info << "Wrong diffusivityModel option: multi-component || lewis-numbers" << endl;
				abort();
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				LewisNumbers.resize(thermodynamicsMapXML->NumberOfSpecies());
				
				const dictionary& lewisNumberDictionary = solverOptions.subDict("LewisNumbers");
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					LewisNumbers(i) = readScalar(lewisNumberDictionary.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
				
				Info << endl;
				Info << "Diffusivities of species will be calculated using the following Lewis numbers: " << endl;
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Info << " " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t" << LewisNumbers(i) << endl;
				Info << endl;
			}
		}
		
		//- Constant pressure reactors
		constPressureBatchReactor = Switch(physicalModelDictionary.lookup(word("constPressureBatchReactor")));
		
		// Catalyst properties
		AlfaCatalyst = readScalar(physicalModelDictionary.lookup("alfaCatalyst"));
		
		
		if(heterogeneousReactions)
		{
			catalyticWalls = List<word>(physicalModelDictionary.lookup("catalyticWalls"));
			
			if ( catalyticWalls.size() == 0 )
			{
				FatalError << "catalytic patch name not defined " << exit(FatalError);
			}
			else
			{
				for ( label i=0;i<catalyticWalls.size();i++)
				{
					label patchID = mesh.boundaryMesh().findPatchID(catalyticWalls[i]);
					if ( patchID == -1 )
					{
						FatalError << catalyticWalls[i] << " not found in mesh " << exit(FatalError);
					}
				}
			}
		}
		
	}

	//- Detect spark
	const dictionary& sparkDictionary = solverOptions.subDict("Spark");
	sparkModel spark;
	spark.Read(sparkDictionary);

	// Batch reactor homogeneous: ode parameters
	const dictionary& odeHomogeneousDictionary = solverOptions.subDict("OdeHomogeneous");
	OpenSMOKE::ODE_Parameters odeParameterBatchReactorHomogeneous;
	{
		//- Mass fractions tolerance
		scalar relTolerance = readScalar(odeHomogeneousDictionary.lookup("relTolerance"));
		scalar absTolerance = readScalar(odeHomogeneousDictionary.lookup("absTolerance"));
		odeParameterBatchReactorHomogeneous.SetRelativeTolerance(relTolerance);
		odeParameterBatchReactorHomogeneous.SetAbsoluteTolerance(absTolerance);
	
		//- Full pivoting (only for OpenSMOKE solver)
		Switch fullPivoting(odeHomogeneousDictionary.lookup("fullPivoting"));
		odeParameterBatchReactorHomogeneous.SetFullPivoting(fullPivoting);
	
		//- Maximum order of integration (only for OpenSMOKE solver)
		label maximumOrder = readLabel(odeHomogeneousDictionary.lookup("maximumOrder"));
		odeParameterBatchReactorHomogeneous.SetMaximumOrder(maximumOrder);
	
		//- 
	
		// Type
		word homogeneousODESolverString(odeHomogeneousDictionary.lookup("odeSolver"));
		if (	homogeneousODESolverString != "OpenSMOKE" 	&& homogeneousODESolverString != "DVODE"  && 
				homogeneousODESolverString != "DLSODE" 		&& homogeneousODESolverString != "DLSODA" && 
				homogeneousODESolverString != "CVODE" 		&& homogeneousODESolverString != "DASPK"  &&
				homogeneousODESolverString != "MEBDF" 		&& homogeneousODESolverString != "RADAU5"  
		   )
		{
			Info << "Wrong homogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
			abort();
		}

		if (homogeneousODESolverString == "OpenSMOKE")	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
		if (homogeneousODESolverString == "DVODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
		if (homogeneousODESolverString == "DLSODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
		if (homogeneousODESolverString == "DLSODA") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
		if (homogeneousODESolverString == "CVODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
		if (homogeneousODESolverString == "DASPK") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
		if (homogeneousODESolverString == "MEBDF") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
		if (homogeneousODESolverString == "RADAU5") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
	
	}

	// Batch reactor heterogeneous: ode parameters
	const dictionary& odeHeterogeneousDictionary = solverOptions.subDict("OdeHeterogeneous");
	OpenSMOKE::ODE_Parameters odeParameterBatchReactorHeterogeneous;
	{
		//- Mass fractions tolerance
		scalar relTolerance = readScalar(odeHeterogeneousDictionary.lookup("relTolerance"));
		scalar absTolerance = readScalar(odeHeterogeneousDictionary.lookup("absTolerance"));
		odeParameterBatchReactorHeterogeneous.SetRelativeTolerance(relTolerance);
		odeParameterBatchReactorHeterogeneous.SetAbsoluteTolerance(absTolerance);
	
		//- Full pivoting (only for OpenSMOKE solver)
		Switch fullPivoting(odeHeterogeneousDictionary.lookup("fullPivoting"));
		odeParameterBatchReactorHeterogeneous.SetFullPivoting(fullPivoting);
	
		//- Maximum order of integration (only for OpenSMOKE solver)
		label maximumOrder = readLabel(odeHeterogeneousDictionary.lookup("maximumOrder"));
		odeParameterBatchReactorHeterogeneous.SetMaximumOrder(maximumOrder);
	
		//- 
	
		// Type
		word heterogeneousODESolverString(odeHeterogeneousDictionary.lookup("odeSolver"));
		if (	heterogeneousODESolverString != "OpenSMOKE" 	&& heterogeneousODESolverString != "DVODE"  && 
				heterogeneousODESolverString != "DLSODE" 		&& heterogeneousODESolverString != "DLSODA" && 
				heterogeneousODESolverString != "CVODE" 		&& heterogeneousODESolverString != "DASPK"  &&
				heterogeneousODESolverString != "MEBDF" 		&& heterogeneousODESolverString != "RADAU5"  
		   )
		{
			Info << "Wrong heterogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
			abort();
		}

		if (heterogeneousODESolverString == "OpenSMOKE")	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
		if (heterogeneousODESolverString == "DVODE") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
		if (heterogeneousODESolverString == "DLSODE") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
		if (heterogeneousODESolverString == "DLSODA") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
		if (heterogeneousODESolverString == "CVODE") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
		if (heterogeneousODESolverString == "DASPK") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
		if (heterogeneousODESolverString == "MEBDF") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
		if (heterogeneousODESolverString == "RADAU5") 	odeParameterBatchReactorHeterogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
	
	}

	#include "checkOdeSolverType.H"

	Info<< "Solver options correctly read\n" << endl;


	#if OPENSMOKE_USE_ISAT == 1
		const dictionary& ISATDictionary = solverOptions.subDict("ISAT");
	    	#include "readOptions_ISAT_Homo.H"
	    	#include "readOptions_ISAT_Hetero.H"
	#endif
