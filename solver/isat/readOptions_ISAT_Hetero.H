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

Eigen::VectorXd scalingFactorsHetero_ISAT;
Eigen::VectorXd scalingErrorsHetero_ISAT;
ISAT *isat_HET;
int luSolverHetero_ISAT = 1;
label numberSubStepsHetero_ISAT = 1;

Switch isatHeteroCheck(ISATDictionary.subDict("Heterogeneous").lookup("ISAT"));

// Homogeneous
if (isatHeteroCheck == true)
{
	scalar epsilon_ISAT = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("tolerance", 1e-4);
	numberSubStepsHetero_ISAT = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<int>("numberSubSteps", 1);

	scalar maxSizeMRU = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<int>("maxSizeMRU", 100);
	scalar maxSizeMFU = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<int>("maxSizeMFU", 100);
	scalar maxSizeBT = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<int>("maxSizeBT", 100000);
	scalar maxSearchMRU = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<int>("maxSearchMRU", 10);
	scalar maxSearchMFU = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<int>("maxSearchMFU", 30);
	Switch searchMRU(ISATDictionary.subDict("Heterogeneous").lookupOrDefault<Switch>("searchMRU", true)); 
	Switch searchMFU(ISATDictionary.subDict("Heterogeneous").lookupOrDefault<Switch>("searchMFU", true));
	Switch searchBF(ISATDictionary.subDict("Heterogeneous").lookupOrDefault<Switch>("searchBruteForce", false));
	Switch clearIfFull(ISATDictionary.subDict("Heterogeneous").lookupOrDefault<Switch>("clearingIfFull", false));
	Switch cleanAndBalance(ISATDictionary.subDict("Heterogeneous").lookupOrDefault<Switch>("cleanAndBalance", true));
	scalar maxTimeOldCoeff = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("maxTimeOldCoeff", 0.7);
	scalar maxGrowCoeff = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("maxGrowCoeff", 0.5);
	scalar maxHeightCoeff = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("maxHeightCoeff", 20.);
	scalar minUsedCoeff = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("minUsedCoeff", 0.01);
	scalar balanceFactorRetrieve = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("balanceFactorRetrieve", 2.);
	scalar balanceFactorAddition = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<double>("balanceFactorAddition", 0.1);
	word   luFactorization = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<word>("luFactorization","Partial");
	word   qrFactorization = ISATDictionary.subDict("Heterogeneous").lookupOrDefault<word>("qrFactorization","Full");

	
	if (luFactorization != "Partial" && luFactorization != "Full")
	{
		Info << "Wrong luFactorization options: Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (luFactorization == "Full") 
			luSolverHetero_ISAT = 0;
		else
			luSolverHetero_ISAT = 1;		
	}

	int qrSolver_ISAT = 0;
	if (qrFactorization != "NoPivoting" && qrFactorization != "Partial" && qrFactorization != "Full")
	{
		Info << "Wrong qrFactorization options: NoPivoting || Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (qrFactorization == "Full") 
			qrSolver_ISAT = 0;
		else if (qrFactorization == "Partial") 
			qrSolver_ISAT = 1;		
		else if (qrFactorization == "NoPivoting")
			qrSolver_ISAT = 2;
	}

	const unsigned int GAS_NC   = thermodynamicsMapXML->NumberOfSpecies();
	const unsigned int SURF_NC  = surfaceThermodynamicsMapXML->number_of_site_species();
	const unsigned int SURF_NP  = surfaceThermodynamicsMapXML->number_of_site_phases(0);
	const unsigned int NEQ = GAS_NC+1+SURF_NP+SURF_NC+1;
	const unsigned int baseIndexSite = GAS_NC+1+SURF_NP;

	//- ISAT scale factor 
	scalingFactorsHetero_ISAT.setConstant(NEQ,1.);
	{
		dictionary scalingFactors(ISATDictionary.subDict("Heterogeneous").subDict("scalingFactors"));

		// Gas species
		for(unsigned int i=0;i<GAS_NC;i++) 
		{
			if(scalingFactors.found(surfaceThermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingFactorsHetero_ISAT(i) = 1./readScalar(scalingFactors.lookup(surfaceThermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingFactorsHetero_ISAT(i) = 1./readScalar(scalingFactors.lookup("others"));
		}
		
		// Mass
		const double maxVolume  = max(mesh.V()).value();
		const double maxDensity = 2.; 
		scalingFactorsHetero_ISAT(GAS_NC) = 1./scalingFactors.lookupOrDefault("mass",1.0);

		// Gamma
		for(unsigned int i=0;i<SURF_NP;i++)
			scalingFactorsHetero_ISAT(GAS_NC+1+i) = 1./scalingFactors.lookupOrDefault("gamma",1.0);


		// Surface species
		for(unsigned int i=0;i<SURF_NC;i++) 
		{
			if(scalingFactors.found(surfaceThermodynamicsMapXML->vector_names_site_species()[i])) 
				scalingFactorsHetero_ISAT(baseIndexSite+i) = 1./readScalar(scalingFactors.lookup(surfaceThermodynamicsMapXML->vector_names_site_species()[i]));
			else 
				scalingFactorsHetero_ISAT(baseIndexSite+i) = 1./readScalar(scalingFactors.lookup("others"));
		}

		// Temperature
		scalingFactorsHetero_ISAT(NEQ-1)  = 1./readScalar(scalingFactors.lookup("T"));
	}

	//- ISAT scaling error factors
	scalingErrorsHetero_ISAT.setConstant(NEQ,1.);
	{
		dictionary scalingErrors(ISATDictionary.subDict("Heterogeneous").subDict("scalingErrors"));

		// Gas species
		for(unsigned int i=0;i<GAS_NC;i++) 
		{
			if(scalingErrors.found(surfaceThermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingErrorsHetero_ISAT(i) = readScalar(scalingErrors.lookup(surfaceThermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingErrorsHetero_ISAT(i) = readScalar(scalingErrors.lookup("others"));
		}
		
		// Mass
		scalingErrorsHetero_ISAT(GAS_NC) = scalingErrors.lookupOrDefault("mass",1.0);

		// Gamma
		for(unsigned int i=0;i<SURF_NP;i++)
			scalingErrorsHetero_ISAT(GAS_NC+1+i) = scalingErrors.lookupOrDefault("gamma",1.0);

		// Surface species
		for(unsigned int i=0;i<SURF_NC;i++) 
		{
			if(scalingErrors.found(surfaceThermodynamicsMapXML->vector_names_site_species()[i])) 
				scalingErrorsHetero_ISAT(baseIndexSite+i) = readScalar(scalingErrors.lookup(surfaceThermodynamicsMapXML->vector_names_site_species()[i]));
			else 
				scalingErrorsHetero_ISAT(baseIndexSite+i) = readScalar(scalingErrors.lookup("others"));
		}

		// Temperature
		scalingErrorsHetero_ISAT(NEQ-1) = readScalar(scalingErrors.lookup("T"));
	}

	//- ISAT HET
	isat_HET = new ISAT(scalingErrorsHetero_ISAT, epsilon_ISAT, NEQ);	

	// - Setting isat_HET param
	isat_HET->setMaxSizeBT(maxSizeBT);
	isat_HET->setMaxSizeMRU(maxSizeMRU);
	isat_HET->setMaxSizeMFU(maxSizeMFU);
	isat_HET->setMaxSearchMRU(maxSearchMRU);
	isat_HET->setMaxSearchMFU(maxSearchMFU);
	isat_HET->setFlagSearchMRU(searchMRU);
	isat_HET->setFlagSearchMFU(searchMFU);
	isat_HET->setFlagClearingIfFull(clearIfFull);
	isat_HET->setMaxGrowCoeff(maxGrowCoeff);
	isat_HET->setMaxHeightCoeff(maxHeightCoeff);
	isat_HET->setMaxTimeOldCoeff(maxTimeOldCoeff);
	isat_HET->setMinUsedCoeff(minUsedCoeff);
	isat_HET->setBalanceFactorRet(balanceFactorRetrieve);
	isat_HET->setBalanceFactorAdd(balanceFactorAddition);
	isat_HET->setFlagCleanAndBalance(cleanAndBalance);
	isat_HET->setQRType(qrSolver_ISAT);
	isat_HET->setFlagCleanAndBalance(cleanAndBalance);

	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

		Info << endl << "ISAT parameters " << endl;
		Info << "   tolerance           : " << epsilon_ISAT << endl;
		Info << "   luFactorization     : " << luFactorization << endl; 	
		Info << "   qrFactorization     : " << qrFactorization << endl; 	

		Info << "   scalingFactors      : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingFactorsHetero_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingFactorsHetero_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(ISATDictionary.subDict("Heterogeneous").subDict("scalingFactors").lookup("others")) << endl;
		Info << "       mass            : " << scalingFactorsHetero_ISAT(GAS_NC) << endl;
		Info << "       gamma           : " << scalingFactorsHetero_ISAT(GAS_NC+1) << endl;
		Info << "       T               : " << scalingFactorsHetero_ISAT(NEQ-1) << endl;


		Info << "   clear if BT full    : " << clearIfFull << endl;

		Info << "   search in MRU       : " << searchMRU << endl;
		Info << "   search in MFU       : " << searchMFU << endl;
		Info << "   search BF           : " << searchBF << endl;

		Info << "   dimension parameters: " << endl;
		Info << "      max size BT      : " << maxSizeBT << endl;	
		if(searchMRU == true) 
		{
			Info << "      max size MRU     : " << maxSizeMRU << endl;
			Info << "      max search MRU   : " << maxSearchMRU << endl;
		}
		if(searchMFU == true) 
		{
			Info << "      max size MFU     : " << maxSizeMRU << endl;
			Info << "      max search MFU   : " << maxSearchMFU << endl;  
		}

		if (cleanAndBalance == true)
		{
			Info << "   balance parameters  : " << endl;
			Info << "      balanceFactorRetrieve   : " << balanceFactorRetrieve << endl;	
			Info << "      balanceFactorAddition   : " << balanceFactorAddition << endl;	
			Info << "      maxHeightCoeff          : " << maxHeightCoeff << endl;	
			Info << "      maxGrowCoeff            : " << maxGrowCoeff << endl;
			Info << "      minUsedCoeff            : " << minUsedCoeff << endl;
			Info << "      maxTimeOldCoeff         : " << maxTimeOldCoeff << endl; 
		}

		Info << "   scaling error       : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingErrorsHetero_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingErrorsHetero_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(ISATDictionary.subDict("Heterogeneous").subDict("scalingErrors").lookup("others")) << endl;
		Info << "       mass            : " << scalingErrorsHetero_ISAT(GAS_NC) << endl;
		Info << "       gamma           : " << scalingErrorsHetero_ISAT(GAS_NC+1) << endl;
		Info << "       T               : " << scalingErrorsHetero_ISAT(NEQ-1) << endl;

		Info << endl;
	}
}

