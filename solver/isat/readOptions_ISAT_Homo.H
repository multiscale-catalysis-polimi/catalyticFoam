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

Eigen::VectorXd scalingFactorsHomo_ISAT;
Eigen::VectorXd scalingErrorsHomo_ISAT;
ISAT *isat_HOM;
int luSolverHomo_ISAT = 1;
label numberSubStepsHomo_ISAT = 1;

Switch isatHomoCheck(ISATDictionary.subDict("Homogeneous").lookup("ISAT"));

// Homogeneous
if (isatHomoCheck == true)
{
	scalar epsilon_ISAT = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("tolerance", 1e-4);
	numberSubStepsHomo_ISAT = ISATDictionary.subDict("Homogeneous").lookupOrDefault<int>("numberSubSteps", 1);

	scalar maxSizeMRU = ISATDictionary.subDict("Homogeneous").lookupOrDefault<int>("maxSizeMRU", 100);
	scalar maxSizeMFU = ISATDictionary.subDict("Homogeneous").lookupOrDefault<int>("maxSizeMFU", 100);
	scalar maxSizeBT = ISATDictionary.subDict("Homogeneous").lookupOrDefault<int>("maxSizeBT", 100000);
	scalar maxSearchMRU = ISATDictionary.subDict("Homogeneous").lookupOrDefault<int>("maxSearchMRU", 10);
	scalar maxSearchMFU = ISATDictionary.subDict("Homogeneous").lookupOrDefault<int>("maxSearchMFU", 30);
	Switch searchMRU(ISATDictionary.subDict("Homogeneous").lookupOrDefault<Switch>("searchMRU", true)); 
	Switch searchMFU(ISATDictionary.subDict("Homogeneous").lookupOrDefault<Switch>("searchMFU", true));
	Switch searchBF(ISATDictionary.subDict("Homogeneous").lookupOrDefault<Switch>("searchBruteForce", false));
	Switch clearIfFull(ISATDictionary.subDict("Homogeneous").lookupOrDefault<Switch>("clearingIfFull", false));
	Switch cleanAndBalance(ISATDictionary.subDict("Homogeneous").lookupOrDefault<Switch>("cleanAndBalance", true));
	scalar maxTimeOldCoeff = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("maxTimeOldCoeff", 0.7);
	scalar maxGrowCoeff = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("maxGrowCoeff", 0.5);
	scalar maxHeightCoeff = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("maxHeightCoeff", 20.);
	scalar minUsedCoeff = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("minUsedCoeff", 0.01);
	scalar balanceFactorRetrieve = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("balanceFactorRetrieve", 2.);
	scalar balanceFactorAddition = ISATDictionary.subDict("Homogeneous").lookupOrDefault<double>("balanceFactorAddition", 0.1);
	word   luFactorization = ISATDictionary.subDict("Homogeneous").lookupOrDefault<word>("luFactorization","Partial");
	word   qrFactorization = ISATDictionary.subDict("Homogeneous").lookupOrDefault<word>("qrFactorization","Full");

	
	if (luFactorization != "Partial" && luFactorization != "Full")
	{
		Info << "Wrong luFactorization options: Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (luFactorization == "Full") 
			luSolverHomo_ISAT = 0;
		else
			luSolverHomo_ISAT = 1;		
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


	//- ISAT scale factor 
	scalingFactorsHomo_ISAT.setConstant(thermodynamicsMapXML->NumberOfSpecies()+1,1.);
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

		massFractionsTol = epsilon_ISAT;

		dictionary scalingFactors(ISATDictionary.subDict("Homogeneous").subDict("scalingFactors"));
		for(unsigned int i=0;i<NC;i++) 
		{	
			if(scalingFactors.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingFactorsHomo_ISAT(i) = 1./readScalar(scalingFactors.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingFactorsHomo_ISAT(i) = 1./readScalar(scalingFactors.lookup("others"));
		}
		scalingFactorsHomo_ISAT(NC) = 1./readScalar(scalingFactors.lookup("T"));
	}

	//- ISAT scaling error factors
	scalingErrorsHomo_ISAT.setConstant(thermodynamicsMapXML->NumberOfSpecies()+1,1.);
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

		dictionary scalingErrors(ISATDictionary.subDict("Homogeneous").subDict("scalingErrors"));

		for(unsigned int i=0;i<NC;i++) 
		{
			if(scalingErrors.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingErrorsHomo_ISAT(i) = readScalar(scalingErrors.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingErrorsHomo_ISAT(i) = readScalar(scalingErrors.lookup("others"));
		}
		scalingErrorsHomo_ISAT(NC) = readScalar(scalingErrors.lookup("T"));
	}

	//- ISAT HOM
	isat_HOM = new ISAT(scalingErrorsHomo_ISAT, epsilon_ISAT, thermodynamicsMapXML->NumberOfSpecies()+1);	

	// - Setting ISAT_HOM param
	isat_HOM->setMaxSizeBT(maxSizeBT);
	isat_HOM->setMaxSizeMRU(maxSizeMRU);
	isat_HOM->setMaxSizeMFU(maxSizeMFU);
	isat_HOM->setMaxSearchMRU(maxSearchMRU);
	isat_HOM->setMaxSearchMFU(maxSearchMFU);
	isat_HOM->setFlagSearchMRU(searchMRU);
	isat_HOM->setFlagSearchMFU(searchMFU);
	isat_HOM->setFlagClearingIfFull(clearIfFull);
	isat_HOM->setMaxGrowCoeff(maxGrowCoeff);
	isat_HOM->setMaxHeightCoeff(maxHeightCoeff);
	isat_HOM->setMaxTimeOldCoeff(maxTimeOldCoeff);
	isat_HOM->setMinUsedCoeff(minUsedCoeff);
	isat_HOM->setBalanceFactorRet(balanceFactorRetrieve);
	isat_HOM->setBalanceFactorAdd(balanceFactorAddition);
	isat_HOM->setFlagCleanAndBalance(cleanAndBalance);
	isat_HOM->setQRType(qrSolver_ISAT);
	isat_HOM->setFlagCleanAndBalance(cleanAndBalance);

	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

		Info << endl << "ISAT parameters " << endl;
		Info << "   tolerance           : " << epsilon_ISAT << endl;
		Info << "   luFactorization     : " << luFactorization << endl; 	
		Info << "   qrFactorization     : " << qrFactorization << endl; 	

		Info << "   scalingFactors      : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingFactorsHomo_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingFactorsHomo_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(ISATDictionary.subDict("Homogeneous").subDict("scalingFactors").lookup("others")) << endl;
		Info << "       T               : " << scalingFactorsHomo_ISAT(NC) << endl;


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
			if(scalingErrorsHomo_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingErrorsHomo_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(ISATDictionary.subDict("Homogeneous").subDict("scalingErrors").lookup("others")) << endl;
		Info << "       T               : " << scalingErrorsHomo_ISAT(NC) << endl;

		Info << endl;
	}
}

