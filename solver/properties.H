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
	Info<< "Properties evaluation..." << endl;

	OpenSMOKE::OpenSMOKEVectorDouble HVector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble Dmixvector(thermodynamicsMapXML->NumberOfSpecies());
	OpenSMOKE::OpenSMOKEVectorDouble tetamixvector(thermodynamicsMapXML->NumberOfSpecies());

	// Internal fields
	{
		const scalarField& TCells = T.ref();
		const scalarField& pCells = p.ref(); 

		scalarField& cTotCells = cTot.ref();
		scalarField&  psiCells =  psi.ref();
		scalarField&  muCells  =  mu.ref();
		scalarField& lambdaCells = lambda.ref();
		scalarField& cpCells = cp.ref();
		scalarField& hCells = h.ref();
		scalarField& MWmixCells = MWmix.ref();

		forAll(TCells, celli)
		{
			thermodynamicsMapXML->SetPressure(pCells[celli]);
			thermodynamicsMapXML->SetTemperature(TCells[celli]);
			
			transportMapXML->SetPressure(pCells[celli]);
			transportMapXML->SetTemperature(TCells[celli]);
				
	
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].ref()[celli];
				
			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),MWmixCells[celli],massFractions.GetHandle());

			cTotCells[celli] = pCells[celli]/PhysicalConstants::R_J_kmol/TCells[celli];
			psiCells[celli]  = cTotCells[celli]*MWmixCells[celli]/pCells[celli];
			hCells[celli] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());		// [J/kmol]
			hCells[celli] /= MWmixCells[celli];															// [J/kg]
			
			muCells[celli] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());

            if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            {
				lambdaCells[celli] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
                
                		cpCells[celli] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
				cpCells[celli] = cpCells[celli]/MWmixCells[celli];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->hMolar_Species(HVector.GetHandle());
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
						HSpecies[i].ref()[celli] = HVector[i+1] / thermodynamicsMapXML->MW(i);
				}
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(), moleFractions.GetHandle());
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].ref()[celli] = Dmixvector[i+1];
			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient = lambdaCells[celli]/(pCells[celli]*psiCells[celli])/cpCells[celli];
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].ref()[celli] = coefficient/LewisNumbers(i);
			}	

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());		
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
					Dsoret[i].ref()[celli] =Dmix[index-1].ref()[celli]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/MWmixCells[celli];
				} 
			}
		}
    }

	// Boundaries
	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pp = p.boundaryFieldRef()[patchi];
		fvPatchScalarField& pT = T.boundaryFieldRef()[patchi];
        
		fvPatchScalarField& pcTot = cTot.boundaryFieldRef()[patchi];
		fvPatchScalarField& ppsi  = psi.boundaryFieldRef()[patchi];
		fvPatchScalarField& pmu   = mu.boundaryFieldRef()[patchi];
        	fvPatchScalarField& plambda = lambda.boundaryFieldRef()[patchi];
       		fvPatchScalarField& pcp = cp.boundaryFieldRef()[patchi];
        	fvPatchScalarField& ph = h.boundaryFieldRef()[patchi];
        	fvPatchScalarField& pMWmix = MWmix.boundaryFieldRef()[patchi];
        
		forAll(pT, facei)
		{
			thermodynamicsMapXML->SetPressure(pp[facei]);
			thermodynamicsMapXML->SetTemperature(pT[facei]);

			transportMapXML->SetPressure(pp[facei]);
			transportMapXML->SetTemperature(pT[facei]);
			
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				massFractions[i+1] = Y[i].boundaryFieldRef()[patchi][facei];

			thermodynamicsMapXML->MoleFractions_From_MassFractions(moleFractions.GetHandle(),pMWmix[facei],massFractions.GetHandle());


			pcTot[facei] = pp[facei]/(PhysicalConstants::R_J_kmol*pT[facei]);
			ppsi[facei]  = pcTot[facei]*pMWmix[facei]/pp[facei];
			ph[facei] = thermodynamicsMapXML->hMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());		// [J/kmol]
			ph[facei] /= pMWmix[facei];																// [J/kg]
			
			pmu[facei] = transportMapXML->DynamicViscosity(moleFractions.GetHandle());

            if (energyEquation == true || diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
            {
				plambda[facei] = transportMapXML->ThermalConductivity(moleFractions.GetHandle());
				
				pcp[facei] = thermodynamicsMapXML->cpMolar_Mixture_From_MoleFractions(moleFractions.GetHandle());			//[J/Kmol/K]
				pcp[facei] = pcp[facei]/pMWmix[facei];
			
				if (iMassDiffusionInEnergyEquation == true)
				{
					thermodynamicsMapXML->hMolar_Species(HVector.GetHandle());
					for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
						HSpecies[i].boundaryFieldRef()[patchi][facei] = HVector[i+1] / thermodynamicsMapXML->MW(i);
				}
			}
			
			if (diffusivityModel == DIFFUSIVITY_MODEL_MULTICOMPONENT)
			{
				transportMapXML->MassDiffusionCoefficients(Dmixvector.GetHandle(),moleFractions.GetHandle());
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryFieldRef()[patchi][facei] = Dmixvector[i+1];
			}
			else if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
			{
				const double coefficient =  plambda[facei]/(pp[facei]*ppsi[facei])/pcp[facei];
				for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
					Dmix[i].boundaryFieldRef()[patchi][facei] = coefficient/LewisNumbers(i);
			}	

			// Thermal diffusion coefficients [-]
			if (soretEffect == true)
			{
				transportMapXML->ThermalDiffusionRatios(tetamixvector.GetHandle(), moleFractions.GetHandle());
				for(int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
				{
					unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
					Dsoret[i].boundaryFieldRef()[patchi][facei] = Dmix[index-1].boundaryFieldRef()[patchi][facei]*tetamixvector[index]*thermodynamicsMapXML->MW(index-1)/pMWmix[facei];
				}
			}		

            for(int i=0; i<indexCatalyst.size(); i++)
            {
                if (patchi == indexCatalyst[i])
                {
                    for(label i=0; i<Ygrad.size(); i++)
						Ygrad[i].boundaryFieldRef()[patchi][facei] *= thermodynamicsMapXML->MW(i)*AlfaCatalyst/Dmix[i].boundaryFieldRef()[patchi][facei]/(pp[facei]*ppsi[facei]); 
			
                    if(energyEquation == true)
                    {
                     	  grad_T.boundaryFieldRef()[patchi][facei] *= AlfaCatalyst/plambda[facei];
                    }
                }
            }

		}
	}
}
