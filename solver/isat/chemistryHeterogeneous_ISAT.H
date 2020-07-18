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
	unsigned int nAddHET  = 0;
	unsigned int nGrowHET = 0;
	unsigned int nRetHET  = 0;

	double cpuTimeRet      = 0.;
	double cpuTimeDI       = 0.;
	double cpuTimeGrowth   = 0.;
	double cpuTimeAddition = 0.; 
	
	double tEnd = 0.;
	double tStart = 0.;
	
	if (totalNumberOfCatalyticReactors != 0 && heterogeneousReactions == true)
	{
		if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			Info <<" * Solving heterogeneous chemistry... "<<endl;
		
			unsigned int counter  = 0;
			
			const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int SURF_NC  = surfaceThermodynamicsMapXML->number_of_site_species();
			const unsigned int SURF_NP  = surfaceThermodynamicsMapXML->number_of_site_phases(0);
			const unsigned int NEQ = NC+1+SURF_NP+SURF_NC+1;
			const unsigned int baseIndexSite = NC+1+SURF_NP;
			
			// Min and max values
			Eigen::VectorXd yMin(NEQ); for (unsigned int i = 0; i<NEQ; i++) yMin(i) = 0.;  yMin(NEQ-1) = 200.;
			Eigen::VectorXd yMax(NEQ); for (unsigned int i = 0; i<NEQ; i++) yMax(i) = 1.;  yMax(NEQ-1) = 5000.;
			Eigen::VectorXd y0(NEQ);
			Eigen::VectorXd yf(NEQ);

			tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
			    //- Allocating initial conditions
				if (isNeighbourCell[celli] > 0)
				{		
					// Gas phase mass fractions
					for(unsigned int i=0;i<NC;i++)
						y0(i) = Y[i].ref()[celli];
					
					// Total mass of the cell
					y0(NC)    = catalyticReactors[counter].volume * rhoCells[celli];
				
					// Phases (TO IMPROVE)
					for(unsigned int i=0;i<SURF_NP;i++)
						y0(NC+1+i) = surfaceThermodynamicsMapXML->matrix_densities_site_phases()[0][i];
					
					// Surface site fractions
					const unsigned int iPatch = catalyticReactors[counter].patchIndex[1];
					const unsigned int iFace = catalyticReactors[counter].faceIndex[1];
					for(unsigned int i=0; i<SURF_NC;i++)
						y0(baseIndexSite+i) = siteFraction[i].boundaryFieldRef()[iPatch][iFace];	
				
					// Temperature
					y0(NEQ-1) = TCells[celli];				

					//- Cycle on number of edges
					double totalCatalyticAreaCell=0.0;
					for (unsigned int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
						totalCatalyticAreaCell += catalyticReactors[counter].faceArea[edge];
					
					if (constPressureBatchReactor == true)
					{
						for(unsigned int i=0;i<NEQ;i++)
							phiISAT_HET(i) = y0(i)*scalingFactorsHetero_ISAT(i);

						chemComp *phi0base = NULL;
						
						// Retrieve
						if(isat_HET->retrieve(phiISAT_HET, phi0base)) 
						{	
							double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
	
							// makes interpolation
							isat_HET->interpol(phiISAT_HET, RphiISAT_HET, phi0base);
		
							//check negative value
							for(unsigned int i=0;i<NEQ;i++)
								yf(i) = std::max(RphiISAT_HET(i), 0.)/scalingFactorsHetero_ISAT(i);

							nRetHET++;	
							
							double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
							cpuTimeRet += (t2-t1);
						} 
						else 
						{
							// Direct integration
							{
								double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();

								// Set reactor data
								batchReactorHeterogeneousConstantPressure.SetReactor(pCells[celli], totalCatalyticAreaCell, AlfaCatalyst);
								batchReactorHeterogeneousConstantPressure.SetEnergyEquation(energyEquation);
								batchReactorHeterogeneousConstantPressure.SetReactionHeatFromHeterogeneousReactions(reactionHeatFromHeterogeneousReactions);
								batchReactorHeterogeneousConstantPressure.SetHomogeneousReactions(homogeneousReactions);
				
								// Set initial conditions
								odeHeterogeneousSolverConstantPressure.SetInitialConditions(t0, y0);
	
								// Additional ODE solver options
								if (nGrowHET == 0 && nAddHET == 0)
								{
									// Set linear algebra options
									odeHeterogeneousSolverConstantPressure.SetLinearAlgebraSolver(odeParameterBatchReactorHeterogeneous.linear_algebra());
									odeHeterogeneousSolverConstantPressure.SetFullPivoting(odeParameterBatchReactorHeterogeneous.full_pivoting());

									// Set relative and absolute tolerances
									odeHeterogeneousSolverConstantPressure.SetAbsoluteTolerances(odeParameterBatchReactorHeterogeneous.absolute_tolerance());
									odeHeterogeneousSolverConstantPressure.SetRelativeTolerances(odeParameterBatchReactorHeterogeneous.relative_tolerance());

									// Set minimum and maximum values
									odeHeterogeneousSolverConstantPressure.SetMinimumValues(yMin);
									odeHeterogeneousSolverConstantPressure.SetMaximumValues(yMax);
								}

								// Solve
								OdeSMOKE::OdeStatus status = odeHeterogeneousSolverConstantPressure.Solve(tf);
								odeHeterogeneousSolverConstantPressure.Solution(yf);	

								// Move the solution from DI to ISAT
								for(unsigned int i=0;i<NEQ;i++)
									RphiISAT_HET(i) = std::max(yf(i), 0.)*scalingFactorsHetero_ISAT(i);

								double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

								cpuTimeDI += (t2-t1);
							}
						
							// Growth or Addition
							{
								double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
								if(isat_HET->grow(phiISAT_HET, RphiISAT_HET, phi0base)) 
								{
									nGrowHET++;
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
									cpuTimeGrowth += (t2-t1);
								}
								else
								{
									t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
					
									// compute mapping gradient
									calcMappingGradient(	phiISAT_HET, RphiISAT_HET, mapGrad_HET, scalingFactorsHetero_ISAT, 
												luSolverHetero_ISAT, (tf-t0), numberSubStepsHetero_ISAT, &odeHeterogeneousSolverConstantPressure);
			
									// add a new leaf 
									bool flag = isat_HET->add(phiISAT_HET, RphiISAT_HET, mapGrad_HET, phi0base); 
									if(flag == false)	
										Info << "ISAT Error - Addition process failed..." << endl;
					
									nAddHET++;
					
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
					
									cpuTimeAddition += (t2-t1);
								}
							} 
						}	
						// check if is necessary to balance
						isat_HET->cleanAndBalance();
					}
					else  // constantVolume
					{

						for(unsigned int i=0;i<NEQ;i++)
							phiISAT_HET(i) = y0(i)*scalingFactorsHetero_ISAT(i);

						chemComp *phi0base = NULL;
						
						// Retrieve
						if(isat_HET->retrieve(phiISAT_HET, phi0base)) 
						{	
							double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
	
							// makes interpolation
							isat_HET->interpol(phiISAT_HET, RphiISAT_HET, phi0base);
		
							//check negative value
							for(unsigned int i=0;i<NEQ;i++)
								yf(i) = std::max(RphiISAT_HET(i), 0.)/scalingFactorsHetero_ISAT(i);

							nRetHET++;	
							
							double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
							cpuTimeRet += (t2-t1);
						} 
						else 
						{
							// Direct integration
							{
								double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();

								// Set reactor data
								batchReactorHeterogeneousConstantVolume.SetReactor(catalyticReactors[counter].volume, totalCatalyticAreaCell, AlfaCatalyst);
								batchReactorHeterogeneousConstantVolume.SetEnergyEquation(energyEquation);
								batchReactorHeterogeneousConstantVolume.SetReactionHeatFromHeterogeneousReactions(reactionHeatFromHeterogeneousReactions);
								batchReactorHeterogeneousConstantVolume.SetHomogeneousReactions(homogeneousReactions);
				
								// Set initial conditions
								odeHeterogeneousSolverConstantVolume.SetInitialConditions(t0, y0);
	
								// Additional ODE solver options
								if (nGrowHET == 0 && nAddHET == 0)
								{
									// Set linear algebra options
									odeHeterogeneousSolverConstantVolume.SetLinearAlgebraSolver(odeParameterBatchReactorHeterogeneous.linear_algebra());
									odeHeterogeneousSolverConstantVolume.SetFullPivoting(odeParameterBatchReactorHeterogeneous.full_pivoting());

									// Set relative and absolute tolerances
									odeHeterogeneousSolverConstantVolume.SetAbsoluteTolerances(odeParameterBatchReactorHeterogeneous.absolute_tolerance());
									odeHeterogeneousSolverConstantVolume.SetRelativeTolerances(odeParameterBatchReactorHeterogeneous.relative_tolerance());

									// Set minimum and maximum values
									odeHeterogeneousSolverConstantVolume.SetMinimumValues(yMin);
									odeHeterogeneousSolverConstantVolume.SetMaximumValues(yMax);
								}

								// Solve
								OdeSMOKE::OdeStatus status = odeHeterogeneousSolverConstantVolume.Solve(tf);
								odeHeterogeneousSolverConstantVolume.Solution(yf);		

								// Move the solution from DI to ISAT
								for(unsigned int i=0;i<NEQ;i++)
									RphiISAT_HET(i) = std::max(yf(i), 0.)*scalingFactorsHetero_ISAT(i);

								double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

								cpuTimeDI += (t2-t1);
							}
						
							// Growth or Addition
							{
								double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
								if(isat_HET->grow(phiISAT_HET, RphiISAT_HET, phi0base)) 
								{
									nGrowHET++;
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
									cpuTimeGrowth += (t2-t1);
								}
								else
								{
									t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
					
									// compute mapping gradient
									calcMappingGradient(	phiISAT_HET, RphiISAT_HET, mapGrad_HET, scalingFactorsHetero_ISAT, 
												luSolverHetero_ISAT, (tf-t0), numberSubStepsHetero_ISAT, &odeHeterogeneousSolverConstantVolume);
			
									// add a new leaf 
									bool flag = isat_HET->add(phiISAT_HET, RphiISAT_HET, mapGrad_HET, phi0base); 
		
									if(flag == false)	
										Info << "ISAT Error - Addition process failed..." << endl;
					
									nAddHET++;
					
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
					
									cpuTimeAddition += (t2-t1);
								}
							} 
						}	
			
						// check if is necessary to balance
						isat_HET->cleanAndBalance();	
					}
					
					// Check mass fractions
					normalizeMassAndSiteFractions(yf, celli, NC, SURF_NP, SURF_NC, massFractionsTol, siteFractionsTol);
				
					// Assign mass fractions
					for(int i=0;i<NC;i++)
						Y[i].ref()[celli] = yf(i);
							
					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf(NEQ-1);
				
					//- Allocating final values: surface
					for (unsigned int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
					{
						const unsigned int iPatch = catalyticReactors[counter].patchIndex[edge];
						const unsigned int iFace  = catalyticReactors[counter].faceIndex[edge];

						// temperature
						if (energyEquation == true && reactionHeatFromHeterogeneousReactions == true)
						{
							T.boundaryFieldRef()[iPatch][iFace] = yf(NEQ-1);
						grad_T.boundaryFieldRef()[iPatch][iFace] = 0.0;
						}

						// mass fractions on boundary
						for(unsigned int i=0;i<NC;i++)
						{
							Y[i].boundaryFieldRef()[iPatch][iFace] = yf(i);
							Ygrad[i].boundaryFieldRef()[iPatch][iFace] = 0.0;
						}

						for(unsigned int i=0;i<SURF_NC;i++)
							siteFraction[i].boundaryFieldRef()[iPatch][iFace] = yf(baseIndexSite+i);  
					}				

					if (counter%(int(0.20*totalNumberOfCatalyticReactors)+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << totalNumberOfCatalyticReactors << endl;

					counter++;
				}
			
			}

			tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			Info << "   Heterogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfCatalyticReactors)*1000. << " ms per reactor)" << endl;
			
			const double cpuTimeIntegration = cpuTimeDI + cpuTimeGrowth + cpuTimeAddition;
		
			if(isatHeteroCheck == true && !Pstream::parRun()) 
			{
				Info << endl;
				Info << " ********* ISAT HET stats **********" << endl;
		
				Info << "   Direct Integration : " << isat_HET->nAdd()+isat_HET->nGrow()  << " (" << nAddHET+nGrowHET << ")" << " (" << (nAddHET+nGrowHET)/double(totalNumberOfCatalyticReactors)*100. << "%)" << endl;
				Info << "      Add             : " << isat_HET->nAdd()  << " (" << nAddHET  << ")" << " (" << nAddHET/double(totalNumberOfCatalyticReactors)*100. << "%)" << endl;
				Info << "      Grow            : " << isat_HET->nGrow() << " (" << nGrowHET << ")" << " (" << nGrowHET/double(totalNumberOfCatalyticReactors)*100. << "%)" << endl;
				Info << "   Retrieve           : " << isat_HET->nUse()  << " (" << nRetHET  << ")" << " (" << nRetHET/double(totalNumberOfCatalyticReactors)*100. << "%)" << endl;
				Info << endl;				

				Info << "   CPU Integration  : " << cpuTimeIntegration  << " (" << cpuTimeIntegration/(tEnd-tStart+1.e-16)*100. << "%)" << endl;
				Info << "     CPU DI         : " << cpuTimeDI           << " (" << cpuTimeDI/(tEnd-tStart+1.e-16)*100.          << "%)" << endl;				
				Info << "     CPU Growth     : " << cpuTimeGrowth       << " (" << cpuTimeGrowth/(tEnd-tStart+1.e-16)*100.         << "%)" << endl;
				Info << "     CPU Addition   : " << cpuTimeAddition     << " (" << cpuTimeAddition/(tEnd-tStart+1.e-16)*100.    << "%)" << endl;
				Info << "   CPU Retrieve     : " << cpuTimeRet          << " (" << cpuTimeRet/(tEnd-tStart+1.e-16)*100.         << "%)" << endl;
				Info << endl;

				Info << " ********* ISAT HET Retrieve stats **********" << endl;
				Info << "      BTS  : " << isat_HET->nBTS()  << endl;
				Info << "      MRU  : " << isat_HET->nMRU()  << endl;
				Info << "      MFU  : " << isat_HET->nMFU()  << endl << endl;
				Info << endl;
			}
		}
		else if (odeParameterBatchReactorHeterogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			Info << "Fatal Error: in the current version of catalyticPimpleFOAM ISAT can be used only together with the native OpenSMOKE ODE solver." << endl;
			Info << "             External ODE Solvers are not supported." << endl;
			abort();
		}
	}
	else 
	{
		Info <<" * No heterogeneous chemistry in the master processor... "<<endl;
	}

	bool checkheterogeneousReactions = heterogeneousReactions;
	label tempTotalNumber = totalNumberOfCatalyticReactors;

	reduce(tempTotalNumber, sumOp<label>());
	reduce(checkheterogeneousReactions, maxOp<bool>());
		
	if (tempTotalNumber != 0 && checkheterogeneousReactions == true && Pstream::parRun())
	{
		unsigned int isatnAdd = isat_HET->nAdd();
		unsigned int isatnGrow = isat_HET->nGrow();
		unsigned int isatnUse = isat_HET->nUse();
		unsigned int isatnBTS = isat_HET->nBTS();
		unsigned int isatnMRU = isat_HET->nMRU();
		unsigned int isatnMFU = isat_HET->nMFU();
		double interval = tEnd-tStart;
	
		reduce(isatnAdd, sumOp<unsigned int>());
		reduce(isatnGrow, sumOp<unsigned int>());
		reduce(isatnUse, sumOp<unsigned int>());
		reduce(isatnBTS, sumOp<unsigned int>());
		reduce(isatnMRU, sumOp<unsigned int>());
		reduce(isatnMFU, sumOp<unsigned int>());
		reduce(interval, maxOp<double>());
		reduce(nAddHET, sumOp<unsigned int>());
		reduce(nGrowHET, sumOp<unsigned int>());
		reduce(nRetHET, sumOp<unsigned int>());
		reduce(cpuTimeDI, maxOp<double>());
		reduce(cpuTimeGrowth, maxOp<double>());
		reduce(cpuTimeAddition, maxOp<double>());
		
		const double cpuTimeIntegration = cpuTimeDI + cpuTimeGrowth + cpuTimeAddition;
		
		if(isatHeteroCheck == true) 
		{
			Info << endl;
			Info << " ********* ISAT HET stats **********" << endl;
		
			Info << "   Direct Integration : " << isatnAdd+isatnGrow  << " (" << nAddHET+nGrowHET << ")" << " (" << (nAddHET+nGrowHET)/double(tempTotalNumber)*100. << "%)" << endl;
			Info << "      Add             : " << isatnAdd  << " (" << nAddHET  << ")" << " (" << nAddHET/double(tempTotalNumber)*100. << "%)" << endl;
			Info << "      Grow            : " << isatnGrow << " (" << nGrowHET << ")" << " (" << nGrowHET/double(tempTotalNumber)*100. << "%)" << endl;
			Info << "   Retrieve           : " << isatnUse  << " (" << nRetHET  << ")" << " (" << nRetHET/double(tempTotalNumber)*100. << "%)" << endl;
			Info << endl;				

			Info << "   CPU Integration  : " << cpuTimeIntegration  << " (" << cpuTimeIntegration/(interval+1.e-16)*100. << "%)" << endl;
			Info << "     CPU DI         : " << cpuTimeDI           << " (" << cpuTimeDI/(interval+1.e-16)*100.          << "%)" << endl;				
			Info << "     CPU Growth     : " << cpuTimeGrowth       << " (" << cpuTimeGrowth/(interval+1.e-16)*100.         << "%)" << endl;
			Info << "     CPU Addition   : " << cpuTimeAddition     << " (" << cpuTimeAddition/(interval+1.e-16)*100.    << "%)" << endl;
			Info << "   CPU Retrieve     : " << cpuTimeRet          << " (" << cpuTimeRet/(interval+1.e-16)*100.         << "%)" << endl;
			Info << endl;

			Info << " ********* ISAT HET Retrieve stats **********" << endl;
			Info << "      BTS  : " << isatnBTS  << endl;
			Info << "      MRU  : " << isatnMRU  << endl;
			Info << "      MFU  : " << isatnMFU  << endl << endl;
			Info << endl;
		}
	}
}
