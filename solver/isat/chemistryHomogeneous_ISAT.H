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

	if (homogeneousReactions == true)
	{
		if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
			
			// Min and max values
			Eigen::VectorXd yMin(NEQ); for(unsigned int i=0;i<NEQ;i++) yMin(i) = 0.;  yMin(NC) = 200.;
			Eigen::VectorXd yMax(NEQ); for(unsigned int i=0;i<NEQ;i++) yMax(i) = 1.;  yMax(NC) = 6000.;
			Eigen::VectorXd y0(NEQ);
			Eigen::VectorXd yf(NEQ);

			Info <<" * Solving homogeneous chemistry... "<<endl;
			{			
				unsigned int counter = 0;
				unsigned int nAddHOM = 0;
				unsigned int nGrowHOM = 0;
				unsigned int nRetHOM = 0;

				double cpuTimeRet      = 0.;
				double cpuTimeDI       = 0.;
				double cpuTimeGrowth   = 0.;
				double cpuTimeAddition = 0.;

				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
				
				forAll(TCells, celli)
				{
					//- Allocating initial conditions
					if (isNeighbourCell[celli] == 0)
					{
						//- Solving for celli:	
						{
							for(unsigned int i=0;i<NC;i++)
								y0(i) = Y[i].ref()[celli];
							y0(NC) = TCells[celli];
							
							if (constPressureBatchReactor == true)
							{
								for(unsigned int i=0;i<NEQ;i++)
									phiISAT_HOM(i) = y0(i)*scalingFactorsHomo_ISAT(i);

								chemComp *phi0base = NULL;
								
								// Retrieve
								if(isat_HOM->retrieve(phiISAT_HOM, phi0base)) 
								{	
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
			
									// makes interpolation
									isat_HOM->interpol(phiISAT_HOM, RphiISAT_HOM, phi0base);
				
									//check negative value
									for(unsigned int i=0;i<NEQ;i++)
										yf(i) = std::max(RphiISAT_HOM(i), 0.)/scalingFactorsHomo_ISAT(i);

									nRetHOM++;	
									
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
									cpuTimeRet += (t2-t1);
								} 
								else 
								{		
									// Direct integration
									{
										double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
	
										// Set reactor data
										batchReactorHomogeneousConstantPressure.SetReactor(pCells[celli]);
										batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
						
										// Set initial conditions
										odeHomogeneousSolverConstantPressure.SetInitialConditions(t0, y0);

										// Additional ODE solver options
										if (nGrowHOM == 0 && nAddHOM == 0)
										{
											// Set linear algebra options
											odeHomogeneousSolverConstantPressure.SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
											odeHomogeneousSolverConstantPressure.SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

											// Set relative and absolute tolerances
											odeHomogeneousSolverConstantPressure.SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
											odeHomogeneousSolverConstantPressure.SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

											// Set minimum and maximum values
											odeHomogeneousSolverConstantPressure.SetMinimumValues(yMin);
											odeHomogeneousSolverConstantPressure.SetMaximumValues(yMax);
										}

										OdeSMOKE::OdeStatus status = odeHomogeneousSolverConstantPressure.Solve(tf);
										odeHomogeneousSolverConstantPressure.Solution(yf);

										// Move the solution from DI to ISAT
										for(unsigned int i=0;i<NEQ;i++)
											RphiISAT_HOM(i) = std::max(yf(i), 0.)*scalingFactorsHomo_ISAT(i);

										double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

										cpuTimeDI += (t2-t1);
									}

									// Growth or Addition
									{
										double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
										if(isat_HOM->grow(phiISAT_HOM, RphiISAT_HOM, phi0base)) 
										{
											nGrowHOM++;
											double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
											cpuTimeGrowth += (t2-t1);
										}
										else
										{
											t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
							
											// compute mapping gradient
											calcMappingGradient(	phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, scalingFactorsHomo_ISAT, 
														luSolverHomo_ISAT, (tf-t0), numberSubStepsHomo_ISAT, &odeHomogeneousSolverConstantPressure);
					
											// add a new leaf 
											bool flag = isat_HOM->add(phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, phi0base); 
				
											if(flag == false)	
												Info << "ISAT Error - Addition process failed..." << endl;
							
											nAddHOM++;
							
											double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
							
											cpuTimeAddition += (t2-t1);
										}
									} 
								}	
			
								// check if is necessary to balance
								isat_HOM->cleanAndBalance();						
							}
							else
							{
								for(unsigned int i=0;i<NEQ;i++)
									phiISAT_HOM(i) = y0(i)*scalingFactorsHomo_ISAT(i);

								chemComp *phi0base = NULL;
								
								// Retrieve
								if(isat_HOM->retrieve(phiISAT_HOM, phi0base)) 
								{	
									double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
			
									// makes interpolation
									isat_HOM->interpol(phiISAT_HOM, RphiISAT_HOM, phi0base);
				
									//check negative value
									for(unsigned int i=0;i<NEQ;i++)
										yf(i) = std::max(RphiISAT_HOM(i), 0.)/scalingFactorsHomo_ISAT(i);

									nRetHOM++;	
									
									double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
									cpuTimeRet += (t2-t1);
								} 
								else 
								{		
									// Direct integration
									{
										double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
	
										// Set reactor
										batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], pCells[celli], rhoCells[celli]);
										batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
								
										// Set initial conditions
										odeHomogeneousSolverConstantVolume.SetInitialConditions(t0, y0);

										// Additional ODE solver options
										if (nGrowHOM == 0 && nAddHOM == 0)
										{
											// Set linear algebra options
											odeHomogeneousSolverConstantVolume.SetLinearAlgebraSolver(odeParameterBatchReactorHomogeneous.linear_algebra());
											odeHomogeneousSolverConstantVolume.SetFullPivoting(odeParameterBatchReactorHomogeneous.full_pivoting());

											// Set relative and absolute tolerances
											odeHomogeneousSolverConstantVolume.SetAbsoluteTolerances(odeParameterBatchReactorHomogeneous.absolute_tolerance());
											odeHomogeneousSolverConstantVolume.SetRelativeTolerances(odeParameterBatchReactorHomogeneous.relative_tolerance());

											// Set minimum and maximum values
											odeHomogeneousSolverConstantVolume.SetMinimumValues(yMin);
											odeHomogeneousSolverConstantVolume.SetMaximumValues(yMax);
										}

										OdeSMOKE::OdeStatus status = odeHomogeneousSolverConstantVolume.Solve(tf);
										odeHomogeneousSolverConstantVolume.Solution(yf);

										// Move the solution from DI to ISAT
										for(unsigned int i=0;i<NEQ;i++)
											RphiISAT_HOM(i) = std::max(yf(i), 0.)*scalingFactorsHomo_ISAT(i);

										double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();

										cpuTimeDI += (t2-t1);
									}

									// Growth or Addition
									{
										double t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
										if(isat_HOM->grow(phiISAT_HOM, RphiISAT_HOM, phi0base)) 
										{
											nGrowHOM++;
											double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
											cpuTimeGrowth += (t2-t1);
										}
										else
										{
											t1 = OpenSMOKE::OpenSMOKEGetCpuTime();
							
											// compute mapping gradient
											calcMappingGradient(	phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, scalingFactorsHomo_ISAT, 
														luSolverHomo_ISAT, (tf-t0), numberSubStepsHomo_ISAT, &odeHomogeneousSolverConstantVolume);
					
											// add a new leaf 
											bool flag = isat_HOM->add(phiISAT_HOM, RphiISAT_HOM, mapGrad_HOM, phi0base); 
				
											if(flag == false)	
												Info << "ISAT Error - Addition process failed..." << endl;
							
											nAddHOM++;
							
											double t2 = OpenSMOKE::OpenSMOKEGetCpuTime();
							
											cpuTimeAddition += (t2-t1);
										}
									} 
								}	
			
								// check if is necessary to balance
								isat_HOM->cleanAndBalance();
							}
						}
		
						// Check mass fractions
						normalizeMassFractions(yf, celli, massFractionsTol);

						// Assign mass fractions
						for(int i=0;i<NC;i++)
							Y[i].ref()[celli] = yf(i);

						//- Allocating final values: temperature
						if (energyEquation == true)
							TCells[celli] = yf(NC);

						if (counter%(int(0.20*totalNumberOfNonCatalyticReactors)+1) == 0)
							Info <<"   Accomplished: " << counter << "/" << totalNumberOfNonCatalyticReactors << endl;

						counter++;
					}
				}
				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfNonCatalyticReactors)*1000. << " ms per reactor)" << endl;

				if(isatHomoCheck == true) 
				{
					Info << endl;
					Info << " ********* ISAT HOM stats **********" << endl;
				
					Info << "   Direct Integration : " << isat_HOM->nAdd()+isat_HOM->nGrow()  << " (" << nAddHOM+nGrowHOM << ")" << " (" << (nAddHOM+nGrowHOM)/double(totalNumberOfNonCatalyticReactors)*100. << "%)" << endl;
					Info << "      Add             : " << isat_HOM->nAdd()  << " (" << nAddHOM  << ")" << " (" << nAddHOM/double(totalNumberOfNonCatalyticReactors)*100. << "%)" << endl;
					Info << "      Grow            : " << isat_HOM->nGrow() << " (" << nGrowHOM << ")" << " (" << nGrowHOM/double(totalNumberOfNonCatalyticReactors)*100. << "%)" << endl;
					Info << "   Retrieve           : " << isat_HOM->nUse()  << " (" << nRetHOM  << ")" << " (" << nRetHOM/double(totalNumberOfNonCatalyticReactors)*100. << "%)" << endl;
					Info << endl;				
		
					const double cpuTimeIntegration = cpuTimeDI + cpuTimeGrowth + cpuTimeAddition;
					Info << "   CPU Integration  : " << cpuTimeIntegration  << " (" << cpuTimeIntegration/(tEnd-tStart+1.e-16)*100. << "%)" << endl;
					Info << "     CPU DI         : " << cpuTimeDI           << " (" << cpuTimeDI/(tEnd-tStart+1.e-16)*100.          << "%)" << endl;				
					Info << "     CPU Growth     : " << cpuTimeGrowth       << " (" << cpuTimeGrowth/(tEnd-tStart+1.e-16)*100.         << "%)" << endl;
					Info << "     CPU Addition   : " << cpuTimeAddition     << " (" << cpuTimeAddition/(tEnd-tStart+1.e-16)*100.    << "%)" << endl;
					Info << "   CPU Retrieve     : " << cpuTimeRet          << " (" << cpuTimeRet/(tEnd-tStart+1.e-16)*100.         << "%)" << endl;
					Info << endl;

					Info << "      BTS  : " << isat_HOM->nBTS()  << endl;
					Info << "      MRU  : " << isat_HOM->nMRU()  << endl;
					Info << "      MFU  : " << isat_HOM->nMFU()  << endl << endl;
					Info << endl;
				}
			}
		}
		else if (odeParameterBatchReactorHomogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			Info << "Fatal Error: in the current version of catalyticPimpleFOAM ISAT can be used only together with the native OpenSMOKE ODE solver." << endl;
			Info << "             External ODE Solvers are not supported." << endl;
			abort();
		}
		
	}
