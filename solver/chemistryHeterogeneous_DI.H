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
    
	if (totalNumberOfCatalyticReactors != 0 && heterogeneousReactions == true)
	{
		if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			Info <<" * Solving heterogeneous chemistry... "<<endl;
		
			int counter = 0;
		
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

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
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
						// Set reactor data
						batchReactorHeterogeneousConstantPressure.SetReactor(pCells[celli], totalCatalyticAreaCell, AlfaCatalyst);
						batchReactorHeterogeneousConstantPressure.SetEnergyEquation(energyEquation);
						batchReactorHeterogeneousConstantPressure.SetReactionHeatFromHeterogeneousReactions(reactionHeatFromHeterogeneousReactions);
						batchReactorHeterogeneousConstantPressure.SetHomogeneousReactions(homogeneousReactions);
				
						// Set initial conditions
						odeHeterogeneousSolverConstantPressure.SetInitialConditions(t0, y0);
	
						// Additional ODE solver options
						if (counter == 0)
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
					}
					else
					{
						// Set reactor data
						batchReactorHeterogeneousConstantVolume.SetReactor(catalyticReactors[counter].volume, totalCatalyticAreaCell, AlfaCatalyst);
						batchReactorHeterogeneousConstantVolume.SetEnergyEquation(energyEquation);
						batchReactorHeterogeneousConstantVolume.SetReactionHeatFromHeterogeneousReactions(reactionHeatFromHeterogeneousReactions);
						batchReactorHeterogeneousConstantVolume.SetHomogeneousReactions(homogeneousReactions);

						// Set initial conditions
						odeHeterogeneousSolverConstantVolume.SetInitialConditions(t0, y0);
	
						// Additional ODE solver options
						if (counter == 0)
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
					}
				
					// Check mass fractions
					normalizeMassAndSiteFractions(yf, celli, NC, SURF_NP, SURF_NC, massFractionsTol, siteFractionsTol);
				
					// Assign mass fractions
					for(int i=0;i<NC;i++)
						Y[i].ref()[celli] = yf(i);
							
					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf(NEQ-1);
				
					//- Allocating finalT values: surface
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

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			Info << "   Heterogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfCatalyticReactors)*1000. << " ms per reactor)" << endl;
		}
		else if (odeParameterBatchReactorHeterogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			Info <<" * Solving heterogeneous chemistry (external ODE solver)... "<<endl;
		
			int counter = 0;
		
			const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int SURF_NC  = surfaceThermodynamicsMapXML->number_of_site_species();
			const unsigned int SURF_NP  = surfaceThermodynamicsMapXML->number_of_site_phases(0);
			const unsigned int NEQ = NC+1+SURF_NP+SURF_NC+1;
			const unsigned int baseIndexSite = NC+1+SURF_NP;
			
			// Min and max values
			OpenSMOKE::OpenSMOKEVectorDouble yMin(NEQ); for (unsigned int i = 1; i<=NEQ; i++) yMin[i] = 0.;  yMin[NEQ] = 200.;
			OpenSMOKE::OpenSMOKEVectorDouble yMax(NEQ); for (unsigned int i = 1; i<=NEQ; i++) yMax[i] = 1.;  yMax[NEQ] = 5000.;
			OpenSMOKE::OpenSMOKEVectorDouble y0(NEQ);
			OpenSMOKE::OpenSMOKEVectorDouble yf(NEQ);

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			forAll(TCells, celli)
			{
			    //- Allocating initial conditions
				if (isNeighbourCell[celli] > 0)
				{		
					// Gas phase mass fractions
					for(unsigned int i=0;i<NC;i++)
						y0[i+1] = Y[i].ref()[celli];
					
					// Total mass of the cell
					y0[NC+1]    = catalyticReactors[counter].volume * rhoCells[celli];
				
					// Phases (TO IMPROVE)
					for(unsigned int i=0;i<SURF_NP;i++)
						y0[NC+2+i] = surfaceThermodynamicsMapXML->matrix_densities_site_phases()[0][i];
					
					// Surface site fractions
					const unsigned int iPatch = catalyticReactors[counter].patchIndex[1];
					const unsigned int iFace = catalyticReactors[counter].faceIndex[1];
					for(unsigned int i=0; i<SURF_NC;i++)
						y0[baseIndexSite+1+i] = siteFraction[i].boundaryFieldRef()[iPatch][iFace];	
				
					// Temperature
					y0[NEQ] = TCells[celli];				

					//- Cycle on number of edges
					double totalCatalyticAreaCell=0.0;
					for (unsigned int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
						totalCatalyticAreaCell += catalyticReactors[counter].faceArea[edge];
					
					if (constPressureBatchReactor == true)
					{

						// Set reactor data
						batchReactorHeterogeneousConstantPressure.SetReactor(pCells[celli], totalCatalyticAreaCell,AlfaCatalyst);
						batchReactorHeterogeneousConstantPressure.SetCheckMassFractions(true);
						batchReactorHeterogeneousConstantPressure.SetEnergyEquation(energyEquation);
						batchReactorHeterogeneousConstantPressure.SetReactionHeatFromHeterogeneousReactions(reactionHeatFromHeterogeneousReactions);
						batchReactorHeterogeneousConstantPressure.SetHomogeneousReactions(homogeneousReactions);

						if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
							{
								#if OPENSMOKE_USE_DVODE == 1
								typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DVODE> Solver;
								SolveOpenSourceSolvers<Solver>(dvodeHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);
								#endif
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DLSODE> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodeHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);					
								#endif
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DLSODA> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodaHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);							
								#endif							
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
							{
								#if OPENSMOKE_USE_DASPK == 1
								typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_DASPK> Solver;
								SolveOpenSourceSolvers<Solver>(daspkHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);								
								#endif							
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
							{
								#if OPENSMOKE_USE_RADAU == 1
								typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_RADAU5> Solver;
								SolveOpenSourceSolvers<Solver>(radauHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);						
								#endif							
							}			
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
							{
								#if OPENSMOKE_USE_MEBDF == 1
								typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_MEBDF> Solver;
								SolveOpenSourceSolvers<Solver>(mebdfHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);							
								#endif							
							}	
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
							{
								#if OPENSMOKE_USE_SUNDIALS == 1
								typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHeterogeneousConstantPressure_ODE_CVODE> Solver;
								SolveOpenSourceSolvers<Solver>(cvodeHeterogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);							
								#endif							
							}	
					}
					else
					{
						
						// Set reactor data
						batchReactorHeterogeneousConstantVolume.SetReactor(catalyticReactors[counter].volume, totalCatalyticAreaCell,AlfaCatalyst);
						batchReactorHeterogeneousConstantPressure.SetCheckMassFractions(true);
						batchReactorHeterogeneousConstantVolume.SetEnergyEquation(energyEquation);
						batchReactorHeterogeneousConstantVolume.SetReactionHeatFromHeterogeneousReactions(reactionHeatFromHeterogeneousReactions);
						batchReactorHeterogeneousConstantVolume.SetHomogeneousReactions(homogeneousReactions);

							if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
							{
								#if OPENSMOKE_USE_DVODE == 1
								typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DVODE> Solver;
								SolveOpenSourceSolvers<Solver>(dvodeHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);
								#endif
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DLSODE> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodeHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);					
								#endif
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DLSODA> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodaHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);							
								#endif							
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
							{
								#if OPENSMOKE_USE_DASPK == 1
								typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_DASPK> Solver;
								SolveOpenSourceSolvers<Solver>(daspkHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);								
								#endif							
							}
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
							{
								#if OPENSMOKE_USE_RADAU == 1
								typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_RADAU5> Solver;
								SolveOpenSourceSolvers<Solver>(radauHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);						
								#endif							
							}			
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
							{
								#if OPENSMOKE_USE_MEBDF == 1
								typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_MEBDF> Solver;
								SolveOpenSourceSolvers<Solver>(mebdfHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);							
								#endif							
							}	
							else if (odeParameterBatchReactorHeterogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
							{
								#if OPENSMOKE_USE_SUNDIALS == 1
								typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHeterogeneousConstantVolume_ODE_CVODE> Solver;
								SolveOpenSourceSolvers<Solver>(cvodeHeterogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHeterogeneous);							
								#endif							
							}
					}
				
					// Check mass fractions
					normalizeMassAndSiteFractions(yf, celli, NC, SURF_NP, SURF_NC, massFractionsTol, siteFractionsTol);
				
					// Assign mass fractions
					for(int i=0;i<NC;i++)
						Y[i].ref()[celli] = yf[i+1];
							
					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf[NEQ];
				
					//- Allocating final values: surface
					for (unsigned int edge=1; edge<=catalyticReactors[counter].numberOfFaces; edge++)
					{
						const unsigned int iPatch = catalyticReactors[counter].patchIndex[edge];
						const unsigned int iFace  = catalyticReactors[counter].faceIndex[edge];

						// temperature
						if (energyEquation == true && reactionHeatFromHeterogeneousReactions == true)
						{
							T.boundaryFieldRef()[iPatch][iFace] = yf[NEQ];
						}

						// mass fractions on boundary
						for(unsigned int i=0;i<NC;i++)
						{
							Y[i].boundaryFieldRef()[iPatch][iFace] = yf[i+1];
						}

						for(unsigned int i=0;i<SURF_NC;i++)
							siteFraction[i].boundaryFieldRef()[iPatch][iFace] = yf[baseIndexSite+1+i];  
					}				

					if (counter%(int(0.20*totalNumberOfCatalyticReactors)+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << totalNumberOfCatalyticReactors << endl;

					counter++;
				}
			
			}

			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
			Info << "   Heterogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfCatalyticReactors)*1000. << " ms per reactor)" << endl;
		}
	}
	else 
	{
		Info <<" * No heterogeneous chemistry (in the master processor)... Wait for the chemistry evaluation in the other processors"<<endl;
	}
	
