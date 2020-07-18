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
								// Set reactor data
								batchReactorHomogeneousConstantPressure.SetReactor(pCells[celli]);
								batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
							
								// Set initial conditions
								odeHomogeneousSolverConstantPressure.SetInitialConditions(t0, y0);

								// Additional ODE solver options
								if (celli == 0)
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
							
								// Solve
								OdeSMOKE::OdeStatus status = odeHomogeneousSolverConstantPressure.Solve(tf);
								odeHomogeneousSolverConstantPressure.Solution(yf);							
							}
							else
							{
								// Set reactor
								batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], pCells[celli], rhoCells[celli]);
								batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
								
								// Set initial conditions
								odeHomogeneousSolverConstantVolume.SetInitialConditions(t0, y0);

								// Additional ODE solver options
								if (celli == 0)
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
								
								// Solve
								OdeSMOKE::OdeStatus status = odeHomogeneousSolverConstantVolume.Solve(tf);
								odeHomogeneousSolverConstantVolume.Solution(yf);
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
			}
		}
		else if (odeParameterBatchReactorHomogeneous.type() != OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE)
		{
			const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
			const unsigned int NEQ = thermodynamicsMapXML->NumberOfSpecies()+1;
			
			// Min and max values
			OpenSMOKE::OpenSMOKEVectorDouble y0(NEQ);
			OpenSMOKE::OpenSMOKEVectorDouble yf(NEQ);

			Info <<" * Solving homogeneous chemistry (external ODE solver)... "<<endl;
			{			
				unsigned int counter = 0;
				
				double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
				forAll(TCells, celli)
				{
					//- Solving for celli
					if (isNeighbourCell[celli] == 0)
					{
						for(unsigned int i=0;i<NC;i++)
							y0[i+1] = Y[i].ref()[celli];
						y0[NC+1] = TCells[celli];
						
						if (constPressureBatchReactor == true)
						{
							// Set reactor pressure
							batchReactorHomogeneousConstantPressure.SetReactor(pCells[celli]);
							batchReactorHomogeneousConstantPressure.SetCheckMassFractions(true);
							batchReactorHomogeneousConstantPressure.SetEnergyEquation(energyEquation);
							
							if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
							{
								#if OPENSMOKE_USE_DVODE == 1
								typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE> Solver;
								SolveOpenSourceSolvers<Solver>(dvodeHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodeHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);					
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodaHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
							{
								#if OPENSMOKE_USE_DASPK == 1
								typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK> Solver;
								SolveOpenSourceSolvers<Solver>(daspkHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);								
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
							{
								#if OPENSMOKE_USE_RADAU == 1
								typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5> Solver;
								SolveOpenSourceSolvers<Solver>(radauHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);						
								#endif							
							}			
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
							{
								#if OPENSMOKE_USE_MEBDF == 1
								typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF> Solver;
								SolveOpenSourceSolvers<Solver>(mebdfHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}	
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
							{
								#if OPENSMOKE_USE_SUNDIALS == 1
								typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE> Solver;
								SolveOpenSourceSolvers<Solver>(cvodeHomogeneousSolverConstantPressure, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}																							
						}
						else
						{
							// Set reactor data
							batchReactorHomogeneousConstantVolume.SetReactor(vCells[celli], pCells[celli], rhoCells[celli]);
							batchReactorHomogeneousConstantVolume.SetCheckMassFractions(true);
							batchReactorHomogeneousConstantVolume.SetEnergyEquation(energyEquation);
							
							if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
							{
								#if OPENSMOKE_USE_DVODE == 1
								typedef OpenSMOKE::OpenSMOKE_DVODE<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE> Solver;
								SolveOpenSourceSolvers<Solver>(dvodeHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODE<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodeHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);					
								#endif
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA)
							{
								#if OPENSMOKE_USE_ODEPACK == 1
								typedef OpenSMOKE::OpenSMOKE_DLSODA<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA> Solver;
								SolveOpenSourceSolvers<Solver>(dlsodaHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
							{
								#if OPENSMOKE_USE_DASPK == 1
								typedef OpenSMOKE::OpenSMOKE_DASPK<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK> Solver;
								SolveOpenSourceSolvers<Solver>(daspkHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);								
								#endif							
							}
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
							{
								#if OPENSMOKE_USE_RADAU == 1
								typedef OpenSMOKE::OpenSMOKE_RADAU<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5> Solver;
								SolveOpenSourceSolvers<Solver>(radauHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);						
								#endif							
							}			
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
							{
								#if OPENSMOKE_USE_MEBDF == 1
								typedef OpenSMOKE::OpenSMOKE_MEBDF<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF> Solver;
								SolveOpenSourceSolvers<Solver>(mebdfHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}	
							else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
							{
								#if OPENSMOKE_USE_SUNDIALS == 1
								typedef OpenSMOKE::OpenSMOKE_CVODE_Sundials<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE> Solver;
								SolveOpenSourceSolvers<Solver>(cvodeHomogeneousSolverConstantVolume, t0, tf, y0, yf, odeParameterBatchReactorHomogeneous);							
								#endif							
							}		
						}
					}
					
					// Check mass fractions
					normalizeMassFractions(yf, celli, massFractionsTol);

					// Assign mass fractions
					for(int i=0;i<NC;i++)
						Y[i].ref()[celli] = yf[i+1];

					//- Allocating final values: temperature
					if (energyEquation == true)
						TCells[celli] = yf[NC+1];

					if (counter%(int(0.20*totalNumberOfNonCatalyticReactors)+1) == 0)
						Info <<"   Accomplished: " << counter << "/" << totalNumberOfNonCatalyticReactors << endl;

					counter++;
				}
				double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
				Info << "   Homogeneous chemistry solved in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(totalNumberOfNonCatalyticReactors)*1000. << " ms per reactor)" << endl;
			}
		}
		
	}
