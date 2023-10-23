/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "math/OpenSMOKEUtilities.h"
#include "ThermodynamicsMap.h"

namespace OpenSMOKE
{
	KineticsMap_Surface_CHEMKIN::KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN& thermo, const unsigned int nSpecies) :
	thermodynamics_(thermo)
	{
		this->number_of_species_ = nSpecies;
		this->T_ = this->P_ = 0.;
	}

	KineticsMap_Surface_CHEMKIN::KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN& thermo, boost::property_tree::ptree& ptree) :
	thermodynamics_(thermo)
	{
		ImportSpeciesFromXMLFile(ptree);
		ImportCoefficientsFromXMLFile(ptree);
		this->T_ = this->P_ = 0.;
	}

	void KineticsMap_Surface_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_old_ = this->T_;
		this->T_ = T;
		this->uT_ = 1./this->T_;
		this->logT_ = log(this->T_);
		Patm_over_RT_ = 101325./PhysicalConstants::R_J_kmol/this->T_;
		log_Patm_over_RT_ = log(Patm_over_RT_);

		if (std::fabs(this->T_-this->T_old_)/this->T_>1.e-14)
		{
			arrhenius_kinetic_constants_must_be_recalculated_ = true;
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
			reaction_h_and_s_must_be_recalculated_ = true;
		}
	}

	void KineticsMap_Surface_CHEMKIN::SetPressure(const double& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;
	//	if (std::fabs(this->P_-this->P_old_)/this->P_>1.e-14)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = true;
		}
	}
	 
	void KineticsMap_Surface_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)
	{
		std::string kinetics_type =  ptree.get<std::string>("opensmoke.Kinetics.<xmlattr>.type");
		std::string kinetics_version = ptree.get<std::string>("opensmoke.Kinetics.<xmlattr>.version");
		if (kinetics_type != "OpenSMOKE" || kinetics_version != "01-02-2014")
			ErrorMessage("void KineticsMap_Surface_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The current kinetic scheme is not supported.");

		unsigned int target_material = 1;
		BOOST_FOREACH( boost::property_tree::ptree::value_type const& node, ptree.get_child( "opensmoke.Kinetics" ) ) 
		{
			boost::property_tree::ptree subtree = node.second;  

			if( node.first == "MaterialKinetics" ) 
			{
				if (target_material == subtree.get<unsigned int>("<xmlattr>.index") )
				{
					this->number_of_species_ = subtree.get<unsigned int>("NumberOfSpecies"); 
					this->number_of_reactions_ = subtree.get<unsigned int>("NumberOfReactions"); 

					// Irreversible reactions
					{
						std::cout << "Reading irreversible..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Irreversible") );  
						Load(indices_of_irreversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_irreversible_reactions_ = indices_of_irreversible_reactions__.size();
					}
		
					// Reversible reactions
					{
						std::cout << "Reading reversible..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Reversible") );  
						Load(indices_of_reversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_reversible_reactions_ = indices_of_reversible_reactions__.size();
					}

					// Thermodynamic Reversible reactions
					{
						std::cout << "Reading reversible thermodynamics..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Reversible-Thermodynamics") );  
						Load(indices_of_thermodynamic_reversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_thermodynamic_reversible_reactions_ = indices_of_thermodynamic_reversible_reactions__.size();
					}

					// Explicit Reversible reactions
					{
						std::cout << "Reading reversible explicit..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Reversible-Explicit") );  
						Load(indices_of_explicitly_reversible_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_explicitly_reversible_reactions_ = indices_of_explicitly_reversible_reactions__.size();
					}

					// Stick reactions
					{
						std::cout << "Reading stick reactions..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Stick") );  
						Load(indices_of_stick_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_stick_reactions_ = indices_of_stick_reactions__.size();
					}

					// Cov dependent reactions
					{
						std::cout << "Reading coverage dependent reactions..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("CoverageDependent") );  
						Load(indices_of_coverage_dependent_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_coverage_dependent_reactions_ = indices_of_coverage_dependent_reactions__.size();
					}

					// Langmuir-Hinshelwood reactions
					{
						std::cout << "Reading Langmuir-Hinshelwood reactions..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Langmuir") );  
						Load(indices_of_langmuir_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_langmuir_reactions_ = indices_of_langmuir_reactions__.size();
					}

					// Lumped reactions
					{
						std::cout << "Reading Lumped reactions..." << std::endl;
						std::stringstream stream;
						stream.str( subtree.get< std::string >("Lumped") );  
						Load(indices_of_lumped_reactions__, stream, OPENSMOKE_FORMATTED_FILE);
						number_of_lumped_reactions_ = indices_of_lumped_reactions__.size();
					}

					// Reading if the kinetic scheme is conventional or UBI-QEP
					{
						std::cout << "Reading type of kinetic scheme..." << std::endl;
						std::string dummy = subtree.get<std::string>("TypeOfKinetics"); 
						dummy.erase(std::remove(dummy.begin(), dummy.end(), '\n'), dummy.end());
					
						if (dummy == "UBIQEP")				type_of_kinetics_ = TYPE_OF_KINETICS_UBI_QEP;
						else if (dummy == "chemkin_conventional")	type_of_kinetics_ = TYPE_OF_KINETICS_CHEMKIN_CONVENTIONAL;
						else
						{
							std::cout << "This type of kinetic mechanism is not available: *" << dummy << "*" << std::endl;
							std::cout << "Available types: chemkin_conventional | UBIQEP" << std::endl;
							std::cout << "Press enter to exit..." << std::endl;						
							getchar();
							exit(OPENSMOKE_FATAL_ERROR_EXIT);
						}
					}

				
					// Reading UBI
					if (type_of_kinetics_ == TYPE_OF_KINETICS_UBI_QEP)
					{
						ubiqep_submechanism_ = new UBIQEP_SubMechanism();
						ubiqep_submechanism_->ReadFromXMLFile(subtree, thermodynamics_.MWs());
					}

					std::cout << " * Reading kinetic parameters of reactions..." << std::endl;	
					{
						// Direct side
						{
							// lnA
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.lnA") );  
			 					Load(lnA__, stream, OPENSMOKE_FORMATTED_FILE);
							}

							// Beta
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.Beta") );  
			 					Load(Beta__, stream, OPENSMOKE_FORMATTED_FILE);

							}

							// E_over_R
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.E_over_R") );  
			 					Load(E_over_R__, stream, OPENSMOKE_FORMATTED_FILE);
							}

							// Global kinetic order of forward reaction
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Direct.ForwardKineticOrder") );  
			 					Load(forward_kinetic_order__, stream, OPENSMOKE_FORMATTED_FILE);
							}
						}

						// Reverse side
						if (number_of_explicitly_reversible_reactions_ != 0)
						{
							// lnA
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Reverse.lnA") );  
			 					Load(lnA_reversible__, stream, OPENSMOKE_FORMATTED_FILE);
							}

							// Beta
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Reverse.Beta") );  
			 					Load(Beta_reversible__, stream, OPENSMOKE_FORMATTED_FILE);

							}

							// E_over_R
							{
								std::stringstream stream;
								stream.str( subtree.get< std::string >("KineticParameters.Reverse.E_over_R") );  
			 					Load(E_over_R_reversible__, stream, OPENSMOKE_FORMATTED_FILE);
							}
						}

						// Stick reactions
						std::cout << " * Reading kinetic parameters of stick reactions..." << std::endl;
						if (number_of_stick_reactions_ != 0)
						{
							stick_power_.resize(number_of_stick_reactions_);
							stick_motz_wise_.resize(number_of_stick_reactions_);
							stick_constant_coefficient_.resize(number_of_stick_reactions_);

							std::stringstream stream;
							stream.str( subtree.get< std::string >("KineticParameters.StickParameters") );  

							for(unsigned int i=0;i<number_of_stick_reactions_;i++)
							{
								std::string dummy;
								stream >> dummy;
								stick_motz_wise_[i] = boost::lexical_cast<bool>(dummy);
								stream >> dummy;
								stick_power_[i] = boost::lexical_cast<double>(dummy);
								stream >> dummy;
								const unsigned int stick_gas_species = boost::lexical_cast<unsigned int>(dummy);

								stream >> dummy;
								const unsigned int n = boost::lexical_cast<unsigned int>(dummy);
							
								double sigma_power = 1.;
								for (unsigned int ii=0;ii<n;ii++)
								{
									stream >> dummy;
									const unsigned int index = boost::lexical_cast<unsigned int>(dummy);
									stream >> dummy;
									const double exponent = boost::lexical_cast<double>(dummy);

								
									const std::string name_site_species = thermodynamics_.vector_names_site_species()[index-1];
									for(unsigned int j=0;j<thermodynamics_.matrix_occupancies_site_species()[0].size();j++)
										for(unsigned int jj=0;jj<thermodynamics_.matrix_occupancies_site_species()[0][j].size();jj++)
											if (thermodynamics_.matrix_names_site_species()[0][j][jj] == name_site_species)
											{	
												const double sigma = thermodynamics_.matrix_occupancies_site_species()[0][j][jj];
												sigma_power *= std::pow(sigma, exponent);
												break;
											}
								}

								stick_constant_coefficient_[i] = sigma_power * sqrt(PhysicalConstants::R_J_kmol/2./PhysicalConstants::pi/thermodynamics_.MW(stick_gas_species-1));
							}
						}

						std::cout << " * Reading kinetic parameters of coverage dependent reactions..." << std::endl;
						if (number_of_coverage_dependent_reactions_ != 0)
						{
							coverage_dependent_species_site_type_.resize(number_of_coverage_dependent_reactions_);
							coverage_dependent_species_index_.resize(number_of_coverage_dependent_reactions_);
							coverage_dependent_eta_.resize(number_of_coverage_dependent_reactions_);
							coverage_dependent_mu_.resize(number_of_coverage_dependent_reactions_);
							coverage_dependent_epsilon_.resize(number_of_coverage_dependent_reactions_);

							std::stringstream stream;
							stream.str( subtree.get< std::string >("KineticParameters.CoverageDependentParameters") );  

							for(unsigned int i=0;i<number_of_coverage_dependent_reactions_;i++)
							{
								std::string dummy;
								stream >> dummy;
								const unsigned int n = boost::lexical_cast<unsigned int>(dummy);

								coverage_dependent_species_site_type_[i].resize(n);
								coverage_dependent_species_index_[i].resize(n);
								coverage_dependent_eta_[i].resize(n);
								coverage_dependent_mu_[i].resize(n);
								coverage_dependent_epsilon_[i].resize(n);

								for(unsigned int j=0;j<n;j++)
								{
									stream >> dummy;
									coverage_dependent_species_site_type_[i][j] = boost::lexical_cast<bool>(dummy);
									stream >> dummy;
									coverage_dependent_species_index_[i][j] = boost::lexical_cast<unsigned int>(dummy);
									stream >> dummy;
									coverage_dependent_eta_[i][j] = boost::lexical_cast<double>(dummy);
									stream >> dummy;
									coverage_dependent_mu_[i][j] = boost::lexical_cast<double>(dummy);
									stream >> dummy;
									coverage_dependent_epsilon_[i][j] = boost::lexical_cast<double>(dummy);
								}
							}
						}

						std::cout << " * Reading kinetic parameters of Langmuir-Hinshelwood reactions..." << std::endl;
						if (number_of_langmuir_reactions_ != 0)
						{
							langmuir_denominator_order_.resize(number_of_langmuir_reactions_);
							langmuir_units_.resize(number_of_langmuir_reactions_);
							langmuir_species_index_.resize(number_of_langmuir_reactions_);
							langmuir_numerator_species_.resize(number_of_langmuir_reactions_);
							langmuir_lnA_.resize(number_of_langmuir_reactions_);
							langmuir_Beta_.resize(number_of_langmuir_reactions_);
							langmuir_H_over_R_.resize(number_of_langmuir_reactions_);
							langmuir_order_.resize(number_of_langmuir_reactions_);

							std::stringstream stream;
							stream.str( subtree.get< std::string >("KineticParameters.LangmuirParameters") );  

							for(unsigned int i=0;i<number_of_langmuir_reactions_;i++)
							{
								std::string dummy;
								stream >> dummy;
								const unsigned int n = boost::lexical_cast<unsigned int>(dummy);

								stream >> dummy;
								langmuir_denominator_order_[i] = boost::lexical_cast<double>(dummy);

								stream >> dummy;
								langmuir_units_[i] = 
									static_cast<PhysicalConstants::UNITS_REACTION_COMPOSITION>(boost::lexical_cast<unsigned int>(dummy));

								langmuir_species_index_[i].resize(n);
								langmuir_numerator_species_[i].resize(n);
								langmuir_lnA_[i].resize(n);
								langmuir_Beta_[i].resize(n);
								langmuir_H_over_R_[i].resize(n);
								langmuir_order_[i].resize(n);

								for(unsigned int j=0;j<n;j++)
								{
									stream >> dummy;
									langmuir_species_index_[i][j] = boost::lexical_cast<unsigned int>(dummy);
								
									stream >> dummy;
									langmuir_numerator_species_[i][j] = boost::lexical_cast<bool>(dummy);
								
									stream >> dummy;
									langmuir_lnA_[i][j] = boost::lexical_cast<double>(dummy);
									if (langmuir_lnA_[i][j] > 0.)
										langmuir_lnA_[i][j] = log(langmuir_lnA_[i][j]);

									stream >> dummy;
									langmuir_Beta_[i][j] = boost::lexical_cast<double>(dummy);
								
									stream >> dummy;
									langmuir_H_over_R_[i][j] = boost::lexical_cast<double>(dummy);
								
									stream >> dummy;
									langmuir_order_[i][j] = boost::lexical_cast<double>(dummy);
								}
							}
						}

						if (number_of_lumped_reactions_ != 0)
						{
							names_of_lumped_functions_.resize(number_of_lumped_reactions_);

							std::stringstream stream;
							stream.str( subtree.get< std::string >("KineticParameters.LumpedParameters") );  

							for (unsigned int i = 0; i<number_of_lumped_reactions_; i++)
								stream >> names_of_lumped_functions_[i];
						}
					}

					// Reactions needing conversion
					{			
						std::stringstream stream;
						stream.str( subtree.get< std::string >("ReactionsNeedingConversion") );  

						std::string dummy;
						stream >> dummy;
						unsigned int number_of_reactions_needing_conversion_ = boost::lexical_cast<unsigned int>(dummy);
						indices_of_reactions_needing_conversion_.resize(number_of_reactions_needing_conversion_);

						for (unsigned int j=0;j<number_of_reactions_needing_conversion_;j++)
						{		
							stream >> dummy;
							indices_of_reactions_needing_conversion_[j] = boost::lexical_cast<unsigned int>(dummy);
						}
					}

					// Reactions with non conservation of sites
					{			
						std::stringstream stream;
						stream.str( subtree.get< std::string >("ConservationOfSites") );  

						std::string dummy;
						stream >> dummy;
						unsigned int number_of_reactions_with_non_conservation_of_sites = boost::lexical_cast<unsigned int>(dummy);
						non_conservation_of_sites_indices_of_reactions_.resize(number_of_reactions_with_non_conservation_of_sites);
						non_conservation_of_sites_phase_of_reactions_.resize(number_of_reactions_with_non_conservation_of_sites);
						non_conservation_of_sites_delta_sigma_.resize(number_of_reactions_with_non_conservation_of_sites);

						for (unsigned int j=0;j<number_of_reactions_with_non_conservation_of_sites;j++)
						{		
							stream >> dummy;
							non_conservation_of_sites_indices_of_reactions_[j] = boost::lexical_cast<unsigned int>(dummy);
							stream >> dummy;
							non_conservation_of_sites_phase_of_reactions_[j] = boost::lexical_cast<unsigned int>(dummy);
							stream >> dummy;
							non_conservation_of_sites_delta_sigma_[j] = boost::lexical_cast<double>(dummy);
						}
					}

					// Thermodynamic kinetic constants
					{
						std::stringstream stream;
						stream.str( subtree.get< std::string >("ThermodynamicReversibleReactions") );  

						std::string dummy;
						stream >> dummy;
						unsigned int n = boost::lexical_cast<unsigned int>(dummy);

						delta_sigma_times_log_Gamma0_.resize(n);
						delta_nu_gas_.resize(n);

						for (unsigned int j=0;j<n;j++)
						{
							stream >> dummy;
							const unsigned int index_phase = boost::lexical_cast<unsigned int>(dummy);
							stream >> dummy;
							const double delta_sigma = boost::lexical_cast<double>(dummy);
							stream >> dummy;
							delta_nu_gas_[j] = boost::lexical_cast<double>(dummy);

							delta_sigma_times_log_Gamma0_[j] = delta_sigma * log(thermodynamics_.matrix_densities_site_phases()[0][index_phase-1]);
						}				
					}

					// Stoichiometry
					{
						std::string stoichiometry_type =  subtree.get<std::string>("Stoichiometry.<xmlattr>.type");
						std::string stoichiometry_version = subtree.get<std::string>("Stoichiometry.<xmlattr>.version");
						if (stoichiometry_type != "OpenSMOKE" || stoichiometry_version != "01-02-2014")
							ErrorMessage("void KineticsMap_Surface_CHEMKIN::ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree)", "The current stoichiometric data are not supported.");

						std::stringstream stream;
						stream.str( subtree.get< std::string >("Stoichiometry") ); 

						stoichiometry_ = new StoichiometricMap(this->number_of_species_, this->number_of_reactions_);
						stoichiometry_->ReadFromASCIIFile(stream);
						changeOfMoles__ = stoichiometry_->ChangeOfMoles();

						{
							OpenSMOKE::OpenSMOKEVectorBool tmp(this->number_of_reactions_);
							tmp = false;
							for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
								tmp[indices_of_thermodynamic_reversible_reactions__[k]] = true;
							stoichiometry_->CompleteChangeOfMoles(tmp.GetHandle());
	 					}
			
		

						// Memory allocation
						aux_vector__.resize(this->number_of_species_);
						std::fill(aux_vector__.begin(), aux_vector__.end(), 0.);

						cSites__.resize(thermodynamics_.number_of_site_species());
						std::fill(cSites__.begin(), cSites__.end(), 0.);

						c__.resize(this->number_of_species_);
						std::fill(c__.begin(), c__.end(), 0.);

						reaction_s_over_R__.resize(this->number_of_reactions_);
						std::fill(reaction_s_over_R__.begin(), reaction_s_over_R__.end(), 0.);

						reaction_h_over_RT__.resize(this->number_of_reactions_);
						std::fill(reaction_h_over_RT__.begin(), reaction_h_over_RT__.end(), 0.);

						kArrhenius__.resize(this->number_of_reactions_);
						std::fill(kArrhenius__.begin(), kArrhenius__.end(), 0.);

						kArrheniusModified__.resize(this->number_of_reactions_);
						std::fill(kArrheniusModified__.begin(), kArrheniusModified__.end(), 0.);

						uKeq__.resize(number_of_thermodynamic_reversible_reactions_);
						std::fill(uKeq__.begin(), uKeq__.end(), 0.);

						kArrhenius_reversible__.resize(number_of_explicitly_reversible_reactions_);
						std::fill(kArrhenius_reversible__.begin(), kArrhenius_reversible__.end(), 0.);

						forwardReactionRates__.resize(this->number_of_reactions_);
						std::fill(forwardReactionRates__.begin(), forwardReactionRates__.end(), 0.);
					
						reverseReactionRates__.resize(this->number_of_reactions_);
						std::fill(reverseReactionRates__.begin(), reverseReactionRates__.end(), 0.);
					
						netReactionRates__.resize(this->number_of_reactions_);
						std::fill(netReactionRates__.begin(), netReactionRates__.end(), 0.);

						isThermodynamicallyReversible__.resize(this->number_of_reactions_);
						std::fill(isThermodynamicallyReversible__.begin(), isThermodynamicallyReversible__.end(), 0);
					
						isExplicitlyReversible__.resize(this->number_of_reactions_);
						std::fill(isExplicitlyReversible__.begin(), isExplicitlyReversible__.end(), 0);

						for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
							isThermodynamicallyReversible__[indices_of_thermodynamic_reversible_reactions__[k]-1] = k+1;
						for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
							isExplicitlyReversible__[indices_of_explicitly_reversible_reactions__[k]-1] = k+1;

			
						// Additional indices for sensitivity analysis
						{						
							local_family_index__.resize(this->number_of_reactions_);
							std::fill(local_family_index__.begin(), local_family_index__.end(), 0);

							type_of_reaction__.resize(this->number_of_reactions_);
							std::fill(type_of_reaction__.begin(), type_of_reaction__.end(), PhysicalConstants::REACTION_SIMPLE);
						}

						std::cout << std::endl;
						std::cout << "----------------------------------------------------------------------------" << std::endl;
						std::cout << " Kinetic Mechanism Summary"<< std::endl;
						std::cout << "----------------------------------------------------------------------------" << std::endl;
						std::cout << " Total number of species:          " << this->number_of_species_ << std::endl;
						std::cout << " Total number of reactions:        " << this->number_of_reactions_ << std::endl;
						std::cout << "   Reversible reactions:           " << number_of_reversible_reactions_ << " (" << number_of_reversible_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
						std::cout << "    * by thermodynamics:           " << number_of_thermodynamic_reversible_reactions_ << " (" << number_of_thermodynamic_reversible_reactions_/std::max(1.,double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
						std::cout << "    * by Arrhenius' law:           " << number_of_explicitly_reversible_reactions_ << " (" << number_of_explicitly_reversible_reactions_/std::max(1.,double(number_of_reversible_reactions_))*100. << "%)" << std::endl;
						std::cout << "   Stick reactions:                " << number_of_stick_reactions_ << " (" << number_of_stick_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
						std::cout << "   Coverage dependent reactions:   " << number_of_coverage_dependent_reactions_ << " (" << number_of_coverage_dependent_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
						std::cout << "   Langmuir-Hinshelwood reactions: " << number_of_langmuir_reactions_ << " (" << number_of_langmuir_reactions_/std::max(1.,double(this->number_of_reactions_))*100. << "%)" << std::endl;
						std::cout << "   Lumped reactions:               " << number_of_lumped_reactions_ << " (" << number_of_lumped_reactions_ / std::max(1., double(this->number_of_reactions_))*100. << "%)" << std::endl;
						std::cout << std::endl;

						stoichiometry_->Summary(std::cout);
					}
				}
			}
		}
	}

	void KineticsMap_Surface_CHEMKIN::ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree)
	{
		try
		{
			this->number_of_species_ = ptree.get<unsigned int>("opensmoke.NumberOfSpecies"); 				
		}
		catch(...)
		{
			ErrorMessage("KineticsMap_Surface_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the number of species.");
		}
	}

	void KineticsMap_Surface_CHEMKIN::ReactionEnthalpiesAndEntropies()
	{
		if (reaction_h_and_s_must_be_recalculated_ == true)
		{
			stoichiometry_->ReactionEnthalpyAndEntropy(	reaction_h_over_RT__, reaction_s_over_R__,
														thermodynamics_.Species_H_over_RT(), thermodynamics_.Species_S_over_R() );

			reaction_h_and_s_must_be_recalculated_ = false;
		}
	}

	void KineticsMap_Surface_CHEMKIN::KineticConstants()
	{
		ReactionEnthalpiesAndEntropies();

		if (arrhenius_kinetic_constants_must_be_recalculated_ == true)
		{
			// Forward kinetic constants (Arrhenius' Law)
			{
				double *pt_lnA = lnA__.data();
				double *pt_Beta = Beta__.data();
				double *pt_E_over_R = E_over_R__.data();
				double *pt_kArrheniusT = kArrhenius__.data();
			
				for(unsigned int j=0;j<this->number_of_reactions_;j++)
					*pt_kArrheniusT++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius__, &kArrhenius__);
			}

			// Equilibrium constants (inverse value)
			{
				for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
				{
					unsigned int j = indices_of_thermodynamic_reversible_reactions__[k];
					uKeq__[k] = -reaction_s_over_R__[j-1] + reaction_h_over_RT__[j-1] - log_Patm_over_RT_ * delta_nu_gas_[k] - delta_sigma_times_log_Gamma0_[k];
				}
				Exp(uKeq__, &uKeq__);
			}

			// Explicit reverse Arrhenius constants
			if (number_of_explicitly_reversible_reactions_ != 0)
			{
				double *pt_lnA = lnA_reversible__.data();
				double *pt_Beta = Beta_reversible__.data();
				double *pt_E_over_R = E_over_R_reversible__.data();
				double *pt_kArrhenius = kArrhenius_reversible__.data();

				for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
					*pt_kArrhenius++ = (*pt_lnA++) + (*pt_Beta++)*this->logT_ - (*pt_E_over_R++)*this->uT_;
			
				Exp(kArrhenius_reversible__, &kArrhenius_reversible__);
			}

			arrhenius_kinetic_constants_must_be_recalculated_ = false;
		}

		//if (nonconventional_kinetic_constants_must_be_recalculated_ == true)
		{
			nonconventional_kinetic_constants_must_be_recalculated_ = false;
		}

		// Conversions
		if (indices_of_reactions_needing_conversion_.size() > 0)
		{
			const double R_times_T = PhysicalConstants::R_J_kmol*this->T_;
			for(unsigned int k=0;k<indices_of_reactions_needing_conversion_.size();k++)
			{
				unsigned int j = indices_of_reactions_needing_conversion_[k];
				kArrhenius__[j-1] *= std::pow(R_times_T, forward_kinetic_order__[j-1]);
			}
		}

		kArrheniusModified__ = kArrhenius__;
	}

	void KineticsMap_Surface_CHEMKIN::ReactionRates(const double* cGas, const double* z, const double* a, const double* Gamma)
	{
		double cTot = 0.;
		for (unsigned int j = 0; j<thermodynamics_.number_of_gas_species(); j++)
			cTot += cGas[j];

		for (unsigned int j = 0; j<thermodynamics_.number_of_site_species(); j++)
			cSites__[j] = z[j]*Gamma[thermodynamics_.vector_site_phases_belonging()[j]] / 
							thermodynamics_.vector_occupancies_site_species()[j];
		
		unsigned int count = 0;
		for(unsigned int j=0;j<thermodynamics_.number_of_gas_species();j++)
			c__[count++] = cGas[j];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_site_species();j++)
			c__[count++] = cSites__[j];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_bulk_species();j++)
			c__[count++] = a[j];
		
		if (type_of_kinetics_ == TYPE_OF_KINETICS_CHEMKIN_CONVENTIONAL)
		{
			// 1. Kinetic constants
			KineticConstants();
		
			// 2. Correct the effective kinetic constants by stick reactions
			if (number_of_stick_reactions_>0)
			{
				double total_site_density = 0.;
				for (unsigned int i = 0;i<thermodynamics_.number_of_site_phases(0);i++)
					total_site_density = Gamma[i];

				for(unsigned int s=1;s<=number_of_stick_reactions_;s++)
				{
					const unsigned int j=indices_of_stick_reactions__[s-1];

					const double gamma = std::min(1., kArrhenius__[j-1]);
					if (stick_motz_wise_[s-1] == false)	kArrheniusModified__[j-1] = gamma;
					else                                kArrheniusModified__[j-1] = gamma/(1.-gamma/2.);

					kArrheniusModified__[j-1] *= stick_constant_coefficient_[s-1]*std::sqrt(this->T_)/std::pow(total_site_density, stick_power_[s-1]);
				}
			}

			// 3. Correct the effective kinetic constants by coverage dependent reactions
			if (number_of_coverage_dependent_reactions_>0)
			{
				for(unsigned int s=1;s<=number_of_coverage_dependent_reactions_;s++)
				{
					double correction = 0.;
					for(unsigned int k=0;k<coverage_dependent_species_site_type_[s-1].size();k++)
					{
						double value;
						if (coverage_dependent_species_site_type_[s-1][k] == true)
							value = z[coverage_dependent_species_index_[s-1][k]-1];
						else
							value = a[coverage_dependent_species_index_[s-1][k]-1];

						const double eps_value = 1.e-20;
						correction +=	PhysicalConstants::ln_10*coverage_dependent_eta_[s-1][k]*value + 
										coverage_dependent_mu_[s-1][k]*log(value+eps_value) -
										coverage_dependent_epsilon_[s-1][k]*value/(PhysicalConstants::R_J_kmol*this->T_);
						
					}
					correction = std::exp(correction);
				
					const unsigned int j=indices_of_coverage_dependent_reactions__[s-1];
					kArrheniusModified__[j-1] *= correction;
				}
			}

			// 4. Correct the effective kinetic constants by Langmuir-Hinshelwood reactions
			if (number_of_langmuir_reactions_>0)
			{
				for(unsigned int s=1;s<=number_of_langmuir_reactions_;s++)
				{
					double correction_denominator = 1.;
					double correction_numerator = 1.;
					for(unsigned int k=0;k<langmuir_species_index_[s-1].size();k++)
					{
						const double  K = std::exp( langmuir_lnA_[s-1][k] + langmuir_Beta_[s-1][k]*log(this->T_) - langmuir_H_over_R_[s-1][k]/this->T_ );

						double value = cGas[langmuir_species_index_[s-1][k]-1];
					
						if (langmuir_units_[s-1] != PhysicalConstants::UNITS_STD)
							value *= PhysicalConstants::R_J_kmol*this->T_;

						correction_denominator += K*std::pow(value, langmuir_order_[s-1][k]);

						if (langmuir_numerator_species_[s-1][k] == true)
							correction_numerator *= K;
					}

					const unsigned int j=indices_of_langmuir_reactions__[s-1];
					kArrheniusModified__[j-1] *= correction_numerator/std::pow(correction_denominator,langmuir_denominator_order_[s-1]);
				}
			}
		}
		
		else if (type_of_kinetics_ == TYPE_OF_KINETICS_UBI_QEP)
		{
			ubiqep_submechanism_->CalculateDissociationEnergies(thermodynamics_.Species_H_over_RT());
			ubiqep_submechanism_->CalculateChemisorptionHeats(this->T_, z);
			ubiqep_submechanism_->CalculateSurfaceEnthalpies(this->T_);
			ubiqep_submechanism_->CalculateActivationEnergies();

			double total_site_density = 0.;
			for (unsigned int i = 0; i<thermodynamics_.number_of_site_phases(0); i++)
				total_site_density = Gamma[i];
			ubiqep_submechanism_->ForwardKineticConstants(this->T_, total_site_density, kArrheniusModified__.data());
		}

		// Calculates the product of conenctrations (for forward and reverse reactions)
		// Be careful: the reverseReactionRates_ vector is defined for all the reactions
		// in the kinetic scheme, not only for the reversible reactions. After calling the
		// function reported below the value of reverseReactionRates_ vector for non reversible
		// reactions is put equal to 1.
		stoichiometry_->ProductOfConcentrations(forwardReactionRates__, reverseReactionRates__, c__.data());

		// Corrects the product of concentrations for reverse reaction by the 
		// thermodynamic equilibrium constant
		for(unsigned int k=0;k<number_of_thermodynamic_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_thermodynamic_reversible_reactions__[k]-1;
			reverseReactionRates__[j] *= uKeq__[k];
		}

		// Corrects the product of concentrations for reverse reaction by the 
		// explicit Arrhenius kinetic parameters (if provided)
		for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k]-1;
			reverseReactionRates__[j] *= kArrhenius_reversible__[k]/kArrhenius__[j];
		}

		// Calculates the net reaction rate
		// Be careful: the netReactionRates_ vector must be multiplied by the effective 
		// forward kinetic constant, to obtain the real reaction rates in [kmol/m3/s]
		netReactionRates__ = forwardReactionRates__;
		for(unsigned int k=0;k<number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k]-1;
			netReactionRates__[j] -= reverseReactionRates__[j];
		}

		// Multiplies the net reaction rate by the effective kinetic constant (accounting for 
		// third-body effects, fall-off, etc.). At the end of this function the netReactionRates_
		// vector contains the net reaction rates of all the reactions in [kmol/m3/s]
		ElementByElementProduct(netReactionRates__.size(), netReactionRates__.data(), kArrheniusModified__.data(), netReactionRates__.data());

		// User defined reaction rates
		if (number_of_lumped_reactions_ != 0)
			UserDefinedReactionRates(cGas, z, a, Gamma);
	}

	void KineticsMap_Surface_CHEMKIN::UserDefinedReactionRates(const double* cGas, const double* z, const double* a, const double* Gamma)
	{
		FatalErrorMessage("KineticsMap_Surface_CHEMKIN::UserDefinedReactionRates: No user defined reaction rates are provided by the user!");
	}

	void KineticsMap_Surface_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa) const
	{
		stoichiometry_->RateOfProductionAnalysis(netReactionRates__.data(), false);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	void KineticsMap_Surface_CHEMKIN::RateOfProductionAnalysis(ROPA_Data& ropa, const double* rf, const double* rb) const
	{
		stoichiometry_->RateOfProductionAnalysis(rf, rb);
		stoichiometry_->WriteRateOfProductionAnalysis(ropa);
	}

	void KineticsMap_Surface_CHEMKIN::ProductionAndDestructionRates(double* P, double* D)
	{
		stoichiometry_->ProductionAndDestructionRatesFromReactionRates(P, D, netReactionRates__.data());
	}

	void KineticsMap_Surface_CHEMKIN::FormationRates(double* R)
	{
		stoichiometry_->FormationRatesFromReactionRates(R, netReactionRates__.data());
	}

	void KineticsMap_Surface_CHEMKIN::FormationRates(double* Rgas, double* Rsite, double* Rbulk, double* RsitePhases)
	{
		std::vector<double> R(thermodynamics_.NumberOfSpecies());
		stoichiometry_->FormationRatesFromReactionRates(R.data(), netReactionRates__.data());
		
		unsigned int count = 0;
		for(unsigned int j=0;j<thermodynamics_.number_of_gas_species();j++)
			Rgas[j] = R[count++];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_site_species();j++)
			Rsite[j] = R[count++];
		
		for(unsigned int j=0;j<thermodynamics_.number_of_bulk_species();j++)
			Rbulk[j] = R[count++];
		
		for (unsigned int i = 0; i<thermodynamics_.number_of_site_phases(0); i++)
			RsitePhases[i] = 0.;
		for(unsigned int j=0;j<non_conservation_of_sites_indices_of_reactions_.size();j++)
			RsitePhases[non_conservation_of_sites_phase_of_reactions_[j]-1] += 
				netReactionRates__[non_conservation_of_sites_indices_of_reactions_[j]-1] * non_conservation_of_sites_delta_sigma_[j];
	}

	double KineticsMap_Surface_CHEMKIN::HeatRelease(const double* Rgas, const double* Rsurface, const double* Rbulk)
	{
		unsigned int k = 0;
		for (unsigned int j = 0; j<thermodynamics_.number_of_gas_species(); j++)
			aux_vector__[k++] = Rgas[j];
		for (unsigned int j = 0; j<thermodynamics_.number_of_site_species(); j++)
			aux_vector__[k++] = Rsurface[j];
		for (unsigned int j = 0; j<thermodynamics_.number_of_bulk_species(); j++)
			aux_vector__[k++] = Rbulk[j];

		return -Dot(aux_vector__.size(), aux_vector__.data(), thermodynamics_.Species_H_over_RT().data()) * PhysicalConstants::R_J_kmol * this->T_;
	}

	const std::vector<double>& KineticsMap_Surface_CHEMKIN::GiveMeReactionRates()
	{
		return netReactionRates__;
	}

	void KineticsMap_Surface_CHEMKIN::GiveMeReactionRates(double* r)
	{
		r = netReactionRates__.data();
	}

	void KineticsMap_Surface_CHEMKIN::GetForwardReactionRates(double* r)
	{
		ElementByElementProduct(forwardReactionRates__.size(), forwardReactionRates__.data(), kArrheniusModified__.data(), r);
	}

	void KineticsMap_Surface_CHEMKIN::GetBackwardReactionRates(double* r)
	{
		for (unsigned int j = 0; j < this->number_of_reactions_; j++)
			r[j] = 0.;

		for(unsigned int k=0;k<number_of_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_reversible_reactions__[k]-1;
			r[j] = reverseReactionRates__[j]*kArrheniusModified__[j];
		}
		for(unsigned int k=0;k<number_of_explicitly_reversible_reactions_;k++)
		{
			unsigned int j = indices_of_explicitly_reversible_reactions__[k]-1;
			r[j] = reverseReactionRates__[j]*kArrheniusModified__[j];
		}                
	}

	void KineticsMap_Surface_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k)
	{			
		thermodynamics_.SetPressure(101325.);
		thermodynamics_.SetTemperature(298.15);
		SetTemperature(298.15);
		SetPressure(101325.);
		ReactionEnthalpiesAndEntropies();
		
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << (reaction_h_over_RT__[k-1]-reaction_s_over_R__[k-1])*PhysicalConstants::R_kcal_mol*this->T_;	// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_h_over_RT__[k-1] *PhysicalConstants::R_kcal_mol*this->T_;;						// [kcal/mol]
		fOut    << std::setw(14) << std::left << std::setprecision(3) << std::fixed << reaction_s_over_R__[k-1] *PhysicalConstants::R_cal_mol;;							// [cal/mol/K]
	}

	void KineticsMap_Surface_CHEMKIN::WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const double conversion_forward, const double conversion_backward)
	{			
		// TODO
	}

	void KineticsMap_Surface_CHEMKIN::FittedReverseKineticConstants(const double* x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters)
	{			
		// TODO
	}

	void KineticsMap_Surface_CHEMKIN::FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters)
	{
		if (isThermodynamicallyReversible__[k-1] != 0)
		{
			const unsigned int j = isThermodynamicallyReversible__[k-1];

			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(fittedKineticParameters(0, j-1));
			if (fittedKineticParameters.rows() == 2)
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << 0.;
			else
				fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << fittedKineticParameters(2, j - 1);
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << fittedKineticParameters(1,j-1)/Conversions::J_from_kcal;
			fOut << std::setw(5)  << "";
		}
		else if (isExplicitlyReversible__[k-1] != 0)
		{
			const unsigned int j = isExplicitlyReversible__[k-1];
				
			fOut << std::setw(18) << std::right << std::scientific << std::setprecision(4) << std::exp(lnA_reversible__[j-1]);
			fOut << std::setw(10) << std::right << std::fixed << std::setprecision(3) << Beta_reversible__[j-1];
			fOut << std::setw(16) << std::right << std::fixed << std::setprecision(2) << E_over_R_reversible__[j-1];
			fOut << std::setw(5)  << "";
		}
	}

	void UBIQEP_SubMechanism::ReadFromXMLFile(boost::property_tree::ptree& ptree, const std::vector<double>& MWs)
	{

		std::cout << "Reading Heats of chemisorption..." << std::endl;
		{
			std::stringstream stream;
			stream.str( ptree.get< std::string >("HeatsOfChemisorption") );  

			std::string dummy;

			stream >> dummy;
			number_of_gas_species_ = boost::lexical_cast<unsigned int>(dummy);

			stream >> dummy;
			number_of_site_species_ = boost::lexical_cast<unsigned int>(dummy);

			stream >> dummy;
			unsigned int number_of_gas_species_chemisorption_heats = boost::lexical_cast<unsigned int>(dummy);

			stream >> dummy;
			reference_temperature_ = boost::lexical_cast<double>(dummy);

			unsigned int size = number_of_site_species_+number_of_gas_species_chemisorption_heats;
			
			chemisorption_heats_constant_coefficient__.resize(size);
			std::fill(chemisorption_heats_constant_coefficient__.begin(), chemisorption_heats_constant_coefficient__.end(), 0.);
			
			chemisorption_heats_temperature_coefficient__.resize(size);
			std::fill(chemisorption_heats_temperature_coefficient__.begin(), chemisorption_heats_temperature_coefficient__.end(), 0.);

			QStar__.resize(size);
			std::fill(QStar__.begin(), QStar__.end(), 0.);

			chemisorption_heats_gas_indices__.resize(number_of_gas_species_chemisorption_heats);
			std::fill(chemisorption_heats_gas_indices__.begin(), chemisorption_heats_gas_indices__.end(), 0);

			chemisorption_heats_coefficients__ = new std::vector<double>[size];
			chemisorption_heats_indices__ = new std::vector<unsigned int>[size];
			for (unsigned int j=1;j<=number_of_site_species_+number_of_gas_species_chemisorption_heats;j++)
			{		
				if (j > number_of_site_species_)
				{
					stream >> dummy;
					chemisorption_heats_gas_indices__[j-number_of_site_species_-1] = boost::lexical_cast<unsigned int>(dummy);
				}

				stream >> dummy;
				chemisorption_heats_temperature_coefficient__[j-1] = boost::lexical_cast<double>(dummy);

				stream >> dummy;
				chemisorption_heats_constant_coefficient__[j-1] = boost::lexical_cast<double>(dummy);

				stream >> dummy;
				unsigned int n = boost::lexical_cast<unsigned int>(dummy);

				chemisorption_heats_coefficients__[j-1].resize(n);
				std::fill(chemisorption_heats_coefficients__[j-1].begin(), chemisorption_heats_coefficients__[j-1].end(), 0.);

				chemisorption_heats_indices__[j-1].resize(n);
				std::fill(chemisorption_heats_indices__[j-1].begin(), chemisorption_heats_indices__[j-1].end(), 0);

				for (unsigned int k=0;k<n;k++)
				{
					stream >> dummy;
					chemisorption_heats_indices__[j-1][k] = boost::lexical_cast<unsigned int>(dummy);

					stream >> dummy;
					chemisorption_heats_coefficients__[j-1][k] = boost::lexical_cast<double>(dummy);
				}
			}
		}

		std::cout << "Reading reaction parameters (direct reactions)" << std::endl;
		{
			std::stringstream stream;
			stream.str( ptree.get< std::string >("UBIParameters.UBIDirect") );  

			std::string dummy;

			stream >> dummy;
			half_number_of_ubiqep_reactions_ = boost::lexical_cast<unsigned int>(dummy);
			
			dissociation_energies__.resize(half_number_of_ubiqep_reactions_);
			std::fill(dissociation_energies__.begin(), dissociation_energies__.end(), 0.);

			surface_enthalpies__.resize(half_number_of_ubiqep_reactions_);
			std::fill(surface_enthalpies__.begin(), surface_enthalpies__.end(), 0.);

			E_forward__.resize(half_number_of_ubiqep_reactions_);
			std::fill(E_forward__.begin(), E_forward__.end(), 0.);

			E_backward__.resize(half_number_of_ubiqep_reactions_);
			std::fill(E_backward__.begin(), E_backward__.end(), 0.);

			adsorption_coefficient__.resize(half_number_of_ubiqep_reactions_);
			std::fill(adsorption_coefficient__.begin(), adsorption_coefficient__.end(), 0.);

			sigma_direct__.resize(half_number_of_ubiqep_reactions_);
			std::fill(sigma_direct__.begin(), sigma_direct__.end(), 0.);

			Beta_direct__.resize(half_number_of_ubiqep_reactions_);
			std::fill(Beta_direct__.begin(), Beta_direct__.end(), 0.);

			lnA_direct__.resize(half_number_of_ubiqep_reactions_);
			std::fill(lnA_direct__.begin(), lnA_direct__.end(), 0.);

			lambda__.resize(half_number_of_ubiqep_reactions_);
			std::fill(lambda__.begin(), lambda__.end(), 0.);
	
			ubiqep_class__.resize(half_number_of_ubiqep_reactions_);
			std::fill(ubiqep_class__.begin(), ubiqep_class__.end(), 0);

			ubiqep_type_direct__.resize(half_number_of_ubiqep_reactions_);

			index_A__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_A__.begin(), index_A__.end(), 0);
			
			index_B__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_B__.begin(), index_B__.end(), 0);

			index_C__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_C__.begin(), index_C__.end(), 0);

			index_D__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_D__.begin(), index_D__.end(), 0);

			index_Star__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_Star__.begin(), index_Star__.end(), 0);

			index_A2__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_A2__.begin(), index_A2__.end(), 0);

			index_AB__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_AB__.begin(), index_AB__.end(), 0);

			index_A2_Chemisorption_Heats__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_A2_Chemisorption_Heats__.begin(), index_A2_Chemisorption_Heats__.end(), 0);
			
			index_AB_Chemisorption_Heats__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_AB_Chemisorption_Heats__.begin(), index_AB_Chemisorption_Heats__.end(), 0);

			index_ABStar__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_ABStar__.begin(), index_ABStar__.end(), 0);

			index_AStar__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_AStar__.begin(), index_AStar__.end(), 0);

			index_BStar__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_BStar__.begin(), index_BStar__.end(), 0);

			index_CStar__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_CStar__.begin(), index_CStar__.end(), 0);

			index_DStar__.resize(half_number_of_ubiqep_reactions_);
			std::fill(index_DStar__.begin(), index_DStar__.end(), 0);
			
			for(unsigned int k=0;k<half_number_of_ubiqep_reactions_;k++)
			{
				stream >> dummy;
				lnA_direct__[k] = boost::lexical_cast<double>(dummy);
				lnA_direct__[k] = lnA_direct__[k] > 0. ? log(lnA_direct__[k]) : lnA_direct__[k];
				stream >> dummy;
				Beta_direct__[k] = boost::lexical_cast<double>(dummy);
				stream >> dummy;
				sigma_direct__[k] = boost::lexical_cast<double>(dummy);

				stream >> dummy;
				ubiqep_class__[k] = boost::lexical_cast<unsigned int>(dummy);
				
				stream >> dummy;
				ubiqep_type_direct__[k] = static_cast<PhysicalConstants::UBIQEP_TYPE>(boost::lexical_cast<unsigned int>(dummy));

				if (ubiqep_type_direct__[k] == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
				{
					stream >> dummy;
					const unsigned int gas_index = boost::lexical_cast<unsigned int>(dummy);
					adsorption_coefficient__[k] = sqrt(PhysicalConstants::R_J_kmol/2/PhysicalConstants::pi/MWs[gas_index-1]);
				}

				if (ubiqep_class__[k] == 0)
				{	
				}
				else if (ubiqep_class__[k] == 1)
				{
					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_Star__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 1.;
				}
				else if (ubiqep_class__[k] == 2 || ubiqep_class__[k] == 3)
				{
					stream >> dummy;
					index_A2__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_Star__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 2.;

					// Gas phase
					for (unsigned int j=1;j<=chemisorption_heats_gas_indices__.size();j++)
						if (index_A2__[k] == chemisorption_heats_gas_indices__[j-1]-1)
						{
							index_A2_Chemisorption_Heats__[k] = number_of_site_species_+j-1;
							break;
						}
				}
				else if (ubiqep_class__[k] == 4)
				{
					stream >> dummy;
					index_AB__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_Star__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_BStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_B__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 2.;

					// Gas phase
					for (unsigned int j=1;j<=chemisorption_heats_gas_indices__.size();j++)
						if (index_AB__[k] == chemisorption_heats_gas_indices__[j-1]-1)
						{
							index_AB_Chemisorption_Heats__[k] = number_of_site_species_+j-1;
							break;
						}
				}
				else if (ubiqep_class__[k] == 5)
				{
					stream >> dummy;
					index_ABStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_Star__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_BStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_B__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AB__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 2.;
				}
				else if (ubiqep_class__[k] == 6)
				{
					stream >> dummy;
					index_CStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_DStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_BStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_B__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_C__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_D__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 2.;
				}
				else if (ubiqep_class__[k] == 7)
				{
					stream >> dummy;
					index_CStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_BStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_B__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_C__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 2.;
				}
				else if (ubiqep_class__[k] == 8)
				{
					stream >> dummy;
					index_CStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_DStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_AStar__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_A__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_C__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					stream >> dummy;
					index_D__[k] = boost::lexical_cast<unsigned int>(dummy)-1;

					lambda__[k] = 2.;
				}
			}
		}

		std::cout << "Reading reaction parameters (reverse reactions)" << std::endl;
		{
			std::stringstream stream;
			stream.str( ptree.get< std::string >("UBIParameters.UBIReverse") );  

			std::string dummy;

			stream >> dummy;
			half_number_of_ubiqep_reactions_ = boost::lexical_cast<unsigned int>(dummy);
			
			sigma_reverse__.resize(half_number_of_ubiqep_reactions_);
			std::fill(sigma_reverse__.begin(), sigma_reverse__.end(), 0.);

			Beta_reverse__.resize(half_number_of_ubiqep_reactions_);
			std::fill(Beta_reverse__.begin(), Beta_reverse__.end(), 0.);

			lnA_reverse__.resize(half_number_of_ubiqep_reactions_);
			std::fill(lnA_reverse__.begin(), lnA_reverse__.end(), 0.);

			ubiqep_type_reverse__.resize(half_number_of_ubiqep_reactions_);

			for(unsigned int k=0;k<half_number_of_ubiqep_reactions_;k++)
			{
				stream >> dummy;
				lnA_reverse__[k] = boost::lexical_cast<double>(dummy);
				lnA_reverse__[k] = lnA_reverse__[k] > 0. ? log(lnA_reverse__[k]) : lnA_reverse__[k];
				stream >> dummy;
				Beta_reverse__[k] = boost::lexical_cast<double>(dummy);
				stream >> dummy;
				sigma_reverse__[k] = boost::lexical_cast<double>(dummy);

				stream >> dummy;	// class of reaction
				stream >> dummy;    // 
				ubiqep_type_reverse__[k] = static_cast<PhysicalConstants::UBIQEP_TYPE>(boost::lexical_cast<unsigned int>(dummy));
			}
		}
	}

	void UBIQEP_SubMechanism::CalculateChemisorptionHeats(const double T, const double* Z)
	{
		const double deltaT_times_R = (T-reference_temperature_)*PhysicalConstants::R_kcal_mol;
		
		for (unsigned int j=0;j<QStar__.size();j++)
		{
			QStar__[j] =	chemisorption_heats_constant_coefficient__[j] - 
							chemisorption_heats_temperature_coefficient__[j] * deltaT_times_R;

			for (unsigned int k=0;k<chemisorption_heats_indices__[j].size();k++)
				QStar__[j] += chemisorption_heats_coefficients__[j][k] * Z[chemisorption_heats_indices__[j][k]-1];
		}
	}

	void UBIQEP_SubMechanism::CalculateDissociationEnergies(const std::vector<double>& H)
	{
		for (unsigned int k=0;k<half_number_of_ubiqep_reactions_;k++)
		{
			switch(ubiqep_class__[k])
			{
			case 1:
				dissociation_energies__[k] = 0;	
				break;
			case 2:
				dissociation_energies__[k] = 2.*H[index_A__[k]] - H[index_A2__[k]];
				break;
			case 3:
				dissociation_energies__[k] = 2.*H[index_A__[k]] - H[index_A2__[k]];
				break;
			case 4:
				dissociation_energies__[k] = H[index_A__[k]] + H[index_B__[k]] - H[index_AB__[k]];
				break;
			case 5:
				dissociation_energies__[k] = H[index_A__[k]] + H[index_B__[k]] - H[index_AB__[k]];
				break;
			case 6:
				dissociation_energies__[k] = H[index_A__[k]] + H[index_B__[k]] - H[index_C__[k]]- H[index_D__[k]];
				break;
			case 7:
				dissociation_energies__[k] = H[index_A__[k]] + H[index_B__[k]] - 2.*H[index_C__[k]];
				break;
			case 8:
				dissociation_energies__[k] = 2.*H[index_A__[k]] - H[index_C__[k]] - H[index_D__[k]];
				break;
			}

			dissociation_energies__[k] = std::fabs(dissociation_energies__[k]);
		}
	}

	void UBIQEP_SubMechanism::CalculateSurfaceEnthalpies(const double T)
	{
		const double R_times_T = PhysicalConstants::R_kcal_mol * T;

		for (unsigned int k=0;k<half_number_of_ubiqep_reactions_;k++)
		{
			switch(ubiqep_class__[k])
			{
			case 1:
				surface_enthalpies__[k] = ( QStar__[index_Star__[k]] - QStar__[index_AStar__[k]] );
				break;
			case 2:
				surface_enthalpies__[k] = 2. * ( QStar__[index_Star__[k]] - QStar__[index_AStar__[k]] );
				break;
			case 3:
				surface_enthalpies__[k] = 2. * ( QStar__[index_Star__[k]] - QStar__[index_AStar__[k]] );
				break;
			case 4:
				surface_enthalpies__[k] = 2.*QStar__[index_Star__[k]] - (QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]]);
				break;
			case 5:
				surface_enthalpies__[k] = QStar__[index_ABStar__[k]] + QStar__[index_Star__[k]] -
											(QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]] );
				break;
			case 6:
				surface_enthalpies__[k] = QStar__[index_CStar__[k]] + QStar__[index_DStar__[k]] -
											(QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]] );
				break;
			case 7:
				surface_enthalpies__[k] = 2.*QStar__[index_CStar__[k]] -
					                     (QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]] );
				break;
			case 8:
				surface_enthalpies__[k] = QStar__[index_CStar__[k]] + QStar__[index_DStar__[k]] - 
										 2.*QStar__[index_AStar__[k]];
				break;
			}

			surface_enthalpies__[k] += dissociation_energies__[k] * R_times_T;
		}
	}

	void UBIQEP_SubMechanism::CalculateActivationEnergies()
	{
		for (unsigned int k = 0; k < half_number_of_ubiqep_reactions_; k++)
		{
			switch (ubiqep_class__[k])
			{
			case 1:
				E_forward__[k] = 0.;
				break;
			case 2:
				E_forward__[k] = 0.;
				break;
			case 3:
				E_forward__[k] = sigma_direct__[k] * (surface_enthalpies__[k] +
					0.5*QStar__[index_AStar__[k]] - QStar__[index_A2_Chemisorption_Heats__[k]]);
				break;
			case 4:
				E_forward__[k] =	sigma_direct__[k] * (surface_enthalpies__[k] +
									QStar__[index_AStar__[k]] * QStar__[index_BStar__[k]] / (QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]]) -
									QStar__[index_AB_Chemisorption_Heats__[k]]);
				break;
			case 5:
				E_forward__[k] = sigma_direct__[k] * (surface_enthalpies__[k] +
					QStar__[index_AStar__[k]] * QStar__[index_BStar__[k]] / (QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]]));
				break;
			case 6:
				E_forward__[k] =	sigma_direct__[k] * (surface_enthalpies__[k] +
									QStar__[index_AStar__[k]] * QStar__[index_BStar__[k]] / (QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]]));
				break;
			case 7:
				E_forward__[k] =	sigma_direct__[k] * (surface_enthalpies__[k] +
									QStar__[index_AStar__[k]] * QStar__[index_BStar__[k]] / (QStar__[index_AStar__[k]] + QStar__[index_BStar__[k]]));
				break;
			case 8:
				E_forward__[k] = sigma_direct__[k] * (surface_enthalpies__[k] + 0.5*QStar__[index_AStar__[k]]);
				break;
			}

			E_backward__[k] = E_forward__[k] - surface_enthalpies__[k];

			// Activation energy check
			if (E_forward__[k] < 0.)
			{
				E_forward__[k] = 0.;
				E_backward__[k] = std::fabs(surface_enthalpies__[k]);
			}
			else if (E_backward__[k] < 0.)
			{
				E_backward__[k] = 0.;
				E_forward__[k] = std::fabs(surface_enthalpies__[k]);
			}
		}
	}

	void UBIQEP_SubMechanism::ForwardKineticConstants(const double T, const double total_site_density, double* kForward)
	{
		const double R_times_T = PhysicalConstants::R_kcal_mol * T;
		const double ln_T_over_T0 = log(T/reference_temperature_);
		const double ln_density = log(total_site_density);
		const double sqrt_T = sqrt(T);

		unsigned int j=0;
		for (unsigned int k=0;k<half_number_of_ubiqep_reactions_;k++)
		{
			if (ubiqep_type_direct__[k] == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
				kForward[j]   = std::exp( lnA_direct__[k]  -lambda__[k]*ln_density + Beta_direct__[k]*ln_T_over_T0 - E_forward__[k]/R_times_T ) * adsorption_coefficient__[k] * sqrt_T;
			else
				kForward[j]   = std::exp( lnA_direct__[k]  + (1.-lambda__[k])*ln_density + Beta_direct__[k]*ln_T_over_T0 - E_forward__[k]/R_times_T );

			if (ubiqep_type_reverse__[k] == PhysicalConstants::UBIQEP_TYPE_ADSORPTION)
				kForward[j+1] = std::exp( lnA_reverse__[k] -lambda__[k]*ln_density + Beta_reverse__[k]*ln_T_over_T0 - E_backward__[k]/R_times_T ) * adsorption_coefficient__[k] * sqrt_T;
			else
				kForward[j+1] = std::exp( lnA_reverse__[k] + (1.-lambda__[k])*ln_density + Beta_reverse__[k]*ln_T_over_T0 - E_backward__[k]/R_times_T );

			j+=2;
		}
	}
}

