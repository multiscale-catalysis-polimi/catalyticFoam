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
|	License                                                               |
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

#include "math/OpenSMOKEStdInclude.h"
#include "math/PhysicalConstants.h"
#include "kernel/thermo/AtomicCompositionTable.h"
#include <boost/algorithm/string.hpp>


namespace OpenSMOKE
{

	template<typename Reactions>
	PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::PreProcessorSolidKineticsPolicy_CHEMKIN() {
	}

	template<typename Reactions>
	PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::PreProcessorSolidKineticsPolicy_CHEMKIN(const PreProcessorSolidKineticsPolicy_CHEMKIN& orig) {
	}

	template<typename Reactions>
	PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::~PreProcessorSolidKineticsPolicy_CHEMKIN() {
	}

	template<typename Reactions>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::ReadFromASCIIFile(const std::string file_name)
	{
		std::cout << " * Reading kinetic file... " << std::endl;
		myKinetics = new InputFileCHEMKIN(file_name);
		
		for (unsigned int j=0;j<myKinetics->good_lines().size();j++)
		{
			size_t found=myKinetics->good_lines()[j].find("MATERIAL");
			if (found!=std::string::npos)
				iMaterialLines.push_back(j);
			found=myKinetics->good_lines()[j].find("material");
			if (found!=std::string::npos)
				iMaterialLines.push_back(j);
		}

		number_of_materials_ = boost::lexical_cast<unsigned int>(iMaterialLines.size());
		iReactionLines.resize(number_of_materials_);
		names_of_materials_.resize(number_of_materials_);

		for (unsigned int j=0;j<number_of_materials_;j++)
		{
			names_of_materials_[j] = myKinetics->good_lines()[iMaterialLines[j]];

			char chars[] = "/ ";
			for (unsigned int i=0; i<strlen(chars); ++i)
				names_of_materials_[j].erase (std::remove(names_of_materials_[j].begin(), names_of_materials_[j].end(), chars[i]), names_of_materials_[j].end());
			names_of_materials_[j].erase(0,8);

			if (names_of_materials_[j].size() == 0)
			{
				std::stringstream jstr; jstr << j+1;
				names_of_materials_[j] = "MATERIAL" + jstr.str();
			}

			// Check for names (TODO)
		}

		iMaterialLines.push_back(boost::lexical_cast<int>(myKinetics->good_lines().size()));

		solid_species_.resize(number_of_materials_);
		solid_species_indices_.resize(number_of_materials_);
		for (unsigned int k=0;k<number_of_materials_;k++)
		{
			std::vector<unsigned int> iSolidLines;
			std::vector<unsigned int> iEndLines;
			unsigned int iAdditionalDataLine = 0;

			for (unsigned int j=iMaterialLines[k]+1;j<iMaterialLines[k+1];j++)
			{
				std::string line = myKinetics->good_lines()[j];

				boost::replace_all(line, "SOLID", "SOLID ");
				boost::replace_all(line, "solid", "SOLID ");
				boost::replace_all(line, "end", "END");
				boost::replace_all(line, "reactions", "REACTIONS ");
				boost::replace_all(line, "additional_data", "ADDITIONAL_DATA ");

				{
					size_t found=line.find("SOLID");
					if (found!=std::string::npos)	iSolidLines.push_back(j);
				}
				{
					size_t found=line.find("END");
					if (found!=std::string::npos)	iEndLines.push_back(j);
				}
				{
					size_t found=line.find("REACTIONS");
					if (found!=std::string::npos)	iReactionLines[k].push_back(j);
				}
				{
					size_t found=line.find("ADDITIONAL_DATA");
					if (found!=std::string::npos)	iAdditionalDataLine = j;
				}
			}
			
			if (iSolidLines.size() == 0)
			{
				std::cout << "Reading solid kinetic scheme. No SOLID section found!" << std::endl;
				return false;
			}
			if (iReactionLines[k].size() == 0)
			{
				std::cout << "Reading surface kinetic scheme. The REACTIONS section is compulsory" << std::endl;
				return false;
			}
			
			if (iReactionLines[k].size() > 1)
			{
				std::cout << "Reading surface kinetic scheme. Only one REACTIONS section must be defined for each material" << std::endl;
				return false;
			}
			
			if (iReactionLines[k][0] <= iSolidLines[iSolidLines.size()-1])
			{
				std::cout << "Reading surface kinetic scheme. The REACTIONS section must be defined only after all the SITE sections" << std::endl;
				return false;
			}

			std::vector<std::string> solid_section(iSolidLines.size());
			for (unsigned int j=0;j<iSolidLines.size();j++)
			{
				unsigned int end_line_index = iReactionLines[k][0];
				if (j!=iSolidLines.size()-1)
					end_line_index = iSolidLines[j+1];

				for (unsigned int i=iSolidLines[j];i<end_line_index;i++)
					solid_section[j] += myKinetics->good_lines()[i] + " ";
			}

			// Solids
			{
				//solid_species[k].resize(iSolidLines.size());
				for (unsigned int j=0;j<iSolidLines.size();j++)
				{
					std::string line = solid_section[j];

					typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
					boost::char_separator<char> sep_blank(" ");
					tokenizer_blank tokens(line, sep_blank);
					std::vector<std::string> separate_elements;
					for (tokenizer_blank::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
						separate_elements.push_back(*tok_iter);

					unsigned int i = 1;
					for(;;)
					{
						if (separate_elements[i] == "END")
						{
							i++;
						}
						else
						{
							solid_species_[k].push_back( separate_elements[i] );
							i++;
						}

						if (i>=separate_elements.size())
							break;
					}
				}

				for(unsigned int ii=0;ii<solid_species_[k].size();ii++)
				{
					const std::size_t n = std::count (solid_species_[k].begin(), solid_species_[k].end(), solid_species_[k][ii]);
					if (n>1)
					{
							std::cout << "The following solid species was specified more than once in the same solid material: " << solid_species_[k][ii] << std::endl;
							return false;
					}
				}
				
				
				for(unsigned int i=0;i<solid_species_[k].size();i++)
				{
					bool iFound = false;
					for(unsigned int ii=0;ii<names_solid_species_.size();ii++)
						if (solid_species_[k][i] == names_solid_species_[ii])
						{
							iFound = true;
							break;
						}

					if (iFound == false)
					{
						names_solid_species_.push_back(solid_species_[k][i]);
					}
				}

				solid_species_indices_[k].resize(solid_species_[k].size());
				for(unsigned int i=0;i<solid_species_[k].size();i++)
					for(unsigned int ii=0;ii<names_solid_species_.size();ii++)
						if (solid_species_[k][i] == names_solid_species_[ii])
							solid_species_indices_[k][i] = ii+1;
			}
		}

		return true;
	}

	template<typename Reactions>
	template<typename PreProcessor>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::KineticsFromASCIIFile(	AtomicCompositionTable& atomicComposition, const PreProcessor& preprocessor_kinetics, std::ostream& flog)
	{
		// Populate the maps
		{
			// Populate the map
			{
				unsigned int count = 0;
				for(unsigned int i=0;i<preprocessor_kinetics.names_species().size();i++)
					map_of_species.insert(std::make_pair(preprocessor_kinetics.names_species()[i],  count++));
				for(unsigned int i=0;i<names_solid_species_.size();i++)
					map_of_species.insert(std::make_pair(names_solid_species_[i],  count++));
			}

			// Populate the global list
			{
				unsigned int count = 0;
				for(unsigned int i=0;i<preprocessor_kinetics.names_species().size();i++)
					names_species_.push_back(preprocessor_kinetics.names_species()[i]);
				for(unsigned int i=0;i<names_solid_species_.size();i++)
					names_species_.push_back(names_solid_species_[i]);
			}

			// Checking the names of the species
			{
				for(unsigned int i=0;i<names_species_.size();i++)
					for(unsigned int j=i+1;j<names_species_.size();j++)
						if (names_species_[i] == names_species_[j])
						{
							std::cout << "The same name was used for more than one species: " << names_species_[j] << std::endl;
							return false;
						}
			}
		}

		reaction_lines.resize(number_of_materials_);
		reactions_.resize(number_of_materials_);

		for(unsigned int k=0;k<number_of_materials_;k++)
		{
			std::cout << "Reading solid-phase Material " << k+1 << " over " << number_of_materials_ << std::endl;
			// Reaction units
			PhysicalConstants::UNITS_REACTION e_units = PhysicalConstants::UNITS_CAL_MOLE;
			PhysicalConstants::UNITS_REACTION a_units = PhysicalConstants::UNITS_MOLES;
			PhysicalConstants::UNITS_REACTION_COMPOSITION composition_units = PhysicalConstants::UNITS_STD;
			{
				typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
				boost::char_separator<char> sep_blank(" ");
				tokenizer_blank tokens(myKinetics->good_lines()[iReactionLines[k][0]], sep_blank);
				const std::size_t n = std::distance (tokens.begin(), tokens.end());
				if (n>6)
				{
					std::cout << "Reading surface kinetic scheme. Wrong number of arguments on the REACTIONS line (line " << myKinetics->indices_of_good_lines()[iReactionLines[k][0]] << ")!" << std::endl;
					return false;
				}

				for (std::size_t i = 2; i <= n; i++)
				{
					std::string tag = *(++tokens.begin());
					
					if (tag == "CAL" || tag == "CAL/MOLE" || tag == "cal" || tag == "cal/mole")
						e_units = PhysicalConstants::UNITS_CAL_MOLE;
					else if (tag == "EVOL" || tag == "EVOLTS" || tag == "evol" || tag == "evolts")
						e_units = PhysicalConstants::UNITS_EVOLTS;
					else if (tag == "JOUL" || tag == "JOULES/MOLE" || tag == "joul" || tag == "joules/mole")
						e_units = PhysicalConstants::UNITS_JOULES_MOLE;
					else if (tag == "KCAL" || tag == "KCAL/MOLE" || tag == "kcal" || tag == "kcal/mole")
						e_units = PhysicalConstants::UNITS_KCAL_MOLE;
					else if (tag == "KJOU" || tag == "KJOULES/MOLE" || tag == "kjou" || tag == "kjoules/mole")
						e_units = PhysicalConstants::UNITS_KJOULES_MOLE;
					else if (tag == "KELV" || tag == "KELVINS" || tag == "kelv" || tag == "kelvins")
						e_units = PhysicalConstants::UNITS_KELVINS;

					else if (tag == "ATM" || tag == "atm")
						composition_units = PhysicalConstants::UNITS_ATM;
					else if (tag == "BAR" || tag == "bar")
						composition_units = PhysicalConstants::UNITS_BAR;
					else if (tag == "DYN" || tag == "dyn" || tag == "DYNES" || tag == "dynes")
						composition_units = PhysicalConstants::UNITS_DYNES;
					else if (tag == "PAS" || tag == "pas" || tag == "PASCALS" || tag == "pascals")
						composition_units = PhysicalConstants::UNITS_PASCALS;
					else if (tag == "TOR" || tag == "tor" || tag == "TORR" || tag == "torr")
						composition_units = PhysicalConstants::UNITS_TORR;
					
					else if (tag == "MOLEC" || tag == "MOLECULES" || tag == "molec" || tag == "molecules")
						a_units = PhysicalConstants::UNITS_MOLECULES;
					else if (tag == "MOLE" || tag == "MOLES" || tag == "mole" || tag == "moles")
						a_units = PhysicalConstants::UNITS_MOLES;

					else
					{
						std::cout << "Reading surface kinetic scheme. Wrong units on the REACTION line (line " << myKinetics->indices_of_good_lines()[iReactionLines[k][0]] << ")!" << std::endl;
						std::cout << "Available options: CAL || CAL/MOLE || EVOL || EVOLTS || JOUL || JOULES/MOLE || KCAL || KCAL/MOLE || KELV || KELVINS || MOLEC || MOLECULES || MOLE || MOLES || ATM || BAR || DYN || DYNES || PAS || PASCALS || TOR || TORR" << std::endl;
						return false;
					}
				}
			}
			
			unsigned int iEndReactionLine = 0;
			for (unsigned int j=iReactionLines[k][0];j<myKinetics->good_lines().size();j++)
			{
				size_t found_END=myKinetics->good_lines()[j].find("END");
				if (found_END!=std::string::npos)
				{
					iEndReactionLine = j+1;
					break;
				}
				else
				{
					size_t found_end=myKinetics->good_lines()[j].find("end");
					if (found_end!=std::string::npos)
					{
						iEndReactionLine = j+1;
						break;
					}
				}

				size_t found_MATERIAL=myKinetics->good_lines()[j].find("MATERIAL");
				if (found_MATERIAL!=std::string::npos)
				{
					std::cout << "Reading surface kinetic scheme. Missing END (end) keyword at the end of the list of the reactions!" << std::endl;
					return false;
				}
				size_t found_material=myKinetics->good_lines()[j].find("material");
				if (found_material!=std::string::npos)
				{
					std::cout << "Reading surface kinetic scheme. Missing END (end) keyword at the end of the list of the reactions!" << std::endl;
					return false;
				}
				
			}
			if (iEndReactionLine == 0)
			{
				std::cout << "Reading surface kinetic scheme. Missing END (end) keyword at the end of the list of the reactions!" << std::endl;
				return false;
			}
			
			unsigned int number_of_reactions = 0;
			for (unsigned int j=iReactionLines[k][0];j<iEndReactionLine;j++)
			{
				size_t found=myKinetics->good_lines()[j].find("=");
				if (found!=std::string::npos)
				{
					reaction_lines[k].push_back(j);
					number_of_reactions++;
				}
			}

			reaction_lines[k].push_back(iEndReactionLine-1);

			reactions_[k].resize(number_of_reactions);
			if (number_of_reactions > 20)
				std::cout << " * Parsing " << number_of_reactions << " reactions: ";

			unsigned int large_error_in_stoichiometries = 0;
			unsigned int small_error_in_stoichiometries = 0;
			for (unsigned int j=0;j<number_of_reactions;j++)
			{
			
				if (number_of_reactions > 5)
				{
					if (j%(number_of_reactions/5) == 0) 
						std::cout << "%";
					if (j==number_of_reactions-1)
						std::cout << std::endl;
				}

				std::vector<std::string> list_of_lines;
				for (unsigned int i=reaction_lines[k][j];i<reaction_lines[k][j+1];i++)
					list_of_lines.push_back(myKinetics->good_lines()[i]);
				try
				{
					reactions_[k][j].SetDefaultUnits();
					reactions_[k][j].SetUnits(a_units, e_units, composition_units);
			
					bool successReading = reactions_[k][j].ReadReactionFromStrings(	list_of_lines, map_of_species, 
																				    boost::lexical_cast<unsigned int>(preprocessor_kinetics.names_species().size()),
																					boost::lexical_cast<unsigned int>(names_solid_species_.size()));
					if (successReading == false)
						throw reaction_lines[k][j];
				
					unsigned int successStoichiometry = atomicComposition.CheckStoichiometry(flog, reactions_[k][j], 1e-3);
					if (successStoichiometry > 0)
					{
						if (successStoichiometry == 1)
						{
							large_error_in_stoichiometries++;
							std::string reaction_string; reactions_[k][j].GetReactionString(names_species_,reaction_string);
							boost::erase_all(reaction_string, " ");
							flog << "Error in reaction (line " << reaction_lines[k][j]+1 << "): " << reaction_string << std::endl;
						}
						if (successStoichiometry == 2)
						{
							small_error_in_stoichiometries++;
							std::string reaction_string; reactions_[k][j].GetReactionString(names_species_,reaction_string);
							boost::erase_all(reaction_string, " ");
							flog << "Warning in reaction (line " << reaction_lines[k][j]+1 << "): " << reaction_string << std::endl;					
						}					
						flog << std::endl;
					}
				}
				catch(unsigned int kk)
				{
					std::cout << "Reading kinetic scheme: error in reaction starting at line " << myKinetics->indices_of_good_lines()[kk] << std::endl;
					return false;
				}
			}
			
			if (large_error_in_stoichiometries > 0)
			{
				std::cout << std::endl;
				std::cout << " ! ERROR MESSAGE: Large inconsistencies were found in the stoichiometries of " << large_error_in_stoichiometries << " reactions." << std::endl;
				std::cout << "                  Please check the log file for additional details." << std::endl;
				std::cout << std::endl;
				return false;
			}

			if (small_error_in_stoichiometries > 0)
			{
				std::cout << std::endl;
				std::cout << " ! WARNING MESSAGE: Small inconsistencies were found in the stoichiometries of " << small_error_in_stoichiometries << " reactions." << std::endl;
				std::cout << "                    Please check the log file for additional details." << std::endl;
				std::cout << std::endl;
			}
			
			// Check for duplicate reactions
			{
				std::cout << " * Looking for duplicate reactions... " << std::endl;
			
				for (unsigned int i=0; i<reactions_[k].size(); i++)
					for (unsigned int j=i+1; j<reactions_[k].size(); j++)
					{
						if ( OpenSMOKE_Utilities::compare_vectors( reactions_[k][i].reactant_nu_indices(), reactions_[k][j].reactant_nu_indices() ) == true)
						{
							if ( OpenSMOKE_Utilities::compare_vectors( reactions_[k][i].product_nu_indices(), reactions_[k][j].product_nu_indices() ) == true)
							{
								if ( reactions_[k][i].Tag() == reactions_[k][j].Tag() )
								{
									if (reactions_[k][i].IsDuplicate() == false || reactions_[k][j].IsDuplicate() == false )
									{					
											std::cout << "The following reactions must be declared as DUPLICATE" << std::endl;
											std::cout << "Reaction " << i+1 << " starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[k][i]] << std::endl;
											std::cout << "Reaction " << j+1 << " starting at line: " << myKinetics->indices_of_good_lines()[reaction_lines[k][j]] << std::endl;
											return false;
									}
								}
							}
						}
					}
			}
		}

		return true;
	}

	template<typename Reactions>
	template<typename Thermodynamics>
	void PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteShortSummaryOnASCIIFile(const std::string file_name, Thermodynamics& thermodynamics) const
	{
		std::cout << " * Writing the summary of solid-phase kinetic mechanism..." << std::endl;

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		thermodynamics.WriteElementTableOnASCIIFile(fOutput);

		const std::size_t number_of_gas_species = names_species_.size() - names_solid_species_.size();

		fOutput << "---------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                                SOLID MECHANISM                                      " << std::endl;
		fOutput << "---------------------------------------------------------------------------------------" << std::endl;

		fOutput << std::endl;
		fOutput << std::setw(7) << " "; fOutput << std::left << "Number of materials:       " << number_of_materials_ << std::endl;
		fOutput << std::setw(7) << " "; fOutput << std::left << "Number of solid species:   " << names_solid_species_.size() << std::endl;


		fOutput << std::endl;
		fOutput << std::setw(7) << " "; fOutput << std::left << "Materials: " << std::endl;
		for(unsigned int k=0;k<number_of_materials_;k++)
		{
			fOutput << std::setw(7) << " ";
			fOutput << std::left << names_of_materials_[k] << std::endl;
		}
		
		if (names_solid_species_.size() != 0)
		{
			for(unsigned int k=0;k<number_of_materials_;k++)
			{
				fOutput << std::endl;
				fOutput << std::setw(7) << " ";
				fOutput << std::left << "------------------------------------------------------------------------------------------------------------------" << std::endl;
				fOutput << std::setw(7) << " ";
				fOutput << std::setw(24) << std::left << names_of_materials_[k];
				fOutput << std::setw(16) << std::left << "global-index";
				fOutput << std::setw(16) << std::left << "solid-index";
				fOutput << std::setw(16) << std::left << "local-index";
				fOutput << std::endl;
				fOutput << std::setw(7) << " ";
				fOutput << std::left << "------------------------------------------------------------------------------------------------------------------" << std::endl;
		
				for(unsigned int i=0;i<solid_species_[k].size();i++)
				{
					fOutput << std::setw(7)  << " "; 
					fOutput << std::setw(24) << std::left << solid_species_[k][i];
					fOutput << std::fixed   << std::setw(16) << std::setprecision(0) << std::left << number_of_gas_species + solid_species_indices_[k][i];
					fOutput << std::fixed   << std::setw(16) << std::setprecision(0) << std::left << solid_species_indices_[k][i];
					fOutput << std::fixed   << std::setw(16) << std::setprecision(0) << std::left << i+1;
					fOutput << std::endl;
				}
			}
		}

		fOutput << std::endl;
		fOutput << "--------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                                  CHEMICAL REACTIONS                                                  " << std::endl;
		fOutput << std::endl;
		fOutput << "                         Units: [kmol, m3, s] and [cal/mol]                                           " << std::endl;
		fOutput << "--------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
		fOutput << std::endl;

		unsigned int count = 1;
		for(unsigned int kk=0;kk<number_of_materials_;kk++)
		{
			for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
			{
				fOutput << std::setw(7) << std::right << count++;
				fOutput << ". ";
				(*it).WriteShortSummary(fOutput, names_species_);
				fOutput << std::endl;
				fOutput << std::endl;
			}
		}

		fOutput.close();
	}

/*	template<typename T>
	void WriteObjectASCIIFileOldStyle(const T& v, std::ostream& fOut)
	{
		fOut << v.size() << std::endl;
		for(int i=0;i<v.size();i++)
			fOut << 	v(i) << " ";
		fOut << std::endl;
	}
	*/
	/*
	template<typename Reactions>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteASCIIFileOldStyle(const std::string file_name) const
	{
		std::cout << " * Writing the interpreted kinetic file in ASCII format for old versions of OpenSMOKE..." << std::endl;

		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			if ( (*it).IsExplicitlyReversible()==true)
				ErrorMessage("PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteASCIIFileOldStyle(const std::string file_name)",
							 "The kinetic mechanism contains one or more REV reactions. Please write them as separate reactions to continue...");

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		// 1. Number of species
		fOutput << names_species_.size() << std::endl;

		// 2. Number of reactions
		fOutput << reactions_.size() << std::endl;

		// 3. Kinetic data for each reaction
		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			(*it).WriteKineticsDataOnASCIIFileOldStyle(fOutput);
	
		// 4. Stoichiometric data
		WriteStoichiometricDataOnASCIIFile(fOutput);

		// 6. Writing reaction strings
		{
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				std::string line_reaction;
				(*it).GetReactionString(names_species_, line_reaction);
				boost::erase_all(line_reaction, " ");
				fOutput << line_reaction << std::endl; 
			}
		}

		// 7. Reaction orders
		{
			Eigen::VectorXd forwardOrders(reactions_.size()); 
			Eigen::VectorXd backwardOrders(reactions_.size()); 
			forwardOrders.setConstant(0.);
			backwardOrders.setConstant(0.);
			for (unsigned int k=0; k<reactions_.size(); k++)
			{
				forwardOrders(k) = reactions_[k].sumLambdaReactants();
				backwardOrders(k) = reactions_[k].sumLambdaProducts();
			}

			WriteObjectASCIIFileOldStyle(forwardOrders, fOutput);
			WriteObjectASCIIFileOldStyle(backwardOrders, fOutput);
		}

		// 8. Old data (TODO reaction orders different than stoichiometric coefficients)
		{
			fOutput << 0 << std::endl;
			fOutput << 0 << std::endl;
		}

		fOutput.close();




		return true;
	}
	*/
	
	template<typename Reactions>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteStoichiometricDataOnASCIIFile(std::ostream& fOutput, const unsigned int kk) const	
	{
		{
			Eigen::VectorXi numDir1(names_species_.size()); numDir1.setConstant(0);
			Eigen::VectorXi numDir2(names_species_.size()); numDir2.setConstant(0);
			Eigen::VectorXi numDir3(names_species_.size()); numDir3.setConstant(0);
			Eigen::VectorXi numDir4(names_species_.size()); numDir4.setConstant(0);
			Eigen::VectorXi numDir5(names_species_.size()); numDir5.setConstant(0);

			Eigen::VectorXi jDir1;
			Eigen::VectorXi jDir2;
			Eigen::VectorXi jDir3;
			Eigen::VectorXi jDir4;
			Eigen::VectorXi jDir5;
			Eigen::VectorXd vDir5;

			Eigen::VectorXi numInvTot1(names_species_.size()); numInvTot1.setConstant(0);
			Eigen::VectorXi numInvTot2(names_species_.size()); numInvTot2.setConstant(0);
			Eigen::VectorXi numInvTot3(names_species_.size()); numInvTot3.setConstant(0);
			Eigen::VectorXi numInvTot4(names_species_.size()); numInvTot4.setConstant(0);
			Eigen::VectorXi numInvTot5(names_species_.size()); numInvTot5.setConstant(0);

			Eigen::VectorXi jInvTot1;
			Eigen::VectorXi jInvTot2;
			Eigen::VectorXi jInvTot3;
			Eigen::VectorXi jInvTot4;
			Eigen::VectorXi jInvTot5;
			Eigen::VectorXd vInvTot5;

			Eigen::VectorXi numInvEq1(names_species_.size()); numInvEq1.setConstant(0);
			Eigen::VectorXi numInvEq2(names_species_.size()); numInvEq2.setConstant(0);
			Eigen::VectorXi numInvEq3(names_species_.size()); numInvEq3.setConstant(0);
			Eigen::VectorXi numInvEq4(names_species_.size()); numInvEq4.setConstant(0);
			Eigen::VectorXi numInvEq5(names_species_.size()); numInvEq5.setConstant(0);

			Eigen::VectorXi jInvEq1;
			Eigen::VectorXi jInvEq2;
			Eigen::VectorXi jInvEq3;
			Eigen::VectorXi jInvEq4;
			Eigen::VectorXi jInvEq5;
			Eigen::VectorXd vInvEq5;
			
			{
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
				{
					for (unsigned int i=0;i<(*it).reactant_nu_indices().size();i++)
					{
						if ( (*it).reactant_nu()[i] == 1. )			numDir1((*it).reactant_nu_indices()[i])++;
						else if ( (*it).reactant_nu()[i] == 2.  )	numDir2((*it).reactant_nu_indices()[i])++;
						else if ( (*it).reactant_nu()[i] == 3.  )	numDir3((*it).reactant_nu_indices()[i])++;
						else if ( (*it).reactant_nu()[i] == 0.5 )	numDir4((*it).reactant_nu_indices()[i])++;
						else 										numDir5((*it).reactant_nu_indices()[i])++;
					}
				}

				jDir1.resize(numDir1.sum()); jDir1.setConstant(0);
				jDir2.resize(numDir2.sum()); jDir2.setConstant(0);
				jDir3.resize(numDir3.sum()); jDir3.setConstant(0);
				jDir4.resize(numDir4.sum()); jDir4.setConstant(0);
				jDir5.resize(numDir5.sum()); jDir5.setConstant(0);
				vDir5.resize(numDir5.sum()); vDir5.setConstant(0.);

				const unsigned int sumDirTot = numDir1.sum()+numDir2.sum()+numDir3.sum()+numDir4.sum()+numDir5.sum();

				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					unsigned int k = -1;
					for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
					{
						k++;

						for (unsigned int i=0;i<(*it).reactant_nu_indices().size();i++)
						{
							if ( (*it).reactant_nu_indices()[i] == j )
							{
								if ((*it).reactant_nu()[i] == 1.)
								{
									jDir1(countGlobal1++) = k+1;
								}
								else if ((*it).reactant_nu()[i] == 2.)
								{
									jDir2(countGlobal2++) = k+1;
								}
								else if ((*it).reactant_nu()[i] == 3.)
								{
									jDir3(countGlobal3++) = k+1;
								}
								else if ((*it).reactant_nu()[i] == 0.5)
								{
									jDir4(countGlobal4++) = k+1;
								}
								else
								{
									jDir5(countGlobal5) = k+1;
									vDir5(countGlobal5) = (*it).reactant_nu()[i];
									countGlobal5++;
								}
								countGlobal++;
							}
						}
						if (countGlobal == sumDirTot)
							break;
					}
				}
			}
			
			
			{
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
				{
					for (unsigned int i=0;i<(*it).product_nu_indices().size();i++)
					{
						if ( (*it).product_nu()[i] == 1. )			numInvTot1((*it).product_nu_indices()[i])++;
						else if ( (*it).product_nu()[i] == 2.  )	numInvTot2((*it).product_nu_indices()[i])++;
						else if ( (*it).product_nu()[i] == 3.  )	numInvTot3((*it).product_nu_indices()[i])++;
						else if ( (*it).product_nu()[i] == 0.5 )	numInvTot4((*it).product_nu_indices()[i])++;
						else 										numInvTot5((*it).product_nu_indices()[i])++;
					}
				}

				jInvTot1.resize(numInvTot1.sum()); jInvTot1.setConstant(0);
				jInvTot2.resize(numInvTot2.sum()); jInvTot2.setConstant(0);
				jInvTot3.resize(numInvTot3.sum()); jInvTot3.setConstant(0);
				jInvTot4.resize(numInvTot4.sum()); jInvTot4.setConstant(0);
				jInvTot5.resize(numInvTot5.sum()); jInvTot5.setConstant(0);
				vInvTot5.resize(numInvTot5.sum()); vInvTot5.setConstant(0.);

				const unsigned int sumInvTot = numInvTot1.sum()+numInvTot2.sum()+numInvTot3.sum()+numInvTot4.sum()+numInvTot5.sum();
		
				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					unsigned int k = -1;
					for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
					{
						k++;
						for (unsigned int i=0;i<(*it).product_nu_indices().size();i++)
						{
							if ( (*it).product_nu_indices()[i] == j )
							{
								if ((*it).product_nu()[i] == 1.)
								{
									jInvTot1(countGlobal1++) = k+1;
								}
								else if ((*it).product_nu()[i] == 2.)
								{
									jInvTot2(countGlobal2++) = k+1;
								}
								else if ((*it).product_nu()[i] == 3.)
								{
									jInvTot3(countGlobal3++) = k+1;
								}
								else if ((*it).product_nu()[i] == 0.5)
								{
									jInvTot4(countGlobal4++) = k+1;
								}
								else
								{
									jInvTot5(countGlobal5) = k+1;
									vInvTot5(countGlobal5) = (*it).product_nu()[i];
									countGlobal5++;
								}
								countGlobal++;
							}
						}
						if (countGlobal == sumInvTot)
							break;
					}
				}
			}
			
			{
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
				{
					if ((*it).IsReversible() == true)
					{
						for (unsigned int i=0;i<(*it).product_nu_indices().size();i++)
						{
							if ( (*it).product_nu()[i] == 1. )			numInvEq1((*it).product_nu_indices()[i])++;
							else if ( (*it).product_nu()[i] == 2.  )	numInvEq2((*it).product_nu_indices()[i])++;
							else if ( (*it).product_nu()[i] == 3.  )	numInvEq3((*it).product_nu_indices()[i])++;
							else if ( (*it).product_nu()[i] == 0.5 )	numInvEq4((*it).product_nu_indices()[i])++;
							else 										numInvEq5((*it).product_nu_indices()[i])++;
						}
					}
				}

				jInvEq1.resize(numInvEq1.sum()); jInvEq1.setConstant(0);
				jInvEq2.resize(numInvEq2.sum()); jInvEq2.setConstant(0);
				jInvEq3.resize(numInvEq3.sum()); jInvEq3.setConstant(0);
				jInvEq4.resize(numInvEq4.sum()); jInvEq4.setConstant(0);
				jInvEq5.resize(numInvEq5.sum()); jInvEq5.setConstant(0);
				vInvEq5.resize(numInvEq5.sum()); vInvEq5.setConstant(0.);

				const unsigned int sumInvEq = numInvEq1.sum()+numInvEq2.sum()+numInvEq3.sum()+numInvEq4.sum()+numInvEq5.sum();
		
				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					unsigned int k = -1;
					for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
					{
						k++;
						if ((*it).IsReversible() == true)
						{
							for (unsigned int i=0;i<(*it).product_nu_indices().size();i++)
							{
								if ( (*it).product_nu_indices()[i] == j )
								{
									if ((*it).product_nu()[i] == 1.)
									{
										jInvEq1(countGlobal1++) = k+1;
									}
									else if ((*it).product_nu()[i] == 2.)
									{
										jInvEq2(countGlobal2++) = k+1;
									}
									else if ((*it).product_nu()[i] == 3.)
									{
										jInvEq3(countGlobal3++) = k+1;
									}
									else if ((*it).product_nu()[i] == 0.5)
									{
										jInvEq4(countGlobal4++) = k+1;
									}
									else
									{
										jInvEq5(countGlobal5) = k+1;
										vInvEq5(countGlobal5) = (*it).product_nu()[i];
										countGlobal5++;
									}
									countGlobal++;
								}
							}
						}
						if (countGlobal == sumInvEq)
							break;
					}
				}
			}
			
			
			WriteObjectASCIIFileOldStyle(numDir1, fOutput);
			WriteObjectASCIIFileOldStyle(numDir2, fOutput);
			WriteObjectASCIIFileOldStyle(numDir3, fOutput);
			WriteObjectASCIIFileOldStyle(numDir4, fOutput);
			WriteObjectASCIIFileOldStyle(numDir5, fOutput);
			
			WriteObjectASCIIFileOldStyle(numInvTot1, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot2, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot3, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot4, fOutput);
			WriteObjectASCIIFileOldStyle(numInvTot5, fOutput);

			WriteObjectASCIIFileOldStyle(numInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq5, fOutput);

			WriteObjectASCIIFileOldStyle(jDir1, fOutput);
			WriteObjectASCIIFileOldStyle(jDir2, fOutput);
			WriteObjectASCIIFileOldStyle(jDir3, fOutput);
			WriteObjectASCIIFileOldStyle(jDir4, fOutput);
			WriteObjectASCIIFileOldStyle(jDir5, fOutput);
			WriteObjectASCIIFileOldStyle(vDir5, fOutput);

			WriteObjectASCIIFileOldStyle(jInvTot1, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot2, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot3, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot4, fOutput);
			WriteObjectASCIIFileOldStyle(jInvTot5, fOutput);
			WriteObjectASCIIFileOldStyle(vInvTot5, fOutput);

			WriteObjectASCIIFileOldStyle(jInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq5, fOutput);
			WriteObjectASCIIFileOldStyle(vInvEq5, fOutput);
			
		}
		
		// 5. Sum of stoichiometric coefficients
		{
			Eigen::VectorXd sumStoichiometricCoefficients(reactions_[kk].size()); 
			sumStoichiometricCoefficients.setConstant(0.);
			
			unsigned int k = 0;
			for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
			{
				sumStoichiometricCoefficients(k) = (*it).sumNuProducts()-(*it).sumNuReactants();
				k++;
			}

			WriteObjectASCIIFileOldStyle(sumStoichiometricCoefficients, fOutput);
		}
		
		return true;
	}
	

	
	template<typename Reactions>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteReactionOrdersOnASCIIFile(std::ostream& fOutput, const unsigned int kk) const	
	{
		bool write_reaction_orders = false;
		for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
			if ((*it).IsFORD() == true || (*it).IsRORD() == true)
			{
				write_reaction_orders = true;
				break;
			}

		fOutput << write_reaction_orders << std::endl;

		if (write_reaction_orders == true)
		{
		{
			Eigen::VectorXi numDir1(names_species_.size()); numDir1.setConstant(0);
			Eigen::VectorXi numDir2(names_species_.size()); numDir2.setConstant(0);
			Eigen::VectorXi numDir3(names_species_.size()); numDir3.setConstant(0);
			Eigen::VectorXi numDir4(names_species_.size()); numDir4.setConstant(0);
			Eigen::VectorXi numDir5(names_species_.size()); numDir5.setConstant(0);

			Eigen::VectorXi jDir1;
			Eigen::VectorXi jDir2;
			Eigen::VectorXi jDir3;
			Eigen::VectorXi jDir4;
			Eigen::VectorXi jDir5;
			Eigen::VectorXd vDir5;

			Eigen::VectorXi numInvEq1(names_species_.size()); numInvEq1.setConstant(0);
			Eigen::VectorXi numInvEq2(names_species_.size()); numInvEq2.setConstant(0);
			Eigen::VectorXi numInvEq3(names_species_.size()); numInvEq3.setConstant(0);
			Eigen::VectorXi numInvEq4(names_species_.size()); numInvEq4.setConstant(0);
			Eigen::VectorXi numInvEq5(names_species_.size()); numInvEq5.setConstant(0);

			Eigen::VectorXi jInvEq1;
			Eigen::VectorXi jInvEq2;
			Eigen::VectorXi jInvEq3;
			Eigen::VectorXi jInvEq4;
			Eigen::VectorXi jInvEq5;
			Eigen::VectorXd vInvEq5;

			{
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
				{
					for (unsigned int i=0;i<(*it).reactant_lambda_indices().size();i++)
					{
						if ( (*it).reactant_lambda()[i] == 1. )			numDir1((*it).reactant_lambda_indices()[i])++;
						else if ( (*it).reactant_lambda()[i] == 2.  )	numDir2((*it).reactant_lambda_indices()[i])++;
						else if ( (*it).reactant_lambda()[i] == 3.  )	numDir3((*it).reactant_lambda_indices()[i])++;
						else if ( (*it).reactant_lambda()[i] == 0.5 )	numDir4((*it).reactant_lambda_indices()[i])++;
						else 											numDir5((*it).reactant_lambda_indices()[i])++;
					}
				}

				jDir1.resize(numDir1.sum()); jDir1.setConstant(0);
				jDir2.resize(numDir2.sum()); jDir2.setConstant(0);
				jDir3.resize(numDir3.sum()); jDir3.setConstant(0);
				jDir4.resize(numDir4.sum()); jDir4.setConstant(0);
				jDir5.resize(numDir5.sum()); jDir5.setConstant(0);
				vDir5.resize(numDir5.sum()); vDir5.setConstant(0.);

				const unsigned int sumDirTot = numDir1.sum()+numDir2.sum()+numDir3.sum()+numDir4.sum()+numDir5.sum();

				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					unsigned int k = -1;
					for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
					{
						k++;
						for (unsigned int i=0;i<(*it).reactant_lambda_indices().size();i++)
						{
							if ( (*it).reactant_lambda_indices()[i] == j )
							{
								if ((*it).reactant_lambda()[i] == 1.)
								{
									jDir1(countGlobal1++) = k+1;
								}
								else if ((*it).reactant_lambda()[i] == 2.)
								{
									jDir2(countGlobal2++) = k+1;
								}
								else if ((*it).reactant_lambda()[i] == 3.)
								{
									jDir3(countGlobal3++) = k+1;
								}
								else if ((*it).reactant_lambda()[i] == 0.5)
								{
									jDir4(countGlobal4++) = k+1;
								}
								else
								{
									jDir5(countGlobal5) = k+1;
									vDir5(countGlobal5) =(*it).reactant_lambda()[i];
									countGlobal5++;
								}
								countGlobal++;
							}
						}
						if (countGlobal == sumDirTot)
							break;
					}
				}
			}
			
			{
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
				{
					if ((*it).IsReversible() == true)
					{
						for (unsigned int i=0;i<(*it).product_lambda_indices().size();i++)
						{
							if ( (*it).product_lambda()[i] == 1. )			numInvEq1((*it).product_lambda_indices()[i])++;
							else if ( (*it).product_lambda()[i] == 2.  )	numInvEq2((*it).product_lambda_indices()[i])++;
							else if ( (*it).product_lambda()[i] == 3.  )	numInvEq3((*it).product_lambda_indices()[i])++;
							else if ( (*it).product_lambda()[i] == 0.5 )	numInvEq4((*it).product_lambda_indices()[i])++;
							else 											numInvEq5((*it).product_lambda_indices()[i])++;
						}
					}
				}

				jInvEq1.resize(numInvEq1.sum()); jInvEq1.setConstant(0);
				jInvEq2.resize(numInvEq2.sum()); jInvEq2.setConstant(0);
				jInvEq3.resize(numInvEq3.sum()); jInvEq3.setConstant(0);
				jInvEq4.resize(numInvEq4.sum()); jInvEq4.setConstant(0);
				jInvEq5.resize(numInvEq5.sum()); jInvEq5.setConstant(0);
				vInvEq5.resize(numInvEq5.sum()); vInvEq5.setConstant(0.);

				const unsigned int sumInvEq = numInvEq1.sum()+numInvEq2.sum()+numInvEq3.sum()+numInvEq4.sum()+numInvEq5.sum();
		
				unsigned int countGlobal1 = 0;
				unsigned int countGlobal2 = 0;
				unsigned int countGlobal3 = 0;
				unsigned int countGlobal4 = 0;
				unsigned int countGlobal5 = 0;
				unsigned int countGlobal = 0;
				for (unsigned int j=0;j<names_species_.size();j++)
				{
					unsigned int k = -1;
					for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[kk].begin(); it != reactions_[kk].end(); ++it)
					{
						k++;
						if ((*it).IsReversible() == true)
						{
							for (unsigned int i=0;i<(*it).product_lambda_indices().size();i++)
							{
								if ( (*it).product_lambda_indices()[i] == j )
								{
									if ((*it).product_lambda()[i] == 1.)
									{
										jInvEq1(countGlobal1++) = k+1;
									}
									else if ((*it).product_lambda()[i] == 2.)
									{
										jInvEq2(countGlobal2++) = k+1;
									}
									else if ((*it).product_lambda()[i] == 3.)
									{
										jInvEq3(countGlobal3++) = k+1;
									}
									else if ((*it).product_lambda()[i] == 0.5)
									{
										jInvEq4(countGlobal4++) = k+1;
									}
									else
									{
										jInvEq5(countGlobal5) = k+1;
										vInvEq5(countGlobal5) = (*it).product_lambda()[i];
										countGlobal5++;
									}
									countGlobal++;
								}
							}
						}
						if (countGlobal == sumInvEq)
							break;
					}
				}
			}
			
			WriteObjectASCIIFileOldStyle(numDir1, fOutput);
			WriteObjectASCIIFileOldStyle(numDir2, fOutput);
			WriteObjectASCIIFileOldStyle(numDir3, fOutput);
			WriteObjectASCIIFileOldStyle(numDir4, fOutput);
			WriteObjectASCIIFileOldStyle(numDir5, fOutput);

			WriteObjectASCIIFileOldStyle(numInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(numInvEq5, fOutput);

			WriteObjectASCIIFileOldStyle(jDir1, fOutput);
			WriteObjectASCIIFileOldStyle(jDir2, fOutput);
			WriteObjectASCIIFileOldStyle(jDir3, fOutput);
			WriteObjectASCIIFileOldStyle(jDir4, fOutput);
			WriteObjectASCIIFileOldStyle(jDir5, fOutput);
			WriteObjectASCIIFileOldStyle(vDir5, fOutput);

			WriteObjectASCIIFileOldStyle(jInvEq1, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq2, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq3, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq4, fOutput);
			WriteObjectASCIIFileOldStyle(jInvEq5, fOutput);
			WriteObjectASCIIFileOldStyle(vInvEq5, fOutput);
		}
		}// write_reaction_orders
		
		return true;
	}
	
	/*
	template<typename Reactions>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteASCIIFile(const std::string file_name) const
	{
		std::cout << " * Writing the interpreted kinetic file in ASCII format..." << std::endl;

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out | std::ios::binary);
		fOutput.setf(std::ios::scientific);

		// 1. Number of species
		fOutput << "number-of-species" << std::endl;
		fOutput << names_species_.size() << std::endl;

		// 2. Number of reactions
		fOutput << "number-of-reactions" << std::endl;
		fOutput << reactions_.size() << std::endl;

		// 3. Kinetic data for each reaction
		unsigned int number_of_irreversible_reactions = 0;
		unsigned int number_of_reversible_reactions = 0;
		unsigned int number_of_explicitly_reversible_reactions = 0;
		unsigned int number_of_thermodynamic_reversible_reactions = 0;
		unsigned int number_of_thirdbody_reactions = 0;
		unsigned int number_of_falloff_reactions = 0;
		unsigned int number_of_cabr_reactions = 0;
		unsigned int number_of_chebyshev_reactions = 0;
		
		unsigned int number_of_pressurelog_reactions = 0;
		unsigned int number_of_fit1_reactions = 0;
		unsigned int number_of_janevlanger_reactions = 0;
		unsigned int number_of_landauteller_reactions = 0;

		unsigned int number_of_falloff_lindemann_reactions = 0;
		unsigned int number_of_falloff_troe_reactions = 0;
		unsigned int number_of_falloff_sri_reactions = 0;
		unsigned int number_of_cabr_lindemann_reactions = 0;
		unsigned int number_of_cabr_troe_reactions = 0;
		unsigned int number_of_cabr_sri_reactions = 0;

		for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
		{
			if ( (*it).IsReversible()==true)
			{
				number_of_reversible_reactions++;
				if ( (*it).IsExplicitlyReversible()==true)	number_of_explicitly_reversible_reactions++;
				else										number_of_thermodynamic_reversible_reactions++;
			}
			else
				number_of_irreversible_reactions++;
			
			if ( (*it).Tag() == PhysicalConstants::REACTION_THIRDBODY) number_of_thirdbody_reactions++;
			else if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF) { number_of_falloff_reactions++; number_of_falloff_lindemann_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_FALLOFF) {number_of_falloff_reactions++; number_of_falloff_troe_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_FALLOFF) {number_of_falloff_reactions++; number_of_falloff_sri_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR) {number_of_cabr_reactions++; number_of_cabr_lindemann_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_CABR) {number_of_cabr_reactions++; number_of_cabr_troe_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_CABR) {number_of_cabr_reactions++; number_of_cabr_sri_reactions++; }
			else if ( (*it).Tag() == PhysicalConstants::REACTION_CHEBYSHEV) number_of_chebyshev_reactions++;
			else if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE)
			{
				if ( (*it).IsPressureLog() == true )		number_of_pressurelog_reactions++;
				else if ( (*it).IsJanevLanger() == true )	number_of_janevlanger_reactions++;
				else if ( (*it).IsFit1() == true )			number_of_fit1_reactions++;
				else if ( (*it).IsLandauTeller() == true )	number_of_landauteller_reactions++;
			}
		}

		std::vector<unsigned int> indices_of_irreversible_reactions(number_of_irreversible_reactions);
		std::vector<unsigned int> indices_of_reversible_reactions(number_of_reversible_reactions);
		std::vector<unsigned int> indices_of_explicitly_reversible_reactions(number_of_explicitly_reversible_reactions);
		std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions(number_of_thermodynamic_reversible_reactions);
		std::vector<unsigned int> indices_of_thirdbody_reactions(number_of_thirdbody_reactions);
		std::vector<unsigned int> indices_of_falloff_reactions(number_of_falloff_reactions);
		std::vector<unsigned int> indices_of_cabr_reactions(number_of_cabr_reactions);
		std::vector<unsigned int> indices_of_chebyshev_reactions(number_of_chebyshev_reactions);
		std::vector<unsigned int> indices_of_pressurelog_reactions(number_of_pressurelog_reactions);
		std::vector<unsigned int> indices_of_janevlanger_reactions(number_of_janevlanger_reactions);
		std::vector<unsigned int> indices_of_fit1_reactions(number_of_fit1_reactions);
		std::vector<unsigned int> indices_of_landauteller_reactions(number_of_landauteller_reactions);

		std::vector<unsigned int> indices_of_falloff_lindemann_reactions(number_of_falloff_lindemann_reactions);
		std::vector<unsigned int> indices_of_cabr_lindemann_reactions(number_of_cabr_lindemann_reactions);
		std::vector<unsigned int> indices_of_falloff_troe_reactions(number_of_falloff_troe_reactions);
		std::vector<unsigned int> indices_of_cabr_troe_reactions(number_of_cabr_troe_reactions);
		std::vector<unsigned int> indices_of_falloff_sri_reactions(number_of_falloff_sri_reactions);
		std::vector<unsigned int> indices_of_cabr_sri_reactions(number_of_cabr_sri_reactions);

		{		
			unsigned int count_irreversible=0;
			unsigned int count_reversible=0;
			unsigned int count_explicit=0;
			unsigned int count_thermodynamic=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).IsReversible()==true) 
				{
						indices_of_reversible_reactions[count_reversible++] = it - reactions_.begin() + 1;
						if ( (*it).IsExplicitlyReversible() == true)  indices_of_explicitly_reversible_reactions[count_explicit++] = it - reactions_.begin() + 1;
						else										  indices_of_thermodynamic_reversible_reactions[count_thermodynamic++] = it - reactions_.begin() + 1;
				}
				else
					indices_of_irreversible_reactions[count_irreversible++] = it - reactions_.begin() + 1;
			}

			fOutput << "irreversible-reactions" << std::endl;
			fOutput << reactions_.size() - number_of_reversible_reactions << std::endl;
			for(unsigned int i=0;i<reactions_.size() - number_of_reversible_reactions;i++)
				fOutput << indices_of_irreversible_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "reversible-reactions" << std::endl;
			fOutput << number_of_reversible_reactions << std::endl;
			for(unsigned int i=0;i<number_of_reversible_reactions;i++)
				fOutput << indices_of_reversible_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "thermodynamic-reversible-reactions" << std::endl;
			fOutput << number_of_thermodynamic_reversible_reactions << std::endl;
			for(unsigned int i=0;i<number_of_thermodynamic_reversible_reactions;i++)
				fOutput << indices_of_thermodynamic_reversible_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "explicitly-reversible-reactions" << std::endl;
			fOutput << number_of_explicitly_reversible_reactions << std::endl;
			for(unsigned int i=0;i<number_of_explicitly_reversible_reactions;i++)
				fOutput << indices_of_explicitly_reversible_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_THIRDBODY) indices_of_thirdbody_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "thirdbody-reactions" << std::endl;
			fOutput << number_of_thirdbody_reactions << std::endl;
			for(unsigned int i=0;i<number_of_thirdbody_reactions;i++)
				fOutput << indices_of_thirdbody_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			unsigned int count_lindemann=0;
			unsigned int count_troe=0;
			unsigned int count_sri=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_falloff_lindemann_reactions[count_lindemann++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_falloff_troe_reactions[count_troe++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_FALLOFF ) 
				{
					indices_of_falloff_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_falloff_sri_reactions[count_sri++] = it - reactions_.begin() + 1;
				}
			}

			fOutput << "falloff-reactions" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << indices_of_falloff_reactions[i] << " ";
			fOutput << std::endl;
			
			fOutput << "falloff-lindemann-reactions" << std::endl;
			fOutput << number_of_falloff_lindemann_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_lindemann_reactions;i++)
				fOutput << indices_of_falloff_lindemann_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "falloff-troe-reactions" << std::endl;
			fOutput << number_of_falloff_troe_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_troe_reactions;i++)
				fOutput << indices_of_falloff_troe_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "falloff-sri-reactions" << std::endl;
			fOutput << number_of_falloff_sri_reactions << std::endl;
			for(unsigned int i=0;i<number_of_falloff_sri_reactions;i++)
				fOutput << indices_of_falloff_sri_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			unsigned int count_lindemann=0;
			unsigned int count_troe=0;
			unsigned int count_sri=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				if ( (*it).Tag() == PhysicalConstants::REACTION_LINDEMANN_CABR ) 
				{
					indices_of_cabr_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_cabr_lindemann_reactions[count_lindemann++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_TROE_CABR ) 
				{
					indices_of_cabr_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_cabr_troe_reactions[count_troe++] = it - reactions_.begin() + 1;
				}
				else if ( (*it).Tag() == PhysicalConstants::REACTION_SRI_CABR ) 
				{
					indices_of_cabr_reactions[count++] = it - reactions_.begin() + 1;
					indices_of_cabr_sri_reactions[count_sri++] = it - reactions_.begin() + 1;
				}
			}

			fOutput << "cabr-reactions" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << indices_of_cabr_reactions[i] << " ";
			fOutput << std::endl;
			
			fOutput << "cabr-lindemann-reactions" << std::endl;
			fOutput << number_of_cabr_lindemann_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_lindemann_reactions;i++)
				fOutput << indices_of_cabr_lindemann_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "cabr-troe-reactions" << std::endl;
			fOutput << number_of_cabr_troe_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_troe_reactions;i++)
				fOutput << indices_of_cabr_troe_reactions[i] << " ";
			fOutput << std::endl;

			fOutput << "cabr-sri-reactions" << std::endl;
			fOutput << number_of_cabr_sri_reactions << std::endl;
			for(unsigned int i=0;i<number_of_cabr_sri_reactions;i++)
				fOutput << indices_of_cabr_sri_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_CHEBYSHEV ) indices_of_chebyshev_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "chebyshev-reactions" << std::endl;
			fOutput << number_of_chebyshev_reactions << std::endl;
			for(unsigned int i=0;i<number_of_chebyshev_reactions;i++)
				fOutput << indices_of_chebyshev_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsPressureLog() == true) indices_of_pressurelog_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "pressurelog-reactions" << std::endl;
			fOutput << number_of_pressurelog_reactions << std::endl;
			for(unsigned int i=0;i<number_of_pressurelog_reactions;i++)
				fOutput << indices_of_pressurelog_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsFit1() == true) indices_of_fit1_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "fit1-reactions" << std::endl;
			fOutput << number_of_fit1_reactions << std::endl;
			for(unsigned int i=0;i<number_of_fit1_reactions;i++)
				fOutput << indices_of_fit1_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsJanevLanger() == true) indices_of_janevlanger_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "janevlanger-reactions" << std::endl;
			fOutput << number_of_janevlanger_reactions << std::endl;
			for(unsigned int i=0;i<number_of_janevlanger_reactions;i++)
				fOutput << indices_of_janevlanger_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			unsigned int count=0;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
				if ( (*it).Tag() == PhysicalConstants::REACTION_SIMPLE && (*it).IsLandauTeller() == true) indices_of_landauteller_reactions[count++] = it - reactions_.begin() + 1;

			fOutput << "landauteller-reactions" << std::endl;
			fOutput << number_of_landauteller_reactions << std::endl;
			for(unsigned int i=0;i<number_of_landauteller_reactions;i++)
				fOutput << indices_of_landauteller_reactions[i] << " ";
			fOutput << std::endl;
		}

		{		
			OpenSMOKEVectorDouble lnA(reactions_.size());
			OpenSMOKEVectorDouble Beta(reactions_.size());
			OpenSMOKEVectorDouble E_over_R(reactions_.size());
			
			unsigned int j=1;
			for (std::vector<ReactionPolicy_CHEMKIN>::const_iterator it = reactions_.begin(); it != reactions_.end(); ++it)
			{
				//lnA[j] = log((*it).A());
				lnA[j] = ((*it).A() == 0.) ? log(1.e-100) : log((*it).A());
				
				Beta[j] = (*it).Beta();
				E_over_R[j] = (*it).E_over_R();
				j++;
			}

			fOutput << "lnA" << std::endl;
			lnA.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "Beta" << std::endl;
			Beta.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "E_over_R" << std::endl;
			E_over_R.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
		}

		{		
			OpenSMOKEVectorDouble lnA_reversible(number_of_explicitly_reversible_reactions);
			OpenSMOKEVectorDouble Beta_reversible(number_of_explicitly_reversible_reactions);
			OpenSMOKEVectorDouble E_over_R_reversible(number_of_explicitly_reversible_reactions);

			unsigned int j=1;
			for (unsigned int i=0;i<number_of_explicitly_reversible_reactions;i++)
			{
				double A = reactions_[indices_of_explicitly_reversible_reactions[i]-1].A_reversible();
				lnA_reversible[j] = (A == 0.) ? log(1.e-100) : log(A);
				Beta_reversible[j] = reactions_[indices_of_explicitly_reversible_reactions[i]-1].Beta_reversible();
				E_over_R_reversible[j] = reactions_[indices_of_explicitly_reversible_reactions[i]-1].E_over_R_reversible();
				j++;
			}

			fOutput << "reversible-lnA" << std::endl;
			lnA_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "reversible-Beta" << std::endl;
			Beta_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;

			fOutput << "reversible-E_over_R" << std::endl;
			E_over_R_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
			fOutput << std::endl;
		}

		{		
			fOutput << "thirdbody-parameters" << std::endl;
			for (unsigned int i=0;i<number_of_thirdbody_reactions;i++)
				reactions_[indices_of_thirdbody_reactions[i]-1].WriteThirdBodyParametersOnASCIIFile(fOutput);
		}

		{		
			fOutput << "falloff-kinetics" << std::endl;

			fOutput << "lnA-falloff-inf" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << log( reactions_[indices_of_falloff_reactions[i]-1].A_inf() )<< " ";
			fOutput << std::endl;

			fOutput << "Beta-falloff-inf" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << reactions_[indices_of_falloff_reactions[i]-1].Beta_inf() << " ";
			fOutput << std::endl;

			fOutput << "E_over_R-falloff-inf" << std::endl;
			fOutput << number_of_falloff_reactions << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				fOutput << reactions_[indices_of_falloff_reactions[i]-1].E_over_R_inf() << " ";
			fOutput << std::endl;
		}

		{		
			fOutput << "falloff-parameters" << std::endl;
			for (unsigned int i=0;i<number_of_falloff_reactions;i++)
				reactions_[indices_of_falloff_reactions[i]-1].WritePressureDependentParametersOnASCIIFile(fOutput);
		}

		{		
			fOutput << "cabr-kinetics" << std::endl;
			fOutput << "lnA-cabr-inf" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << log( reactions_[indices_of_cabr_reactions[i]-1].A_inf() )<< " ";
			fOutput << std::endl;

			fOutput << "Beta-cabr-inf" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << reactions_[indices_of_cabr_reactions[i]-1].Beta_inf()<< " ";
			fOutput << std::endl;

			fOutput << "E_over_R-cabr-inf" << std::endl;
			fOutput << number_of_cabr_reactions << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				fOutput << reactions_[indices_of_cabr_reactions[i]-1].E_over_R_inf() << " ";
			fOutput << std::endl;
		}

		{		
			fOutput << "cabr-parameters" << std::endl;
			for (unsigned int i=0;i<number_of_cabr_reactions;i++)
				reactions_[indices_of_cabr_reactions[i]-1].WritePressureDependentParametersOnASCIIFile(fOutput);
		}

		{
			fOutput << "chebyshev-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_chebyshev_reactions;j++)
				reactions_[indices_of_chebyshev_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "pressurelog-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_pressurelog_reactions;j++)
				reactions_[indices_of_pressurelog_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "fit1-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_fit1_reactions;j++)
				reactions_[indices_of_fit1_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "janevlanger-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_janevlanger_reactions;j++)
				reactions_[indices_of_janevlanger_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}
		{
			fOutput << "landauteller-parameters " << std::endl;
			for(unsigned int j=0;j<number_of_landauteller_reactions;j++)
				reactions_[indices_of_landauteller_reactions[j]-1].WriteAdditionalDataOnASCIIFile(fOutput);
		}

		{		
			fOutput << "stoichiometric-coefficients" << std::endl;
			WriteStoichiometricDataOnASCIIFile(fOutput);
			WriteReactionOrdersOnASCIIFile(fOutput);
		}

		return true;

	}
	*/

/*	template<typename T>
	void FormatXML(std::stringstream& xml_string, const std::string tag, std::vector<T>& v, bool close=true)
	{
			xml_string<< "<" << tag << ">" << std::endl;
			xml_string << v.size() << std::endl;
			if (v.size()>0)
			{
				for(unsigned int i=0;i<v.size();i++)
				{
					xml_string << v[i] << " ";
					if ((i+1)%30==0) xml_string << std::endl;
				}
				xml_string << std::endl;
			}
			if (close == true) xml_string << "</" << tag << ">" << std::endl;
	}
*/

	template<typename Reactions>
	bool PreProcessorSolidKineticsPolicy_CHEMKIN<Reactions>::WriteXMLFile(std::stringstream& fOutput) const
	{
		std::cout << " * Writing the interpreted kinetic file in XML format..." << std::endl;

		// Number of materials
		fOutput << "<NumberOfMaterials>" << std::endl;
		fOutput << number_of_materials_ << std::endl;
		fOutput << "</NumberOfMaterials>" << std::endl;

		// Number of solid species
		fOutput << "<NumberOfSolidSpecies>" << std::endl;
		fOutput << names_solid_species_.size() << std::endl;
		fOutput << "</NumberOfSolidSpecies>" << std::endl;

		// List of site species
		fOutput << "<SolidSpecies>" << std::endl;
		for(unsigned int i=0;i<names_solid_species_.size();i++)
			fOutput << names_solid_species_[i] << std::endl;
		fOutput << "</SolidSpecies>" << std::endl;

		for(unsigned int k=0;k<number_of_materials_;k++)
		{
			fOutput << "<MaterialDescription index=\"" << k+1 << "\" name=\"" << names_of_materials_[k] << "\">" << std::endl;

			fOutput << solid_species_[k].size() << std::endl;
			for(unsigned int i=0;i<solid_species_[k].size();i++)
				fOutput << solid_species_[k][i] << " " << solid_species_indices_[k][i] << std::endl;

			fOutput << "</MaterialDescription>" << std::endl;
		}

		
		fOutput << "<Kinetics type=\"OpenSMOKE\" version=\"01-02-2014\">" << std::endl;
		for(unsigned int k=0;k<number_of_materials_;k++)
		{
			fOutput << "<MaterialKinetics index=\"" << k+1 << "\" name=\"" << names_of_materials_[k] << "\">" << std::endl;

			// 1. Number of species
			fOutput << "<NumberOfSpecies>" << std::endl;
			fOutput << names_species_.size() << std::endl;
			fOutput << "</NumberOfSpecies>" << std::endl;

			// 2. Number of reactions
			fOutput << "<NumberOfReactions>" << std::endl;
			fOutput << reactions_[k].size() << std::endl;
			fOutput << "</NumberOfReactions>" << std::endl;
			
			// 3. Kinetic data for each reaction
			unsigned int number_of_irreversible_reactions = 0;
			unsigned int number_of_reversible_reactions = 0;
			unsigned int number_of_explicitly_reversible_reactions = 0;
			unsigned int number_of_thermodynamic_reversible_reactions = 0;

		
			for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[k].begin(); it != reactions_[k].end(); ++it)
			{
				if ( (*it).IsReversible()==true)
				{
					number_of_reversible_reactions++;
					if ( (*it).IsExplicitlyReversible()==true)	number_of_explicitly_reversible_reactions++;
					else										number_of_thermodynamic_reversible_reactions++;
				}
				else
					number_of_irreversible_reactions++;
			
				if ( (*it).Tag() == PhysicalConstants::REACTION_SURFACE_SIMPLE)
				{
				}
			}

			std::vector<unsigned int> indices_of_irreversible_reactions(number_of_irreversible_reactions);
			std::vector<unsigned int> indices_of_reversible_reactions(number_of_reversible_reactions);
			std::vector<unsigned int> indices_of_explicitly_reversible_reactions(number_of_explicitly_reversible_reactions);
			std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions(number_of_thermodynamic_reversible_reactions);

			{		
				unsigned int count_irreversible=0;
				unsigned int count_reversible=0;
				unsigned int count_explicit=0;
				unsigned int count_thermodynamic=0;
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[k].begin(); it != reactions_[k].end(); ++it)
				{
					if ( (*it).IsReversible()==true) 
					{
							indices_of_reversible_reactions[count_reversible++] = boost::lexical_cast<int>(it - reactions_[k].begin() + 1);
							if ((*it).IsExplicitlyReversible() == true)  indices_of_explicitly_reversible_reactions[count_explicit++] = boost::lexical_cast<int>(it - reactions_[k].begin() + 1);
							else										  indices_of_thermodynamic_reversible_reactions[count_thermodynamic++] = boost::lexical_cast<int>(it - reactions_[k].begin() + 1);
					}
					else
						indices_of_irreversible_reactions[count_irreversible++] = boost::lexical_cast<int>(it - reactions_[k].begin() + 1);
				}

				FormatXML(fOutput, "Irreversible", indices_of_irreversible_reactions);
				FormatXML(fOutput, "Reversible", indices_of_reversible_reactions);
				FormatXML(fOutput, "Reversible-Thermodynamics", indices_of_thermodynamic_reversible_reactions);
				FormatXML(fOutput, "Reversible-Explicit", indices_of_explicitly_reversible_reactions);

			}	

			// Reactions Needing Conversion
			{
				unsigned int count=0;
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[k].begin(); it != reactions_[k].end(); ++it)
					if ((*it).composition_units() != PhysicalConstants::UNITS_STD)
						count++;

				fOutput << "<ReactionsNeedingConversion>" << std::endl;
				fOutput << count << std::endl;
				unsigned int j=0;
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[k].begin(); it != reactions_[k].end(); ++it)
				{
					j++;
					if ((*it).composition_units() != PhysicalConstants::UNITS_STD)
						fOutput << j << std::endl;
				}
				fOutput << "</ReactionsNeedingConversion>" << std::endl;							
			}

			// Type of kinetics
			{
				fOutput << "<TypeOfKinetics>" << std::endl;
				fOutput << "chemkin_conventional" << std::endl;
				fOutput << "</TypeOfKinetics>" << std::endl;
			}

			// Thermodynamic reversible reactions
			{
				fOutput << "<ThermodynamicReversibleReactions>" << std::endl;
				fOutput << number_of_thermodynamic_reversible_reactions << std::endl;

				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[k].begin(); it != reactions_[k].end(); ++it)
					if ( (*it).IsReversible()==true && (*it).IsExplicitlyReversible() == false)
						fOutput << (*it).delta_nu_gas() << std::endl;
				
				fOutput << "</ThermodynamicReversibleReactions>" << std::endl;							
			}

			// Conventional parameters
			{		
				OpenSMOKEVectorDouble lnA(boost::lexical_cast<int>(reactions_[k].size()));
				OpenSMOKEVectorDouble Beta(boost::lexical_cast<int>(reactions_[k].size()));
				OpenSMOKEVectorDouble E_over_R(boost::lexical_cast<int>(reactions_[k].size()));
				OpenSMOKEVectorDouble forward_kinetic_order(boost::lexical_cast<int>(reactions_[k].size()));
			
				unsigned int j=1;
				for (std::vector<ReactionPolicy_Solid_CHEMKIN>::const_iterator it = reactions_[k].begin(); it != reactions_[k].end(); ++it)
				{
					//lnA[j] = log((*it).A());
					lnA[j] = ((*it).A() == 0.) ? std::log(1.e-100) : std::log((*it).A());
					Beta[j] = (*it).Beta();
					E_over_R[j] = (*it).E_over_R();
					forward_kinetic_order[j] = (*it).forward_gas_kinetic_order();
					j++;
				}

				fOutput << "<KineticParameters>" << std::endl;

				fOutput << "<Direct>" << std::endl;

				fOutput << "<lnA>" << std::endl;
				fOutput.precision(6);
				lnA.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</lnA>" << std::endl;

				fOutput << "<Beta>" << std::endl;
				Beta.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</Beta>" << std::endl;

				fOutput << "<E_over_R>" << std::endl;
				E_over_R.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</E_over_R>" << std::endl;

				fOutput << "<ForwardKineticOrder>" << std::endl;
				forward_kinetic_order.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</ForwardKineticOrder>" << std::endl;

				fOutput << "</Direct>" << std::endl;	
			}

			if (number_of_explicitly_reversible_reactions != 0)
			{		
				OpenSMOKEVectorDouble lnA_reversible(number_of_explicitly_reversible_reactions);
				OpenSMOKEVectorDouble Beta_reversible(number_of_explicitly_reversible_reactions);
				OpenSMOKEVectorDouble E_over_R_reversible(number_of_explicitly_reversible_reactions);

				unsigned int j=1;
				for (unsigned int i=0;i<number_of_explicitly_reversible_reactions;i++)
				{
					double A = reactions_[k][indices_of_explicitly_reversible_reactions[i]-1].A_reversible();
					lnA_reversible[j] = (A == 0.) ? std::log(1.e-100) : std::log(A);
					Beta_reversible[j] = reactions_[k][indices_of_explicitly_reversible_reactions[i]-1].Beta_reversible();
					E_over_R_reversible[j] = reactions_[k][indices_of_explicitly_reversible_reactions[i]-1].E_over_R_reversible();
					j++;
				}

				fOutput << "<Reverse>" << std::endl;

				fOutput << "<lnA>" << std::endl;
				lnA_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</lnA>" << std::endl;

				fOutput << "<Beta>" << std::endl;
				Beta_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</Beta>" << std::endl;

				fOutput << "<E_over_R>" << std::endl;
				E_over_R_reversible.Save(fOutput, OPENSMOKE_FORMATTED_FILE);
				fOutput << std::endl;
				fOutput << "</E_over_R>" << std::endl;

				fOutput << "</Reverse>" << std::endl;
			}

			fOutput << "</KineticParameters>" << std::endl;
		
			// Write stoichiometric data
			{		
				fOutput << "<Stoichiometry type=\"OpenSMOKE\" version=\"01-02-2014\">" << std::endl;
				WriteStoichiometricDataOnASCIIFile(fOutput,k);
				WriteReactionOrdersOnASCIIFile(fOutput,k);
				fOutput << false << std::endl; // TOREMOVE
				fOutput << "</Stoichiometry>" << std::endl;
			}

			// Global reactions (TODO)
			{
				fOutput << "<GlobalReactions type=\"none\">" << std::endl;
				fOutput << "</GlobalReactions>" << std::endl;
			}
			
			fOutput << "</MaterialKinetics>" << std::endl;
		
		}	
		
		fOutput << "</Kinetics>" << std::endl;
		
		return true;

	}
}
