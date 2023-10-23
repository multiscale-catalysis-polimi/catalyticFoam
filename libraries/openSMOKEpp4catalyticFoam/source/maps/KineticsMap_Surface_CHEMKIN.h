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

#ifndef OpenSMOKE_KineticsMap_Surface_CHEMKIN_CHEMKIN_H
#define OpenSMOKE_KineticsMap_Surface_CHEMKIN_CHEMKIN_H

#include "KineticsMap.h"
#include "StoichiometricMap.h"

namespace OpenSMOKE
{
	class UBIQEP_SubMechanism
	{
	public:

		void ReadFromXMLFile(boost::property_tree::ptree& ptree, const std::vector<double>& MWs);
		void CalculateChemisorptionHeats(const double T, const double* Z);
		void CalculateDissociationEnergies(const std::vector<double>& H);
		void CalculateSurfaceEnthalpies(const double T);
		void CalculateActivationEnergies();
		void ForwardKineticConstants(const double T, const double total_site_density, double* kForward);

	private:
		unsigned int number_of_gas_species_;
		unsigned int number_of_site_species_;
		unsigned int half_number_of_ubiqep_reactions_;
		double reference_temperature_;

		std::vector<unsigned int> chemisorption_heats_gas_indices__;
		std::vector<double> chemisorption_heats_constant_coefficient__;
		std::vector<double> chemisorption_heats_temperature_coefficient__;
		std::vector<double>* chemisorption_heats_coefficients__;
		std::vector<unsigned int>* chemisorption_heats_indices__;
		std::vector<double> QStar__;
		std::vector<double> dissociation_energies__;
		std::vector<double> surface_enthalpies__;
		std::vector<double> E_forward__;
		std::vector<double> E_backward__;
		std::vector<double> adsorption_coefficient__;


		std::vector<double> sigma_direct__;
		std::vector<double> Beta_direct__;
		std::vector<double> lnA_direct__;

		std::vector<double> sigma_reverse__;
		std::vector<double> Beta_reverse__;
		std::vector<double> lnA_reverse__;

		std::vector<double> lambda__;
			
		std::vector<unsigned int> ubiqep_class__;
		std::vector<PhysicalConstants::UBIQEP_TYPE> ubiqep_type_direct__;
		std::vector<PhysicalConstants::UBIQEP_TYPE> ubiqep_type_reverse__;

		std::vector<unsigned int> index_A__;
		std::vector<unsigned int> index_B__;
		std::vector<unsigned int> index_C__;
		std::vector<unsigned int> index_D__;
		std::vector<unsigned int> index_Star__;
		std::vector<unsigned int> index_A2__;
		std::vector<unsigned int> index_AB__;
		std::vector<unsigned int> index_ABStar__;
		std::vector<unsigned int> index_AStar__;
		std::vector<unsigned int> index_BStar__;
		std::vector<unsigned int> index_CStar__;
		std::vector<unsigned int> index_DStar__;

		std::vector<unsigned int> index_A2_Chemisorption_Heats__;
		std::vector<unsigned int> index_AB_Chemisorption_Heats__;
	};

	enum TYPE_OF_KINETICS { TYPE_OF_KINETICS_CHEMKIN_CONVENTIONAL, TYPE_OF_KINETICS_UBI_QEP };

	//!  A class to efficiently evaluate the reaction and formation rates, to be used in production codes
	/*!
		 This class provides the tools to calculate in a very efficient way the reaction rates and the
		 formation rates. In order to ensure a good efficiency a map is created to store all the data
		 depending on the temperature. Inthis way they are recalculated only if strictly needed, i.e. only
		 if the temperature changes
	*/

	class KineticsMap_Surface_CHEMKIN : public KineticsMap
	{

	public:

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties (obsolete, TOREMOVE)
		*@param thermo the thermodynamic map
		*@param nSpecies number of species 
		*/
		KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN& thermo, const unsigned int nSpecies);

		/**
		*@brief Creates a thermodynamic map for the evaluation of thermodynamic properties
		*@param thermo the thermodynamic map
		*@param doc xml file  
		*/
		KineticsMap_Surface_CHEMKIN(ThermodynamicsMap_Surface_CHEMKIN& thermo, boost::property_tree::ptree& ptree);

		/**
		*@brief Set the temperature at which the properties have to be evaluated
		*@param T the temperature value in K
		*/
		virtual void SetTemperature(const double& T);

		/**
		*@brief Set the pressure at which the properties have to be evaluated
		*@param P the pressure value in Pa
		*/
		virtual void SetPressure(const double& P);		

		/**
		*@brief Returns the names of the species
		*/
		const std::vector<std::string>& NamesOfSpecies() const { return thermodynamics_.names(); }

		/**
		*@brief Imports the kinetic schemes from a file in XML format
		*/
		virtual void ImportCoefficientsFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Imports the list of species from a file in XML format
		*/
		virtual void ImportSpeciesFromXMLFile(boost::property_tree::ptree& ptree);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const double* x_bath, const unsigned int nparameters, Eigen::MatrixXd& fittedKineticParameters);

		/**
		*@brief Calculates the kinetic constants of the reverse reactions
		*/
		void FittedReverseKineticConstants(const unsigned int k, std::ostream& fOut, Eigen::MatrixXd& fittedKineticParameters);

		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k, const double* c_bath, const double conversion_forward=1., const double conversion_backward=1.);
		
		/**
		*@brief Write the data for the reaction tables
		*/
		void WriteKineticData(std::ostream& fOut, const unsigned int k);

		/**
		*@brief Calculates the formation rates for all the species in the kinetic mechanism
		*/
		void FormationRates(double* Rgas, double* Rsite, double* Rbulk, double* RsitePhases);
		void FormationRates(double* R);

		/**
		*@brief Calculates the heat release
		*/
		double HeatRelease(const double* Rgas, const double* Rsurface, const double* Rbulk);

		/**
		*@brief Calculates the production and the destruction rates for all the species in the kinetic mechanism
		*/
		void ProductionAndDestructionRates(double* P, double* D);

		void RateOfProductionAnalysis(ROPA_Data& ropa) const;
		void RateOfProductionAnalysis(ROPA_Data& ropa, const double* rf, const double* rb) const;


		/**
		*@brief Returns the forward reaction rates for all the reactions in the kinetic scheme
		*/
		void GetForwardReactionRates(double* r);

		/**
		*@brief Returns the backward reaction rates for all the reactions in the kinetic scheme
		        If a reaction is irreversible, it returns zero
		*/
		void GetBackwardReactionRates(double* r);

		/**
		*@brief Calculates the reaction rates for all the reactions in the kinetic scheme
		*/
		void ReactionRates(const double* c, const double* z, const double* a, const double* gamma);

		/**
		*@brief Calculates the reaction rates for all the lumped reactions in the kinetic scheme
		*/
		virtual void UserDefinedReactionRates(const double* c, const double* z, const double* a, const double* gamma);

		/**
		*@brief Returns the indices of the reversible reactions
		*/
		const std::vector<unsigned int>& IndicesOfReversibleReactions() const { return indices_of_reversible_reactions__; }

		/**
		*@brief Calculates the reaction enthalpies and entropies (to be used for the kinetic constants)
		*/
		void ReactionEnthalpiesAndEntropies();

		/**
		*@brief Calculates the kinetic constants
		*/
		void KineticConstants();

		/**
		*@brief Return the net reaction rates in [kmol/m2/s]
		*/
		const std::vector<double>& GiveMeReactionRates();

		/**
		*@brief Return the net reaction rates in [kmol/m2/s]
		*/
		void GiveMeReactionRates(double* r);

		const std::vector<double>& KArrheniusModified() const { return kArrheniusModified__; }
		const std::vector<double>& KArrhenius() const { return kArrhenius__; }

		StoichiometricMap& stoichiometry() { return *stoichiometry_; }

	protected:

		ThermodynamicsMap_Surface_CHEMKIN& thermodynamics_;		//!< reference to the thermodynamics

		std::vector<double> cSites__;
		std::vector<double> c__;
		std::vector<double> aux_vector__;

		std::vector<unsigned int> indices_of_irreversible_reactions__;				//!< indices of irreversible reactions
		std::vector<unsigned int> indices_of_reversible_reactions__;				//!< indices of reversible reactions
		std::vector<unsigned int> indices_of_thermodynamic_reversible_reactions__;	//!< indices of reversible (thermodynamic) reactions
		std::vector<unsigned int> indices_of_explicitly_reversible_reactions__;		//!< indices of reversible (explicit) reactions
		std::vector<unsigned int> indices_of_stick_reactions__;						//!< indices of stick reactions
		std::vector<unsigned int> indices_of_coverage_dependent_reactions__;		//!< indices of coverage dependent reactions
		std::vector<unsigned int> indices_of_langmuir_reactions__;					//!< indices of coverage dependent reactions
		std::vector<unsigned int> indices_of_lumped_reactions__;					//!< indices of lumped reactions

		unsigned int number_of_irreversible_reactions_;
		unsigned int number_of_reversible_reactions_;
		unsigned int number_of_thermodynamic_reversible_reactions_;
		unsigned int number_of_explicitly_reversible_reactions_;
		unsigned int number_of_stick_reactions_;
		unsigned int number_of_coverage_dependent_reactions_;
		unsigned int number_of_langmuir_reactions_;
		unsigned int number_of_lumped_reactions_;

		std::vector<double> lnA__;							//!< frequency factors (log)
		std::vector<double> Beta__;							//!< temperature exponents
		std::vector<double> E_over_R__;						//!< activation temperatures
		std::vector<double> forward_kinetic_order__;		//!< global kinetic order for forward reactions

		std::vector<unsigned int> indices_of_reactions_needing_conversion_;

		std::vector<double> lnA_reversible__;			//!< frequency factors (log) for explicitly reversible reactions
		std::vector<double> Beta_reversible__;			//!< temperature exponents for explicitly reversible reactions
		std::vector<double> E_over_R_reversible__;		//!< activation temperatures for explicitly reversible reactions

		std::vector<double> changeOfMoles__;			//!< list of change of moles

		StoichiometricMap* stoichiometry_;			//!< pointer to the stoichiometry

		bool arrhenius_kinetic_constants_must_be_recalculated_;
		bool nonconventional_kinetic_constants_must_be_recalculated_;
		bool reaction_h_and_s_must_be_recalculated_;

		std::vector<double> reaction_s_over_R__;
		std::vector<double> reaction_h_over_RT__;
		std::vector<double> kArrheniusModified__;
		std::vector<double> kArrhenius__;
		std::vector<double> kArrhenius_reversible__;
		std::vector<double> uKeq__;

		std::vector<double> forwardReactionRates__;
		std::vector<double> reverseReactionRates__;
		std::vector<double> netReactionRates__;

		double Patm_over_RT_;
		double log_Patm_over_RT_;

		std::vector<unsigned int> isThermodynamicallyReversible__;		//!< vector containing the local index of thermodynamically reversible reactions
		std::vector<unsigned int> isExplicitlyReversible__;				//!< vector containing the local index of explicitly reversible reactions

		VectorReactionTags type_of_reaction__;
		std::vector<unsigned int> local_family_index__;

		// Stick reactions
		std::vector<double> stick_constant_coefficient_;
		std::vector<double> stick_power_;
		std::vector<bool> stick_motz_wise_;

		// Coverage dependent reactions
		std::vector< std::vector<bool> >			coverage_dependent_species_site_type_;
		std::vector< std::vector<unsigned int> >	coverage_dependent_species_index_;
		std::vector< std::vector<double> >			coverage_dependent_eta_;
		std::vector< std::vector<double> >			coverage_dependent_mu_;
		std::vector< std::vector<double> >			coverage_dependent_epsilon_;

		// Langmuir-Hinshelwood reactions
		std::vector< double >											langmuir_denominator_order_;
		std::vector< PhysicalConstants::UNITS_REACTION_COMPOSITION >	langmuir_units_;
		std::vector< std::vector<unsigned int> >						langmuir_species_index_;
		std::vector< std::vector<bool> >								langmuir_numerator_species_;
		std::vector< std::vector<double> >								langmuir_lnA_;
		std::vector< std::vector<double> >								langmuir_Beta_;
		std::vector< std::vector<double> >								langmuir_H_over_R_;
		std::vector< std::vector<double> >								langmuir_order_;

		// Lumped reactions
		std::vector< std::string >										names_of_lumped_functions_;

		// Non conservation of sites
		std::vector<unsigned int>	non_conservation_of_sites_indices_of_reactions_;
		std::vector<unsigned int>	non_conservation_of_sites_phase_of_reactions_;
		std::vector<double>			non_conservation_of_sites_delta_sigma_;

		// Thermodynamic reversible reactions
		std::vector<double> delta_sigma_times_log_Gamma0_;
		std::vector<double> delta_nu_gas_;

		// UBI
		TYPE_OF_KINETICS type_of_kinetics_;
		UBIQEP_SubMechanism* ubiqep_submechanism_;
	};


}

#include "KineticsMap_Surface_CHEMKIN.hpp"

#endif /* OpenSMOKE_KineticsMap_Surface_CHEMKIN_H */
