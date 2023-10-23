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
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of catalyticFOAM framework.                         |
|                                                                         |
|	License                                                           |
|                                                                         |
|   Copyright(C) 2014-2011, A.Cuoci, M.Maestri,                           |
|                2014-2013, S.Rebughini                                   |
|                     2013, T.Maffei                                      |
|                     2013, G.Gentile, F.Manelli                          |
|                     2012, M.Calonaci, F.Furnari                         |
|                     2011, S.Goisis, A.Osio                              |
|                                                                         |
|   catalyticFOAM is free software: you can redistribute it and/or modify |
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
|   along with catalyticFOAM. If not, see <http://www.gnu.org/licenses/>. |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "BatchReactorHeterogeneousConstantPressure.H"

BatchReactorHeterogeneousConstantPressure::BatchReactorHeterogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
																						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
																						OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, 
																						OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap),
	thermodynamicsSurfaceMap_(thermodynamicsSurfaceMap), 
	kineticsSurfaceMap_(kineticsSurfaceMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		SURF_NP_ = thermodynamicsSurfaceMap_.number_of_site_phases(0);
		SURF_NC_ = thermodynamicsSurfaceMap_.number_of_site_species();
		NE_ = NC_ + SURF_NC_ + SURF_NP_ + 1 + 1;
		
		QRGas_     = 0.;
		QRSurface_ = 0.;
				
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &RfromGas_, true);
		ChangeDimensions(NC_, &RfromSurface_, true);
		ChangeDimensions(NC_, &MWs_, true);
		
		ChangeDimensions(SURF_NC_, &Z_, true);
		ChangeDimensions(SURF_NP_, &Gamma_, true);
		ChangeDimensions(SURF_NC_, &Rsurface_, true);
		ChangeDimensions(SURF_NP_, &RsurfacePhases_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;
		reactionHeatFromHeterogeneousReactions_ = true;
		homogeneousReactions_ = false;
		
		MWs_.CopyFrom(&thermodynamicsMap_.MWs()[0]);
	}

void BatchReactorHeterogeneousConstantPressure::SetReactor( const double P0, const double A, const double alfaCatalyst)
{
	P0_    = P0;
	A_     = A;
	alfaCatalyst_ = alfaCatalyst;
}

int BatchReactorHeterogeneousConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	// Recover unknowns
	{
		unsigned int k=1;
		if (checkMassFractions_ == true)
		{
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = std::max(y[k++], 0.);
			mass_ = std::max(y[k++], 0.);
			for(unsigned int i=1;i<=SURF_NP_;++i)
				Gamma_[i] = std::max(y[k++], 0.);
			for(unsigned int i=1;i<=SURF_NC_;++i)
				Z_[i] = std::max(y[k++], 0.);
		}
		else
		{
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = y[k++];
			mass_ = y[k++];
			for(unsigned int i=1;i<=SURF_NP_;++i)
				Gamma_[i] = y[k++];
			for(unsigned int i=1;i<=SURF_NC_;++i)
				Z_[i] = y[k++];
		}
		
		T_ = y[k++];
	}
	
	// Calculates the volume and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
	const double cTot_ = P0_/(PhysicalConstants::R_J_kmol * T_);
	Product(cTot_, x_, &c_);
	const double rho_ = cTot_*MW_;
	const double V_ = mass_/rho_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);
		
	// Calculates homogeneous kinetics
	{
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P0_);

		if ( homogeneousReactions_ == true)
		{
			kineticsMap_.ReactionRates(c_.GetHandle());
			kineticsMap_.FormationRates(RfromGas_.GetHandle());
		}
	}

	// Calculates heterogeneous kinetics
    thermodynamicsSurfaceMap_.SetPressure(P0_);
	thermodynamicsSurfaceMap_.SetTemperature(T_);
	kineticsSurfaceMap_.SetPressure(P0_);
	kineticsSurfaceMap_.SetTemperature(T_);
    
	kineticsSurfaceMap_.ReactionEnthalpiesAndEntropies();
	kineticsSurfaceMap_.KineticConstants();
	kineticsSurfaceMap_.ReactionRates(c_.GetHandle(), Z_.GetHandle(), dummy.GetHandle(), Gamma_.GetHandle());
	kineticsSurfaceMap_.FormationRates(RfromSurface_.GetHandle(), Rsurface_.GetHandle(), dummy.GetHandle(), RsurfacePhases_.GetHandle());
	
	double dm_over_dt = A_*alfaCatalyst_*Dot(RfromSurface_, MWs_);

	// Recovering residuals
	{
		unsigned int k=1;
		
		// Gas phase species
		for (unsigned int i=1;i<=NC_;++i)	
		{
			dy[k++] = thermodynamicsMap_.MW(i-1)*RfromGas_[i]/rho_ +
					 (-omega_[i]*dm_over_dt + A_*alfaCatalyst_*RfromSurface_[i]*thermodynamicsMap_.MW(i-1))/mass_;
		
		}
		
		// Total mass
		dy[k++] = dm_over_dt;
		
		// Phases
		for (unsigned int i=1;i<=SURF_NP_;++i)	
			dy[k++] = RsurfacePhases_[i];
		
		// Surface site species
		for (unsigned int i=1;i<=SURF_NC_;++i)	
		{
			const unsigned int index_phase = thermodynamicsSurfaceMap_.vector_site_phases_belonging()[i-1]+1;
			dy[k++] = (thermodynamicsSurfaceMap_.vector_occupancies_site_species()[i-1]*Rsurface_[i] - 
									Z_[i]*dy[NC_+1+index_phase])/Gamma_[index_phase];
                                   
		}
		
		// Energy equation
		if (energyEquation_ == true)
		{
			double CpMixMolar; 
			CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
			const double CpMixMass_ = CpMixMolar / MW_;
			
			if ( homogeneousReactions_ == true)
				QRGas_ = kineticsMap_.HeatRelease(RfromGas_.GetHandle());

			if (reactionHeatFromHeterogeneousReactions_ == true)
				QRSurface_ = kineticsSurfaceMap_.HeatRelease(RfromSurface_.GetHandle(), Rsurface_.GetHandle(), dummy.GetHandle());
			
			dy[k++] = (V_*QRGas_+ A_*alfaCatalyst_*QRSurface_) / (mass_*CpMixMass_);
		}
		else
		{
			dy[k++] = 0.;
		}
	}
	
 	//std::cout << "\nVelocità di Reazione H2 = " << RfromSurface_[1] << std::endl;
	
	return 0;
}

int BatchReactorHeterogeneousConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	return 0;
}
