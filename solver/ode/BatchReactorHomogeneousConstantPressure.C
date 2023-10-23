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

#include "BatchReactorHomogeneousConstantPressure.H"

BatchReactorHomogeneousConstantPressure::BatchReactorHomogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NE_ = NC_+1;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;
	}

void BatchReactorHomogeneousConstantPressure::SetReactor( const double P0 )
{
	P0_    = P0;
}

int BatchReactorHomogeneousConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	// Recover mass fractions
	if (checkMassFractions_ == true)
	{	for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = std::max(y[i], 0.);
	}
	else
	{
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = y[i];
	}

	// Recover temperature
	T_ = y[NC_+1];

	// Calculates the pressure and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
	cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
    rho_ = cTot_*MW_;
	Product(cTot_, x_, &c_);

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);
	
	// Calculates kinetics
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P0_);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_.GetHandle());
	kineticsMap_.FormationRates(R_.GetHandle());
	

	// Species equations
	for (unsigned int i=1;i<=NC_;++i)	
		dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho_;
		
	// Energy equation
   	dy[NC_+1] = 0.;     
    	if (energyEquation_ == true)
   	{
		double CpMixMolar; 
		CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_.GetHandle());
		dy[NC_+1]  = QR_ / (rho_*CpMixMass_);
	}

	return 0;
}

int BatchReactorHomogeneousConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
	return 0;
}
