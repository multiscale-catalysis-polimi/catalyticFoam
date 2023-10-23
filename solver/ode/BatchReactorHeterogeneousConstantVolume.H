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

#ifndef BatchReactorHeterogeneousConstantVolume_H
#define	BatchReactorHeterogeneousConstantVolume_H

#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"

//!  A class to solve constant volume batch reactors with heterogeneous reactions
/*!
	This class provides the equations governing a batch reactor with homogeneous and
        heterogeneous reactions
*/

class BatchReactorHeterogeneousConstantVolume
{
public:

	/**
	* Default constructor
	*/
	BatchReactorHeterogeneousConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap,
							OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN& thermodynamicsSurfaceMap, 
							OpenSMOKE::KineticsMap_Surface_CHEMKIN& kineticsSurfaceMap);
	/**
	*@brief Sets the volume, the geometric area, and the ratio between catalytic and geometric area
	*/
	void SetReactor( const double V0, const double A, const double alfaCatalyst);
	
	/**
	*@brief Turns on or off the energy equation
	*/
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }

	/**
	*@brief Turns on or off the homogeneous reactions
	*/
	void SetHomogeneousReactions(const bool flag) { homogeneousReactions_ = flag; }

	/**
	*@brief Turns on or off the reaction heat from surface reaction in the energy equation
	*/
	void SetReactionHeatFromHeterogeneousReactions(const bool flag) { reactionHeatFromHeterogeneousReactions_ = flag; }
	
	/**
	*@brief Turns on or off the control of minimum values of mass fractions of species
	*/	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

    /**
	*@brief Turns on or off the use of ML surrogates for the TOF
	*/	
	void SetMLSurrogate( const bool flag ) { MLSurrogate_ = flag; };

	/**
	*@brief Returns the total number of equations
	*/
	unsigned int NumberOfEquations() { return NE_; }

	/**
	*@brief Returns the governing equations
	*/
	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	/**
	*@brief This function is called at the end of each internal time step
	*/
	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

private:

	double T_;
	double A_;
	double alfaCatalyst_;
	double V0_;
	double MW_;
	double mass_;

	unsigned int NC_;
	unsigned int SURF_NP_;
	unsigned int SURF_NC_;
	unsigned int NE_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&			thermodynamicsMap_;				//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN&					kineticsMap_;					//!< kinetic map
	OpenSMOKE::ThermodynamicsMap_Surface_CHEMKIN&	thermodynamicsSurfaceMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_Surface_CHEMKIN&			kineticsSurfaceMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble MWs_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromGas_;
	OpenSMOKE::OpenSMOKEVectorDouble RfromSurface_;
	OpenSMOKE::OpenSMOKEVectorDouble Z_;
	OpenSMOKE::OpenSMOKEVectorDouble Gamma_;
	OpenSMOKE::OpenSMOKEVectorDouble Rsurface_;
	OpenSMOKE::OpenSMOKEVectorDouble RsurfacePhases_;
	OpenSMOKE::OpenSMOKEVectorDouble dummy;
	
	double QRGas_;
	double QRSurface_;
	
	bool checkMassFractions_;
	bool energyEquation_;
	bool reactionHeatFromHeterogeneousReactions_;
	bool homogeneousReactions_;
    bool MLSurrogate_;
};

#include "BatchReactorHeterogeneousConstantVolume.C"

#endif // BatchReactorHeterogeneousConstantVolume_H

