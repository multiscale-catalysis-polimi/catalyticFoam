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

#ifndef BatchReactorHomogeneousConstantPressure_H
#define	BatchReactorHomogeneousConstantPressure_H

#include "OpenSMOKEpp"
#include "maps/Maps_CHEMKIN"

//!  A class to solve constant pressure batch reactors with homogeneous reactions
/*!
	This class provides the equations governing a batch reactor with homogeneous reactions
*/

class BatchReactorHomogeneousConstantPressure
{
public:

	/**
	* Default constructor
	*/
	BatchReactorHomogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	/**
	*@brief Sets the pressure
	*/
	void SetReactor( const double P0 );
	
	/**
	*@brief Turns on or off the energy equation
	*/
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	/**
	*@brief Turns on or off the control of minimum values of mass fractions of species
	*/	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

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
	double P0_;
	double rho_;

	double cTot_;
	double MW_;
	double CpMixMass_;

	unsigned int NC_;
	unsigned int NE_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	
	bool checkMassFractions_;
	bool energyEquation_;
};

#endif // BatchReactorHomogeneousConstantPressure_H

