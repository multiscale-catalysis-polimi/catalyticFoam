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

// ODE Solver (constant pressure)
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_OpenSMOKE> denseHomogeneousOdeConstantPressure;
typedef OdeSMOKE::MethodGear<denseHomogeneousOdeConstantPressure> methodGearHomogeneousConstantPressure;
OdeSMOKE::MultiValueSolver<methodGearHomogeneousConstantPressure> odeHomogeneousSolverConstantPressure;
odeHomogeneousSolverConstantPressure.SetReactor(&batchReactorHomogeneousConstantPressure);

// ODE Solver (constant volume)
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_OpenSMOKE> denseHomogeneousOdeConstantVolume;
typedef OdeSMOKE::MethodGear<denseHomogeneousOdeConstantVolume> methodGearHomogeneousConstantVolume;
OdeSMOKE::MultiValueSolver<methodGearHomogeneousConstantVolume> odeHomogeneousSolverConstantVolume;	
odeHomogeneousSolverConstantVolume.SetReactor(&batchReactorHomogeneousConstantVolume);

#if OPENSMOKE_USE_DVODE == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE *dvodeHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE::GetInstance();
	dvodeHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE *dvodeHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE::GetInstance();
	dvodeHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_ODEPACK == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE *dlsodeHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE::GetInstance();
	dlsodeHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE *dlsodeHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE::GetInstance();
	dlsodeHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA *dlsodaHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA::GetInstance();
	dlsodaHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA *dlsodaHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA::GetInstance();
	dlsodaHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE *cvodeHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE::GetInstance();
	cvodeHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE *cvodeHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE::GetInstance();
	cvodeHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_DASPK == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK *daspkHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK *daspkHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_RADAU == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5 *radauHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radauHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5::GetInstance();
	radauHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5 *radauHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radauHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5::GetInstance();
	radauHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_MEBDF == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF *mebdfHomogeneousSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfHomogeneousSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF::GetInstance();
	mebdfHomogeneousSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF *mebdfHomogeneousSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfHomogeneousSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF::GetInstance();
	mebdfHomogeneousSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif
