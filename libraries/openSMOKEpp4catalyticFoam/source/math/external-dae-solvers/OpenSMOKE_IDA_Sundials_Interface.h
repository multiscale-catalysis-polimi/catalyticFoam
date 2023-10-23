/*-----------------------------------------------------------------------*\
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
|   Copyright(C) 2016, 2015 Alberto Cuoci                                 |
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


#define COMPLETE_DAESOLVERINTERFACE_IDA_Sundials(classname)\
bool classname::instanceFlag = false;\
classname* classname::m_this = NULL;\
unsigned int classname::n_; \
bool* classname::algebraic_;

#define DEFINE_DAESOLVERINTERFACE_IDA_Sundials(classname)\
private:\
	\
	static bool instanceFlag;\
	static classname* m_this;\
	static unsigned int n_;\
	static bool* algebraic_;\
	classname()\
	{\
	}\
	\
public:\
	   \
	static classname* GetInstance(const unsigned int n, const bool* algebraic)\
	{\
		if(instanceFlag == false)\
		{\
			m_this = new classname();\
			instanceFlag = true;\
			m_this->n_ = n; \
			m_this->algebraic_ = new bool[m_this->n_]; \
			for(unsigned int i=0;i<m_this->n_;i++)\
				m_this->algebraic_[i] = algebraic[i];\
			return m_this;\
		}\
		else\
		{\
			return m_this;\
		}\
	}\
	\
    ~classname()\
	{\
		instanceFlag = false;\
	}\
	\
	static int GetSystemFunctionsStatic(realtype t, N_Vector y, N_Vector yp, N_Vector f, void *user_data)\
	{\
		return m_this->GetSystemFunctionsCallBack(t,y,yp,f,user_data);\
	}\
	\
	static int GetAnalyticalJacobianStatic(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\
	{\
		return m_this->GetAnalyticalJacobianCallBack(N,t,y,fy,J,user_data,tmp1, tmp2, tmp3);\
	}\
	\
	static int GetWriteFunctionStatic(realtype t, N_Vector y)\
	{\
		return m_this->GetWriteFunctionCallBack(t,y);\
	}\
	\
	int GetSystemFunctionsCallBack(realtype t, N_Vector y, N_Vector yp, N_Vector f, void *user_data)\
	{\
		double* ydata  = N_VGetArrayPointer(y);\
		double* ypdata = N_VGetArrayPointer(yp);\
		double* fdata  = N_VGetArrayPointer(f);\
		int flag = GetSystemFunctions(t, ydata, fdata);\
		for(unsigned int i=0;i<m_this->n_;i++)\
			if (m_this->algebraic_[i]==false)	fdata[i] -= ypdata[i];\
		return(flag);\
	}\
	\
	int GetAnalyticalJacobianCallBack(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\
	{\
		return(0);\
	}\
	\
	int GetWriteFunctionCallBack(realtype t, N_Vector y)\
	{\
		double* ydata = N_VGetArrayPointer(y);\
		int flag = GetWriteFunction(t, ydata);\
		return flag; \
	}\

