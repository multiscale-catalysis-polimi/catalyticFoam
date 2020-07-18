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

namespace OpenSMOKE
{
	const char *const NonLinearFunctionSolver_Robust::OPENSMOKE_ERROR = "NonLinearFunctionSolver_Robust";
	const int NonLinearFunctionSolver_Robust::MAX_ITER = 500;
	const double NonLinearFunctionSolver_Robust::TOL_REL = 100. * OpenSMOKE::MachEps();
	const double NonLinearFunctionSolver_Robust::TOL_ABS = 1.e-30;
	const double NonLinearFunctionSolver_Robust::DT_RATIO = .1;
	const double NonLinearFunctionSolver_Robust::DT_MIN = .1;
	const double NonLinearFunctionSolver_Robust::TOL_REL_Y = sqrt(OpenSMOKE::MachEps());
	const double NonLinearFunctionSolver_Robust::TOL_ABS_Y = 1.e-30;
	const double NonLinearFunctionSolver_Robust::U_F_RATIO = 0.1;

	void NonLinearFunctionSolver_Robust::InitializeMono(void)
	{
		lambdaMax = -OPENSMOKE_BIG;
		iterTotal = 0;
		iterTotalFeasible = 0;
		iterTotalInConstructor = 0;
		iterFeasibleInConstructor = 0;
		iter = 0;
		nCritic = 1;
		control = -1;
		maxIter = MAX_ITER;
		unimodal = 0;
		tTolRel = TOL_REL;
		tTolAbs = TOL_ABS;
		yTol = TOL_ABS_Y;
		dtTol = .1 * DT_MIN;
		oneMunfeasibleFeasibleRatio = 1. - unfeasibleFeasibleRatio;

		if(tMin > tMax)
			Swap(&tMin,&tMax);
		if(t0 < tMin || t0 > tMax)
			OpenSMOKE::FatalErrorMessage("NonLinearFunctionSolver_Robust Class: Wrong interval specified by the user.");
		ChangeDimensions(0,&y, true);
		ChangeDimensions(0,&t, true);
		ChangeDimensions(0,&tCriticInf, true);
		ChangeDimensions(0,&tCriticSup, true);
		numPoints = 1;
		t.Insert(1,t0);
		y.Insert(1,y0);
		if(tMin == tMax)
		{
			control = 1;
			tCriticInf.Insert(1,tMin);
			tCriticSup.Insert(1,tMax);
			return;
		}

		if(t0 == tMin)
		{
			tCriticInf.Insert(1,tMin);
			tNew = tMax;
			unfeasible_ = 0;
			iterTotal++;
			yNew = Y(tNew);
			if(unfeasible_ == 1)
			{
				FindtMax(t0,tNew);
				tCriticSup.Insert(1,tMax);
			}
			else
			{
				iterTotalFeasible++; iter++;
				numPoints++;
				if(yNew < ySolution)
				{
					tSolution = tNew;
					ySolution = yNew;
				}
				t.Insert(2,tNew);
				y.Insert(2,yNew);
				tCriticSup.Insert(1,tMax);
				tNew = .5 * (tMin + tMax);
				unfeasible_ = 0;
				iterTotal++;
				yNew = Y(tNew);
				if(unfeasible_ == 1)
				{
					FindtCritic(tMin,tNew,tMax);
				}
				else
				{
					iterTotalFeasible++; iter++;
					numPoints++;
					if(yNew < ySolution)
					{
						tSolution = tNew;
						ySolution = yNew;
					}
					t.Insert(2,tNew);
					y.Insert(2,yNew);
				}
			}
		}
		else if(t0 > tMin && t0 < tMax)
		{
			tNew = tMax;
			unfeasible_ = 0;
			iterTotal++;
			yNew = Y(tNew);
			if(unfeasible_ == 1)
			{
				FindtMax(t0,tNew);
				tCriticSup.Insert(1,tMax);
			}
			else
			{
				iterTotalFeasible++; iter++;
				tCriticSup.Insert(1,tMax);
				numPoints++;
				if(yNew < ySolution)
				{
					tSolution = tNew;
					ySolution = yNew;
				}
				t.Insert(2,tNew);
				y.Insert(2,yNew);
			}
			tNew = tMin;
			unfeasible_ = 0;
			iterTotal++;
			yNew = Y(tNew);
			if(unfeasible_ == 1)
			{
				FindtMin(tNew,t0); // TODO
				tCriticInf.Insert(1,tMin);
			}
			else
			{
				iterTotalFeasible++; iter++;
				tCriticInf.Insert(1,tMin);
				numPoints++;
				if(yNew < ySolution)
				{
					tSolution = tNew;
					ySolution = yNew;
				}
				t.Insert(1,tNew);
				y.Insert(1,yNew);
			}
		}
		else
		{
			tCriticSup.Insert(1,tMax);
			tNew = tMin;
			unfeasible_ = 0; 
			iterTotal++;
			yNew = Y(tNew);
			if(unfeasible_ == 1)
			{
				FindtMin(tNew,tMax);
				tCriticInf.Insert(1,tMin);
			}
			else
			{
				iterTotalFeasible++; iter++;
				tCriticInf.Insert(1,tMin);
				numPoints++;
				if(yNew < ySolution)
				{
					tSolution = tNew;
					ySolution = yNew;
				}
				t.Insert(1,tNew);
				y.Insert(1,yNew);
				tNew = .5 * (tMin + tMax);
				unfeasible_ = 0; 
				iterTotal++;
				yNew = Y(tNew);
				if(unfeasible_ == 1)
				{
					FindtCritic(tMin,tNew,tMax);
				}
				else
				{
					iterTotalFeasible++; iter++;
					numPoints++;
					if(yNew < ySolution)
					{
						tSolution = tNew;
						ySolution = yNew;
					}
					t.Insert(2,tNew);
					y.Insert(2,yNew);
				}
			}
		}

		iterTotalInConstructor = iterTotal;
		iterFeasibleInConstructor = iterTotalFeasible;
	}

	double NonLinearFunctionSolver_Robust::Y(const double t)
	{
		return std::fabs(yMono(t));
	}

	NonLinearFunctionSolver_Robust::NonLinearFunctionSolver_Robust()
	{
		unfeasibleFeasibleRatio = U_F_RATIO;
		mono = 0;
		lambdaMax = -OPENSMOKE_BIG;
	}

	NonLinearFunctionSolver_Robust::NonLinearFunctionSolver_Robust(const NonLinearFunctionSolver_Robust &rval)
	{
		OpenSMOKE::FatalErrorMessage("NonLinearFunctionSolver_Robust Class: copy constructor not available.");
	}

	NonLinearFunctionSolver_Robust::NonLinearFunctionSolver_Robust(const double tt0, const double yy0, double (*yy)(double t), const double tmin, const double tmax)
	{
		const double yy00 = std::fabs(yy0);
		unfeasibleFeasibleRatio = U_F_RATIO;
		mono = 0;
		(*this)(tt0,yy00,yy,tmin,tmax);
	}

	void NonLinearFunctionSolver_Robust::operator()(const double tt0, const double yy0, double (*yy)(double t), const double tmin, const double tmax)
	{
		yMono = yy;
		tSolution = t0 = tt0;
		ySolution = y0 = std::fabs(yy0);
		tMin = tmin;
		tMax = tmax;
		InitializeMono();
	}

	NonLinearFunctionSolver_Robust::NonLinearFunctionSolver_Robust(const double tt0, double (*yy)(double t), const double tmin, const double tmax)
	{
		unfeasibleFeasibleRatio = U_F_RATIO;
		mono = 0;
		(*this)(tt0,yy,tmin,tmax);
	}

	void NonLinearFunctionSolver_Robust::operator()(const double tt0, double (*yy)(double t), const double tminimum, const double tmaximum)
	{
		double tmin = tminimum;
		double tmax = tmaximum;
		int i;
		int itt,itf;
		itt = itf = 0;
		yMono = yy;
		unfeasible_ = 0;
		t0 = tt0;
		y0 = std::fabs(yMono(t0));itt++;
		if(tmin > tmax)
			Swap(&tmin,&tmax);
		if(unfeasible_ != 0)
		{
			t0 = tmin;
			unfeasible_ = 0;
			y0 = std::fabs(yMono(t0));itt++;
			if(unfeasible_ != 0)
			{
				int n = 10000;
				double dt = tmax - tmin;
				if(dt == 0.)
				{
					std::cout << "Warning Message: NonLinearFunctionSolver_Robust Class: tMin == tMax" << std::endl;
					return;
				}
				dt /= double(n);
				for(i = 1;i < n;i++)
				{
					tmin = t0;
					t0 +=  dt;
					unfeasible_ = 0;
					y0 = std::fabs(yMono(t0));itt++;

					if(unfeasible_ != 0)
						continue;
					itf++;
					break;
				}
				if(i > 10000)
				{
					std::cout << "Warning Message: NonLinearFunctionSolver_Robust Class: No feasible region has been located" << std::endl;
					return;
				}
				i++;
			}
			else
			{
				itf++;
				i = 2;
			}
		}
		else
		{
			itf++;
		}
		tSolution = t0;
		ySolution = y0;
		tMin = tmin;
		tMax = tmax;
		InitializeMono();
		iterTotal += itt;
		iterTotalFeasible += itf;
		iterFeasibleInConstructor += itf;
		iterTotalInConstructor += itt;
	}

	NonLinearFunctionSolver_Robust::NonLinearFunctionSolver_Robust(double (*yy)(double t), const double tmin, const double tmax)
	{
		unfeasibleFeasibleRatio = U_F_RATIO;
		mono = 0;
		(*this)(.5*(tmin + tmax),yy,tmin,tmax);
	}

	void NonLinearFunctionSolver_Robust::operator()(double (*yy)(double t), const double tmin, const double tmax)
	{
		(*this)(.5*(tmin + tmax),yy,tmin,tmax);
	}


	char NonLinearFunctionSolver_Robust::operator()(int ni)
	{
		// control < 0 : continue
		// control = 0 : found ySolution < di yTol
		// control > 0 : calculations over
		if(control >= 0)
			return control;
		control = -1;
		iter = 0;
		nIter = ni;
		return MonoSearch();
	}

	char NonLinearFunctionSolver_Robust::MonoSearch(void)
	{
		double lalala;
		int iCritic,imin,imax,i1,i2,i3,j;
		double delta,deltam;
		double dtReduced;

		restart_critic:
		control = -1;
		for(iCritic = 1;iCritic <= nCritic;iCritic++)
		{
			imin = 1 + t.LocateInFirstNSortedElements(numPoints,tCriticInf[iCritic]);
			imax = 1 + t.LocateInFirstNSortedElements(numPoints,tCriticSup[iCritic]);


			if(imin >= imax - 1)
			{
				if(imin == imax)
					continue;
				else
				{
					i1 = imin;
					i2 = imin + 1;
					t1 = t[i1];
					t2 = t[i2];
					delta = tTolAbs + tTolRel * std::fabs(t1);
					dt21 = t2 - t1;
					if(dt21 > delta)
					{
						tNew = .5 * (t1 + t2);
						unfeasible_ = 0; 
						iterTotal++;
						yNew = Y(tNew);

						if(unfeasible_ == 1)
						{
							FindtCritic(t1,tNew,t2);
							goto restart_critic;
						}
						else
						{
							iterTotalFeasible++; iter++;
							if(yNew < ySolution)
							{
								tSolution = tNew;
								ySolution = yNew;
							}
							j = t.InsertElementInFirstNSortedElements(numPoints,tNew);
							y.Insert(j,yNew);
							numPoints++;

							if(ySolution < yTol)
							{
								control = 0;
								return control;
							}
							if(iter >= nIter)
							{
								control = -1;
								return control;
							}
							goto restart_critic;
						}
					}

					else if(unimodal == 1 && y[i1] < y[i2] && y[i1] == ySolution)
					{
						control = 2;
						return control;
					}
				}
			}
			else
			{
				i1 = imin;
				i2 = imin + 1;
				t1 = t[i1];	y1 = y[i1];
				t2 = t[i2];	y2 = y[i2];
				delta = tTolAbs + tTolRel * std::fabs(t1);
				dt21 = t2 - t1;

				if(dt21 > delta && y1 <= y2 && y2 > y1 + TOL_REL_Y * (std::fabs(y1) + std::fabs(y2)))
				{
					deltam = delta * .965;
					tNew = t1 + deltam;
					unfeasible_ = 0; 
					iterTotal++;
					yNew = Y(tNew);

					if(unfeasible_ == 1)
					{
						FindtCritic(t1,tNew,t2);
						goto restart_critic;
					}
					else
					{
						iterTotalFeasible++; iter++;
						if(yNew < ySolution)
						{
							tSolution = tNew;
							ySolution = yNew;
						}
					
						j = t.InsertElementInFirstNSortedElements(numPoints,tNew);
						y.Insert(j,yNew);
						numPoints++;

						if(ySolution < yTol)
						{
							control = 0;
							return control;
						}
						if(iter >= nIter)
						{
							control = -1;
							return control;
						}
						goto restart_critic;
					}
				}
				else if(unimodal == 1 && y1 <= y2  && y1 == ySolution)
				{
					control = 2;
					return control;
				}

				i2 = imax - 1;
				i3 = imax;
				t2 = t[i2];	y2 = y[i2];
				t3 = t[i3];	y3 = y[i3];
				delta = tTolAbs + tTolRel * std::fabs(t3);
				dt32 = t3 - t2;

				if(dt32 > delta && y3 <= y2 && y2 > y3 + TOL_REL_Y * (std::fabs(y2) + std::fabs(y3)))
				{
					deltam = delta * .965;
					tNew = t3 - deltam;
					unfeasible_ = 0;
					iterTotal++;
					yNew = Y(tNew);
					if(unfeasible_ == 1)
					{
						FindtCritic(t2,tNew,t3);
						goto restart_critic;
					}
					else
					{
						iterTotalFeasible++; iter++;
						if(yNew < ySolution)
						{
							tSolution = tNew;
							ySolution = yNew;
						}
						j = t.InsertElementInFirstNSortedElements(numPoints,tNew);
						y.Insert(j,yNew);
						numPoints++;

						if(ySolution < yTol)
						{
							control = 0;
							return control;
						}
						if(iter >= nIter)
						{
							control = -1;
							return control;
						}
						goto restart_critic;
					}
				}
				else if(unimodal == 1 && y3 <= y2  && y3 == ySolution) 
				{
					control = 2;
					return control;
				} 

				for(i2 = imin + 1;i2 < imax;i2++)
				{
					i1 = i2 - 1;
					i3 = i2 + 1;
					t1 = t[i1];	y1 = y[i1];
					t2 = t[i2];	y2 = y[i2];
					t3 = t[i3];	y3 = y[i3];
					delta = tTolAbs + tTolRel * std::fabs(t2);
					deltam = delta * .965;
					dt21 = t2 - t1; dt32 = t3 - t2;
					if(dt21 <= delta && dt32 <= delta)
					{
						if(unimodal == 1 && y[i2] <= y[i1] && y[i2] <= y[i3] && y[i2] == ySolution) 
						{
							control = 2;
							return control;
						} 
						continue;
					}
				
					dt31 = t3 - t1;
					dy21 = y2 - y1; dy32 = y3 - y2;
					if(dy21 >= 0. || dy32 < 0.)
						continue;
					if(dt32 == 0. || dt31 == 0.)
						lambda = 0.;
					else
						lambda = (dy32 * dt21 - dy21 * dt32) / (dt32 * dt31);
					if(lambda <= 0.)
						continue;
					if(dt21 != 0.)
					{
						lalala = lambda / dt21;
						if(lalala > lambdaMax)
							lambdaMax = lalala;
					}
					if(std::fabs(dy21) < TOL_REL_Y * (std::fabs(y1) + std::fabs(y2)) && std::fabs(dy32) < TOL_REL_Y * (std::fabs(y2) + std::fabs(y3)))
					{
						if(unimodal == 1) 
						{
							control = 2;
							return control;
						} 
						continue;
					}

					tNew = .5 * (t2 + t1 - dy21 / lambda);

					if(tNew <= t1 || tNew >= t3)
						continue;
					if(tNew <= t2)
					{
						if(dt21 <= delta)
						{
							tNew = t2 + .5 * dt32;
						}
						else if(dt21 < dt32)
						{
							tNew = t2 + dt32 * DT_RATIO;
						}
						else
						{
							dtReduced = dt21 * DT_RATIO;
							if(dtReduced < delta)
								dtReduced = deltam;
							if(tNew < t1 + dtReduced)
								tNew = t1 + dtReduced;
							if(tNew > t2 - dtReduced)
								tNew = t2 - dtReduced;
						}
					}
					else
					{
						if(dt32 <= delta)
						{
							tNew = t2 - dt21 * DT_RATIO *.5 ;
						}
						else if(dt32 < dt21)
						{
							tNew = t2 - dt21 * DT_RATIO;
						}
						else
						{
							dtReduced = dt32 * DT_RATIO;
							if(dtReduced < delta)
								dtReduced = deltam;
							if(tNew < t2 + dtReduced)
								tNew = t2 + dtReduced;
							if(tNew > t3 - dtReduced)
								tNew = t3 - dtReduced;
						}
					}
					unfeasible_ = 0;
					iterTotal++;
					yNew = Y(tNew);

					if(unfeasible_ == 1)
					{
						if(tNew < t2)
							FindtCritic(t1,tNew,t2);
						else
							FindtCritic(t2,tNew,t3);
						goto restart_critic;
					}
					else
					{
						iterTotalFeasible++; iter++;
						if(yNew < ySolution)
						{
							tSolution = tNew;
							ySolution = yNew;
						}
						j = t.InsertElementInFirstNSortedElements(numPoints,tNew);
						y.Insert(j,yNew);
						numPoints++;

						if(ySolution < yTol)
						{
							control = 0;
							return control;
						}
						if(iter >= nIter)
						{
							control = -1;
							return control;
						}
						goto restart_critic;
					}
				}
			}
		}

		ChangeDimensions(0,&tAux, true);
		ChangeDimensions(0,&yAux, true);
		ChangeDimensions(0,&tAux1, true);
		if(dtTol == DT_MIN && y.Size() < 7)
			dtTol *= .5;
		dtReduced = dtTol * (tMax - tMin);
		j = 1;
		for(iCritic = 1;iCritic <= nCritic;iCritic++)
		{
			imin = 1 + t.LocateInFirstNSortedElements(numPoints,tCriticInf[iCritic]);
			imax = 1 + t.LocateInFirstNSortedElements(numPoints,tCriticSup[iCritic]);
		
			for(i2 = imin + 1;i2 <= imax;i2++)
			{
				i1 = i2 - 1;
				t1 = t[i1];
				t2 = t[i2];
				dt21 = t2 - t1;
				if(dt21 > dtReduced)
				{
					tNew = .5 * (t1 + t2);
					unfeasible_ = 0;
					iterTotal++;
					yNew = Y(tNew);

					if(unfeasible_ == 1)
					{
						FindtCritic(t1,tNew,t2);
						goto restart_critic;
					}
					else
					{
						iterTotalFeasible++; iter++;
						if(yNew < ySolution)
						{
							tSolution = tNew;
							ySolution = yNew;
						}
						tAux.Insert(j,tNew);
						yAux.Insert(j++,yNew);
					}
				}

				if(ySolution < yTol)
				{
					for(int i = 1;i <= tAux.Size();i++)
					{
						j = t.InsertElementInFirstNSortedElements(numPoints,tAux[i]);
						y.Insert(j,yAux[i]);
						numPoints++;
					}
					control = 0;
					return control;
				}
				if(iter >= nIter)
				{
					for(int i = 1;i <= tAux.Size();i++)
					{
						j = t.InsertElementInFirstNSortedElements(numPoints,tAux[i]);
						y.Insert(j,yAux[i]);
						numPoints++;
					}
					control = -1;
					return control;
				}
			}
		}

		if(tAux.Size() == 0)
		{
			control = 1;
			return control;
		}
		else
		{
			for(int i = 1;i <= tAux.Size();i++)
			{
				j = t.InsertElementInFirstNSortedElements(numPoints,tAux[i]);
				y.Insert(j,yAux[i]);
				numPoints++;
			}
			goto restart_critic;
		}

		control = 1;
		return control;
	}

	void NonLinearFunctionSolver_Robust::FindtMax(double tt1,double ttUnfeasible)
	{
		int j;
		double tf,tu,delta;
		double ttmax = tMax;
		tf = tt1;
		tu = ttUnfeasible;
		while(1)
		{

			tNew = unfeasibleFeasibleRatio * tf + oneMunfeasibleFeasibleRatio * tu;
			if(tNew == tf || tNew == tu)
				break;

			delta = tTolAbs + tTolRel * std::fabs(tNew);
			delta *= .1;
			if(tu - tf <= delta)
				break;
			if(iterTotalFeasible >= maxIter)
				break;

			unfeasible_ = 0;
			iterTotal++;
			yNew = Y(tNew);
			if(unfeasible_ == 1)
			{
				tu = tNew;
			}
			else
			{
				tMax = tf = tNew;
				iterTotalFeasible++; iter++;
				if(yNew < ySolution)
				{
					tSolution = tNew;
					ySolution = yNew;
				}
				j = t.InsertElementInFirstNSortedElements(numPoints,tNew);
				y.Insert(j,yNew);
				numPoints++;
			}
		}
		if(ttmax == tMax)
		{
			tMax = tt1;
		}
	}

	void NonLinearFunctionSolver_Robust::FindtMin(double ttUnfeasible,double ttn)
	{
		int j;
		double tf,tu,delta;
		double ttmin = tMin;
		tf = ttn;
		tu = ttUnfeasible;
		while(1)
		{

			tNew = unfeasibleFeasibleRatio * tf + oneMunfeasibleFeasibleRatio * tu;
			if(tNew == tf || tNew == tu)
				break;

			delta = tTolAbs + tTolRel * std::fabs(tNew);
			delta *= .1;
			if(tf - tu <= delta)
				break;
			if(iterTotalFeasible >= maxIter)
				break;
			unfeasible_ = 0;
			iterTotal++;
			yNew = Y(tNew);
			if(unfeasible_ == 1)
			{
				tu = tNew;
			}
			else
			{
				tMin = tf = tNew;
				iterTotalFeasible++; iter++;
				if(yNew < ySolution)
				{
					tSolution = tNew;
					ySolution = yNew;
				}

				j = t.InsertElementInFirstNSortedElements(numPoints,tNew);
				y.Insert(j,yNew);
				numPoints++;
			}
		}
		if(ttmin == tMin)
		{
			tMin =  ttn;
		}
	}

	void NonLinearFunctionSolver_Robust::FindtCritic(const double tt1, const double ttUnfeasible, const double ttn)
	{
		double tMinMemo = tMin;
		double tMaxMemo = tMax;
		FindtMax(tt1,ttUnfeasible);
		tCriticSup.InsertElementInSortedVector(tMax);
		tMax = tMaxMemo;
		FindtMin(ttUnfeasible,ttn);
		tCriticInf.InsertElementInSortedVector(tMin);
		nCritic++;
		tMin = tMinMemo;
	}

	void NonLinearFunctionSolver_Robust::GetSolutions(OpenSMOKE::OpenSMOKEVectorDouble *tSolutions, OpenSMOKE::OpenSMOKEVectorDouble *ySolutions)
	{
		int i,iCritic,imin,imax;
		ChangeDimensions(0,tSolutions, true);
		ChangeDimensions(0,ySolutions, true);
		for(iCritic = 1;iCritic <= nCritic;iCritic++)
		{
			imin = 1 + t.LocateInFirstNSortedElements(numPoints,tCriticInf[iCritic]);
			imax = 1 + t.LocateInFirstNSortedElements(numPoints,tCriticSup[iCritic]);
			if(imin >= imax - 1)
			{
				for(i = imin;i <= imax;i++)
				{
					tSolutions->Append(t[i]);
					ySolutions->Append(y[i]);
				}
			}
			else
			{
				for(i = imin;i <= imax;i++)
				{
					if(i == imin)
					{
						if(y[i] <= y[i + 1])
						{
							tSolutions->Append(t[i]);
							ySolutions->Append(y[i]);
						}
					}
					else if(i == imax)
					{
						if(y[i] <= y[i - 1])
						{
							tSolutions->Append(t[i]);
							ySolutions->Append(y[i]);
						}
					}
					else
					{
						if(y[i] <= y[i - 1] && y[i] <= y[i + 1])
						{
							tSolutions->Append(t[i]);
							ySolutions->Append(y[i]);
						}
					}
				}
			}
		}
	}

	void NonLinearFunctionSolver_Robust::SetTolRel(const double tr)
	{
		control = -1;
		if(tr <= 0.)
			tTolRel = TOL_REL;
		else
			tTolRel = std::fabs(tr);
	}
	void NonLinearFunctionSolver_Robust::SetTolAbs(const double ta)
	{
		control = -1;
		if(ta <= 0.)
			tTolAbs = TOL_ABS;
		else
			tTolAbs = std::fabs(ta);
	}
	void NonLinearFunctionSolver_Robust::SetTolY(const double yt)
	{
		control = -1;
		if(yt <= 0.)
			yTol = TOL_ABS_Y;
		else
		 yTol = yt;
	}

	void NonLinearFunctionSolver_Robust::SetMaxIter(const int mxit)
	{
		if(mxit <= 0)
			maxIter = MAX_ITER;
		else
			maxIter = mxit;
	}

}