/*-----------------------------------------------------------------------*\
|   Main Author: Mauro Bracconi                                           |
|                                                                         |
|   Contacts: Mauro Bracconi                                              |
|   email: mauro.bracconi@polimi.it                                       |
|   Department of Energy                                                  |
|   Politecnico di Milano                                                 |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of ISATLib library.                                 |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2014 Mauro Bracconi, Alberto Cuoci, Matteo Maestri       |
|   ISATLib is free software: you can redistribute it and/or modify       |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   ISATLib is distributed in the hope that it will be useful,            |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with ISATLib. If not, see <http://www.gnu.org/licenses/>.       |
|                                                                         |
\*-----------------------------------------------------------------------*/

/*------------------------------------------------------------------------\
|                                                                         |
|   chemComp class:                                                       |
|                                                                         |
|   Binary Tree leaf instance                                             |
|   Procedure to set up, modify EOA of the leaf                           |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "chemComp.h"
#include <iostream>
/*
 * chemComp: costructor
 * The costructor create the leaf and inizialize the EOA
 */
chemComp::chemComp(const unsigned int index, const VectorXd &phi, const VectorXd &Rphi, const MatrixXd &A, const VectorXd &scaleFactor, const double epsTol, const unsigned int nSpec, const unsigned int qrType, binaryNode *nodeParent) { 
    index_ = index;
    phi_ = phi;
    Rphi_ = Rphi;
    node_ = nodeParent;
    mapGrad = A;
    scaleFactor_ = scaleFactor;
    epsTol_ = epsTol;
    nSpec_ = nSpec;

    // inizialize counter grown
    nGrown_ = 0;
    // inizialize counter used
    nUsed_ = 0;
    
    // set creation time
    lastTimeUse_ = double(clock())/CLOCKS_PER_SEC;
  
    // creation of the matrix scale factor
    B_ = scaleFactor.asDiagonal();
    
    // set qrType
    qrType_ = qrType;
    
    // creation of the EOA and rmin/rmax
    LT_.setZero(nSpec_,nSpec_);
    inizializeEOA(A);
    
}

chemComp::chemComp(chemComp& phi0) { 
    index_ = phi0.getIndex();
    phi_ = phi0.getPhi();
    Rphi_ = phi0.getRphi();
    node_ = phi0.getNode();
    mapGrad = phi0.getMapGrad();
    scaleFactor_ = phi0.getScaleFactor();
    epsTol_ = phi0.getEpsTol();
    nSpec_ = phi0.getSpec();
    
    nGrown_ = phi0.getGrown();
    nUsed_ = phi0.getUsed();
    
    LT_ = phi0.getLT();
    rmin_ = phi0.getRmin();
    rmax_ = phi0.getRmax();
    
    qrType_=phi0.getQRType();
    
    B_ = phi0.getScaleMatrix();
       
}
/*
 * growEOA: modify and update the EOA ellipsoid and the semaxes rmin, rmax
 * input:
 *      - phiQ: query point
 * 
 * First step is to compute the vector phiQ in the space of the EOA:
 * p' = LT*(phiQ-phi)
 * Next step is to compute the rank-one decomposition:
 * G = I + gamma.p'.p'^T, with gamma = (1/|p'|-1)*1/|p'|^2 
 * Next step is to compute L':
 * L'L'^T = (L.G)(L.G)^T, 
 * Last step is to compute the QR decomposition of the new matrix and to save 
 * matrix R which is L^T
 * and to update che rmin and rmax using the SVD of that matrix
 */
void chemComp::growEOA(const VectorXd &phiQ) {
    // p'    
    VectorXd phiQtilde = getLT()*(phiQ - getPhi());
    MatrixXd G = getLT().transpose();
    // rank-one decomposition
    double gamma = 0.0;
    double module = phiQtilde.norm();
    
    gamma = (1./module - 1.)/(std::pow(module,2.));
  
    //modify matrix
    MatrixXd LG = G+(gamma*G*phiQtilde)*phiQtilde.transpose();
  
    // compute QR decomposition
    // L^T = R
    if (qrType_ == 0)
    {
        Eigen::FullPivHouseholderQR<MatrixXd> qr(LG.transpose());
        LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
    } 
    else if (qrType_ == 1)
    {
        Eigen::ColPivHouseholderQR<MatrixXd> qr(LG.transpose());
        LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
    } 
    else if (qrType_ == 2)
    {
        Eigen::HouseholderQR<MatrixXd> qr(LG.transpose());
        LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
    }
   
    // recompute the ellipsoids axes in the way to use it in the next EOA control
    if (nSpec_< 100000000)
	{
		Eigen::JacobiSVD<MatrixXd> svd(LT_);
		VectorXd diag = svd.singularValues();
		rmin_ = 1./diag(0);
		rmax_ = 1./diag(nSpec_-1);
	}
	else
	{
		 // disable floating point traps for BDCSVD
	        int fpe = fegetexcept();
	        fedisableexcept(FE_ALL_EXCEPT);
                                                                                                                 	
		Eigen::BDCSVD<MatrixXd> svd(LT_);
		VectorXd diag = svd.singularValues();
		
		// enable floating point traps after BDCSVD
		feenableexcept(fpe);
		
		rmin_ = 1./diag(0);
		rmax_ = 1./diag(nSpec_-1);
	}
    
    nGrown_++;
}

/*
 * inEOA : check if the query is inside the EOA 
 * input:
 *      - phiQ: query point
 * 
 * The function evaluates the distance of the query point from the center 
 * dphi = || phiQ - phi0 || and check if this distance is larger or smaller 
 * than the max and min semiaxes of the ellipsoid.
 * If dphi < rmin the query point is inside the EOA, instaed if dphi > rmax
 * the query point is certainly outside the EOA.
 * If rmin < dphi < rmax it is necessary to check using the definition of the 
 * EOA, which means || L * dphi || <= 1.
 * The output of the function is a boolean where true means that the query in 
 * inside the EOA.
 */
bool chemComp::inEOA(const VectorXd &phiQ) {
    
    VectorXd dphi = phiQ-getPhi();
    double dist = dphi.norm();
    
    if (dist > getRmax()) { 
        return false;
    } else if (dist < getRmin()) {
        return true;
    } else {  // rmin<||dphi||<rmax -> check with definition
        bool isInEOA = false;
        VectorXd check = getLT()*dphi;
         
        // compute module
        double module = check.norm();
        (module<=1.0) ? isInEOA = true : isInEOA = false;
        return isInEOA;
    }

}

/* 
 * canGrowEOA: check if the reaction mapping got from direct integration is 
 * inside the EOA
 * input:
 *      - phiQ: composition vector
 *      - RphiQ: mapping gradient
 * 
 * The error is defined as eps=||B*(R-Rl)||, where R - Rl = RphiQ - Rphi - A*(phiQ-phi)
 * It is evaluate the error and it is compared with the tollerance; if 
 * eps < epsTol the EOA can grow, instaed if eps > epsTol the EOA cannot grow.
 * The output of the function is a boolean where true means that the EOA can 
 * grow.
 */ 
bool chemComp::canGrowEOA(const VectorXd &phiQ, const VectorXd &RphiQ){
    bool grow = false;
    
    VectorXd epsGrow = getScaleMatrix()*(RphiQ-(getRphi()+getMapGrad()*(phiQ-getPhi())));
    
    // compute module
    double eps = epsGrow.norm();

    // control error<epsTol
    (eps < getEpsTol()) ? grow = true : grow = false;
    return grow;
}

/*
 * inizializeEOA: inizialize the EOA and evaluete the semiaxes of the 
 * input:
 *      - L: Cholesky lower matrix (output)
 *      - rmin: minor semiaxes of the ellipsoid (output)
 *      - rmax: max semiaxes of the ellipsoid (output)
 * The first step is to make the SVD of the mapping gradient, then if the 
 * singular values are smaller of 1/0.1*epsTol^2 they are substituted with 1/0.1*epsTol^2. 
 * After that the mapping gradient is approximated with the SVD of it and it is 
 * moltiplicated for the scale factor.
 * It is, now, used the QR decomposition decomposition to inizialize the matrix L which
 * represents the EOA; the singular values of the L matrix are evaluated in order 
 * to obtain the semiaxes of ellipsiod.
 */ 

void chemComp::inizializeEOA(const MatrixXd &A) {

	if (nSpec_< 100000000)
	{
		Eigen::JacobiSVD<MatrixXd> svd(nSpec_,nSpec_);
		MatrixXd Atilde(nSpec_,nSpec_);	
		
		// compute SVD decomposition
		svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		VectorXd sv = svd.singularValues();

		double svmin = getEpsTol()/(0.1*std::pow(getEpsTol(),0.5));
		double radius = 0.;
		// check singular value
		for(unsigned int i = 0; i<nSpec_; i++) {
			radius = getEpsTol()/std::max(sv(i),svmin);
			sv(i) = 1./radius;
		}

		MatrixXd matSv = sv.asDiagonal();
		Atilde = getScaleMatrix()*(svd.matrixU()*matSv*svd.matrixV().transpose());

		// QR decomposition
		// L^T = R 
		if (qrType_ == 0)
		{
			Eigen::FullPivHouseholderQR<MatrixXd> qr(Atilde.transpose());
			LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
		} 
		else if (qrType_ == 1)
		{
			Eigen::ColPivHouseholderQR<MatrixXd> qr(Atilde.transpose());
			LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
		} 
		else if (qrType_ == 2)
		{
			Eigen::HouseholderQR<MatrixXd> qr(Atilde.transpose());
			LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
		}
		// ensure that element of the diagonal are positive
		for(unsigned int i = 0; i<nSpec_; i++) {
			if(LT_(i,i) < 0) {
				LT_(i,i) = -LT_(i,i);
			}
		}
	  
		// compute max and min radius of EOA
		svd.compute(LT_);
		sv = svd.singularValues();
		rmin_ = 1./sv(0);
		rmax_ = 1./sv(nSpec_-1);
	}
	else
	{
		//Eigen::BDCSVD<MatrixXd> svd(nSpec_,nSpec_);
		MatrixXd Atilde(nSpec_,nSpec_);
		
		// disable floating point traps for BDCSVD
		int fpe = fegetexcept();
		fedisableexcept(FE_ALL_EXCEPT);

		// compute SVD decomposition
		//svd.compute(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::BDCSVD<MatrixXd> svd(A,Eigen::ComputeFullU|Eigen::ComputeFullV);

		// enable floating point traps after BDCSVD
	        feenableexcept(fpe);	
			
		VectorXd sv = svd.singularValues();
		
		double svmin = getEpsTol()/(0.1*std::pow(getEpsTol(),0.5));
		double radius = 0.;
		// check singular value
		for(unsigned int i = 0; i<nSpec_; i++) {
			radius = getEpsTol()/std::max(sv(i),svmin);
			sv(i) = 1./radius;
		}
		
		MatrixXd matSv = sv.asDiagonal();
		Atilde = getScaleMatrix()*(svd.matrixU()*matSv*svd.matrixV().transpose());

		// QR decomposition
		// L^T = R 
		if (qrType_ == 0)
		{
			Eigen::FullPivHouseholderQR<MatrixXd> qr(Atilde.transpose());
			LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
		} 
		else if (qrType_ == 1)
		{
			Eigen::ColPivHouseholderQR<MatrixXd> qr(Atilde.transpose());
			LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
		} 
		else if (qrType_ == 2)
		{
			Eigen::HouseholderQR<MatrixXd> qr(Atilde.transpose());
			LT_ = qr.matrixQR().triangularView<Eigen::Upper>();
		}
		// ensure that element of the diagonal are positive
		for(unsigned int i = 0; i<nSpec_; i++) {
			if(LT_(i,i) < 0) {
				LT_(i,i) = -LT_(i,i);
			}
		}
		
		// disable floating point traps for BDCSVD
                fedisableexcept(FE_ALL_EXCEPT);
			  	
		// compute max and min radius of EOA
		svd.compute(LT_);
		sv = svd.singularValues();
		
		// enable floating point traps after BDCSVD
		feenableexcept(fpe);
		
		rmin_ = 1./sv(0);
		rmax_ = 1./sv(nSpec_-1);
	}
     
}

chemComp::~chemComp() {
    phi_.resize(0);
    Rphi_.resize(0);
    LT_.resize(0,0);
    mapGrad.resize(0,0);
    B_.resize(0,0);
    scaleFactor_.resize(0);
}
