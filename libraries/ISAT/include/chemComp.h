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

#ifndef CHEM_COMP_H
#define CHEM_COMP_H

#if ISAT_USE_MKL == 1
    #define EIGEN_USE_MKL_ALL 
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include "Eigen/Core"
#include "Eigen/QR"
#include "Eigen/SVD"
#include <fenv.h>


using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;

class binaryNode;

class chemComp {
    private:
                // leaf index
        unsigned int index_;
                
        // number of species
        unsigned int nSpec_;

        // vector of composition
        VectorXd phi_;

        // reaction mapping
        VectorXd Rphi_;

        // mapping gradient
        MatrixXd mapGrad;

        // transpose of the Cholesky lower matrix of EOA
        MatrixXd LT_;

        // matrix scaleFactor
        MatrixXd B_;

        // vector scaleFactor
        VectorXd scaleFactor_;

        // min semiaxes of EOA
        double rmin_;

        // max semiaxes of EOA
        double rmax_;

        // epsTol
        double epsTol_;

        // parent node
        binaryNode *node_;

        // number of times used
        unsigned int nUsed_;

        // number of times grown
        unsigned int nGrown_;

        // last time used
        double lastTimeUse_;
        
        // typeQR
        unsigned int qrType_;

    public:
        // constructor 
        chemComp(const unsigned int index, const VectorXd &phi, const VectorXd &Rphi, const MatrixXd &A, const VectorXd &scaleFactor, const double epsTol, const unsigned int nSpec, const unsigned int qrType, binaryNode *nodeParent = NULL);
        
        // constructr using another leaf as input
        chemComp(chemComp &phi0);

        // grow the EOA
        void growEOA(const VectorXd &phiQ);

        // check if a point is in EOA
        bool inEOA(const VectorXd &phiQ);

        // check if the EOA can grow
        bool canGrowEOA(const VectorXd &phiQ, const VectorXd &RphiQ);

        virtual ~chemComp();

        //----------------Access Function-------------------------------
        // return the pointer to the chemComp node
        inline binaryNode *getNode() const {
            return node_;
        }

        // return phi
        inline VectorXd getPhi() const{
            return phi_;
        }

        // return Rphi
        inline VectorXd getRphi() const {
            return Rphi_;
        }

        // return mappingGradient
        inline MatrixXd getMapGrad() const {
            return mapGrad;
        }

        // return scaleFactor
        inline VectorXd getScaleFactor() const {
            return scaleFactor_;
        }

        // return epsTol
        inline double getEpsTol() const {
            return epsTol_;
        }

        // return nSpec
        inline unsigned int getSpec() const{
            return nSpec_;
        }

        // return nUsed
        unsigned int getUsed() const{
	    return nUsed_;
        }

        // return nGrown
        inline unsigned int getGrown() const {
            return nGrown_;
        }
                
                // return index
        inline unsigned int getIndex() const {
            return index_;
        }

        // return B
        inline MatrixXd getScaleMatrix() const {
            return B_;
        }

        // return L
        inline MatrixXd getLT() const {
            return LT_;
        }

        // return rmax
        inline double getRmax() const {
            return rmax_;
        }

        // return rmax
        inline double getRmin() const {
            return rmin_;
        }

        // return lastTimeUse
        inline double getLastTimeUse() const {
            return lastTimeUse_;
        }
        
        // return QRType
        inline double getQRType() const {
            return qrType_;
        }

        //----------------------- Set -----------------------------

        inline void setLastTimeUse(double lastTimeUse){
            lastTimeUse_ = lastTimeUse;
        }

        inline void setNode(binaryNode *node) {
            node_ = node;
        }

        // increment nUsed
        inline void nUsedPP() {
            nUsed_++;
        }

        // increment nGrown
        inline void nGrownPP() {
            nGrown_++;
        }
    //--------------------------------------------------------------
    private:
        // inizialize the EOA 
        void inizializeEOA(const MatrixXd &A);

};


#endif
