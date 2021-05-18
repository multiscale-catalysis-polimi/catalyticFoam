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
|   binaryNode class:                                                     |
|                                                                         |
|   Binary Tree node instance                                             |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "binaryNode.h"

/*
 * binaryNode: NULL costructor
 */
binaryNode::binaryNode() { 
        left_=NULL;
        right_=NULL;
        parent_=NULL;
        dLeft_=NULL;
        dRight_=NULL;
}

/*
 * binaryNode: costructor
 */        
binaryNode::binaryNode(chemComp *compInz, chemComp *compNew, binaryNode *parentN) {
            
        left_=NULL;
        right_=NULL;
        parent_=parentN;
        nSpec=compInz->getSpec();
        v_.resize(nSpec);
        calc(compInz, compNew, v_, a_);
        dLeft_=compInz;
        dRight_=compNew;
}

// destructor
binaryNode::~binaryNode() {
}

/*
 * calc: evaluate the vector v and the double a which separates the two ellipsods
 * input:
 *      - left: left leaf
 *      - rigth: rigth leaf
 *      - v: hyperplane vector (output)
 *      - a: scalar (output)
 * 
 * from "Algorithm for ellipsoids" - Pope
 * 
 *     L1*L1^T(phiQ-phi0)     
 * v=-------------------- , phih = (phi0+phiQ)/2 ; phi0 left, phiQ rigth, L1 left leaf
 *    ||L1*L1^T(phiQ-phi0)||
 */

void binaryNode::calc(chemComp *left, chemComp *right, VectorXd &v, double &a) {

    
    VectorXd phih = 0.5*(left->getPhi() + right->getPhi());
    MatrixXd L1 = left->getLT();

    v = L1*L1.transpose()*(right->getPhi() - left->getPhi());

    double norm = v.norm();

    v /= norm;
    a = v.dot(phih);

    return;
}

