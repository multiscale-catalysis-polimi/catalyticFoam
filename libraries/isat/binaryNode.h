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

#ifndef BYNARY_NODE_H
#define BYNARY_NODE_H

#if ISAT_USE_MKL == 1
    #define EIGEN_USE_MKL_ALL 
#endif

#include <stdlib.h>
#include <stdio.h>
#include "chemComp.h"
#include "Eigen/Core"

using namespace std;

using Eigen::VectorXd;

class binaryNode {
    private:
        // left and rigth leaf
        chemComp *dLeft_;
        chemComp *dRight_;
        
        // pointer to parent node and to children node
        binaryNode *left_, *right_, *parent_;

        // navigate elements
        VectorXd v_;
        double a_;

        // number of specie
                int nSpec;

    public:    
                // costructor empty node
        binaryNode();
            
                // constructor of a new node with two leaf
        binaryNode(chemComp *compInz, chemComp *compNew, binaryNode *parentN);
       
        //----------------Access Function-------------------------------
                // return the pointer to the left node
                inline binaryNode *&getLeft()  {
                    return left_;
                }
                
                // return the pointer to the rigth node
                inline binaryNode *&getRight() {
                    return right_;
                }
                
                // return the pointer to the parent node
                inline binaryNode *&getParent() {
                    return parent_;
                }
                
                // return the pointer to the left leaf
                inline chemComp *&getElLeft() {
                    return dLeft_;
                }
                
                // return the pointer to the rigth leaf
                inline chemComp *&getElRight() {
                    return dRight_;
                }

         // return vector v
                inline VectorXd v() const{
                    return v_;
                }
                
                // return scalar a
                inline double a() const{
                    return a_;
                }
                
        //----------------------- Set -----------------------------
        inline void setLeft(binaryNode *left){
                    left_ = left;
                }

        inline void setRight(binaryNode *right){
                    right_ = right;
                }

        inline void setParent(binaryNode *parent){
                    parent_ = parent;
                }

        inline void setElLeft(chemComp *elLeft){
                    dLeft_ = elLeft;
                }

        inline void setElRight(chemComp *elRight){
                    dRight_ = elRight;
                }

        // destructor
                virtual ~binaryNode();
    private:
        // calculate vector v and scalar a
        void calc(chemComp *left, chemComp *right, VectorXd &v, double &a);
};

#include "binaryNode.C"

#endif
