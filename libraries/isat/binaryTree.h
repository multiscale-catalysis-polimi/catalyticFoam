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
|   binaryTree class:                                                     |
|                                                                         |
|   Binary Tree creation and management                                   |
|   Binary Tree searching and adding procedure                            |
|   Binary Tree balance procedure                                         |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef BINARY_TREE_H
#define    BINARY_TREE_H

#if ISAT_USE_MKL == 1
    #define EIGEN_USE_MKL_ALL 
#endif

#include <stdio.h>
#include <vector>
#include "binaryNode.h"
#include "Eigen/Core"

using namespace std;

//class binaryNode;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class binaryTree {
    
    private:
        // pointer to the tree root
        binaryNode *root_;
        
        // size of the tree (number of leaf)
        unsigned int size_;
        
        // leaf progressive index
        unsigned int progIndex_;
        
    public:
        // costructor
        binaryTree();
        
        // insert new leaf in the tree
        void insertNewLeaf(const VectorXd &phi, const VectorXd &Rphi, const MatrixXd &A, const VectorXd &scaleFactor, double epsTol, const unsigned int nSpec, const unsigned int qrType, chemComp *&phi0);
        
	// search the nearest leaf starting from a composition
        void searchTreeLeaf(const VectorXd &query, binaryNode *node, chemComp *&nearestComp);
               
        // delete a leaf from the tree and reshape the tree
        void deleteLeaf(chemComp *&leaf);
        
        void balance();
        
        // clear the tree
        void clear();
        
        // most left leaf of the subtree
        chemComp *treeMin(binaryNode *subTree);
        
        // return the next leaf starting from one known
        chemComp *treeNextLeaf(chemComp *leaf);
        
        // 
        int getSizeCrossingTree() {
            chemComp *leaf = getTreeMin();
            int cont = 0;
            while ( leaf != NULL) {
                cont++;
                leaf=treeNextLeaf(leaf);
            }
            return cont;
        }
        
        // destructor
        virtual ~binaryTree();
        
        
        // access to class element
        // return the pointer to the tree root
        inline binaryNode *getRoot() const {
            return root_;
        }
        
        // return the size of the tree
        inline unsigned int getSize() const {
            return size_;
        }
        
        // return the most left leaf of the tree     
        inline chemComp *getTreeMin()
        {
            return treeMin(root_);
        }

        // return the difference between the depth of the left side and the depth of the rigth side
        inline int getHeight()
        {
            return fabs(depth(root_->getLeft())-depth(root_->getRight()));
        }

        // return the depth of the tree
        inline int getDepth()
        {
            return depth(root_);
        }
        
        private:
            // insert a new node
            void insertNode(binaryNode *&newNode, chemComp *&compInz);
            // remove a subtree
            void deleteSubTree(binaryNode *subTree);      
            // perform transplant 
            void transplant(binaryNode* u, binaryNode* v);
            // remove all the node of the BT
            void deleteAllNode(binaryNode* x);
            // evaluate depth of BT
            int depth(binaryNode *subTree);    
            // find the nearest leaf to a certain leaf
            chemComp* chemCompSibling(chemComp* x);
            // find the nearest node (not parent) to a certain leaf
            binaryNode* nodeSibling(chemComp* x);
            // count number of node
            int countNode(binaryNode *tree);
            // delete pointer to leaf
            template<typename T>
            inline void deleteData(T &dt) {
                if(dt) {
                    delete dt;
                    dt = NULL;
                }
            }
            
};

#endif    /* BINARYTREE_H */

