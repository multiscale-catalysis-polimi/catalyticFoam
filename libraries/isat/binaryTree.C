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

#include "binaryTree.h"
#include <iostream>

/*
 * binaryTree: costructor
 * Inizialize the tree and set to 0 the size of the tree
 */
binaryTree::binaryTree() {
        root_ = new binaryNode();
        size_ = 0;
        progIndex_ = 0;
}


/*
 * inserNewLeaf: a new leaf is added to the tree
 * input:
 *     - phi: vector of composition
 *     - Rphi: vector of reaction mapping
 *     - A: mapping gradient
 *     - scaleFactor: vectof of scale factor
 *     - epsTol: tollerance
 *     - nSpec: number of species
 *     - phi0: closest leaf (NULL if is not know), set to left of the node
 */

void binaryTree::insertNewLeaf(const VectorXd &phi, const VectorXd &Rphi, const MatrixXd &A, const VectorXd &scaleFactor, double epsTol,const unsigned int nSpec, const unsigned int qrType, chemComp *&phi0) {
  
    if(size_==0) {
        chemComp *data = new chemComp(progIndex_, phi,Rphi,A,scaleFactor,epsTol,nSpec,qrType,root_);
        root_->setElLeft(data);
        phi0 = data;
    } else {
    
        if(phi0 == NULL) {
            chemComp *near;
            searchTreeLeaf(phi,root_,near);
            phi0=dynamic_cast<chemComp*>(near);
        }
              
        binaryNode *parentNode = phi0->getNode();
        chemComp *data = new chemComp(progIndex_,phi,Rphi,A,scaleFactor,epsTol,nSpec,qrType);
        binaryNode *newNode;
                       
        if(size_ > 1) {  
            // initialize a new node
            newNode = new binaryNode(phi0, data, parentNode);
            // insert the node in the tree
            insertNode(newNode, phi0);
        } else { // size == 1
            // initialize a new node which will have two leaf
                        deleteData(root_);
            newNode = new binaryNode(phi0, data, NULL);
            root_ = newNode;
        }
        phi0->setNode(newNode);
        data->setNode(newNode);
    }
    size_++;
    progIndex_++;
}

 /*
 * insertNode: substitute one leaf with a node and set the pointer 
 * to the left and right son leaf
 * input:
 *     - newNode: new node
 *     - compInz: initial composition
 */
void binaryTree::insertNode(binaryNode *&newNode, chemComp *&compInz) { 
    if(compInz==compInz->getNode()->getElRight()) {
        compInz->getNode()->setRight(newNode);
        compInz->getNode()->setElRight(NULL);
    } else {
        compInz->getNode()->setLeft(newNode);
        compInz->getNode()->setElLeft(NULL);
    }    
}

/*
 * searchTreeLeaf: a new leaf is added to the tree
 * input:
 *     - query: vector of composition query
 *     - node: tree/subtree of searching
 *     - nearestComp: leaf near query point (output)
 * 
 * if v^T * query > a go on the rigth, else go on the left
 */
void binaryTree::searchTreeLeaf(const VectorXd &query, binaryNode *node, chemComp *&nearestComp) {

    if(size_ > 1) {
        double conf = 0.0;

        conf = query.dot(node->v());    
        
        if(conf > node->a()) { // go to the rigth
            if(node->getRight() != NULL) { // continue search in rigth subtree
                searchTreeLeaf(query, node->getRight(), nearestComp);
            } else { 
                nearestComp=node->getElRight();
            }
        } else {
            if(node->getLeft() != NULL) { // continue search in left subtree
                searchTreeLeaf(query, node->getLeft(), nearestComp);
            } else { // sei su una foglia prendi elemento 
                nearestComp=node->getElLeft();
            }
        }
    } else if(size_ == 1) { 
        nearestComp = node->getElLeft();
    } else { 
        nearestComp = NULL;
    }    
    
}

/*
 * deleteLeaf: remove a leaf from the tree
 * input:
 *     - phi0: leaf
 * 
 * delete a leaf and reshape the tree
 */
void binaryTree::deleteLeaf(chemComp *&phi0)
{

    if(size_ == 1) //only one point is stored
    {
        deleteData(phi0);
        deleteData(root_);
    }
    else if (size_ > 1)
    {
        binaryNode* z = phi0->getNode();
        binaryNode* x;
        chemComp* siblingPhi0 = chemCompSibling(phi0);
              
        if (siblingPhi0 != NULL)//the sibling of phi0 is a leaf
        {
            if (z->getParent() == NULL) //z was root (only two leaf in the tree)
            {
                root_ = new binaryNode();
                root_->setElLeft(siblingPhi0);
                siblingPhi0->setNode(root_);
                
            }
            else if (z==z->getParent()->getLeft())
            {
                z->getParent()->setElLeft(siblingPhi0);
                z->getParent()->setLeft(NULL);
                siblingPhi0->setNode(z->getParent());
            }
            else
            {
                z->getParent()->setElRight(siblingPhi0);
                z->getParent()->setRight(NULL);
                siblingPhi0->setNode(z->getParent());
            }
        }
        else
        {
            x = nodeSibling(phi0);
            transplant(z,x);
        }
        
        deleteData(phi0);
        deleteData(z);
    }
    size_--;
}

/*
 * transplant
 * input:
 *     - u: binary node
 *     - v: binary node
 * 
 */
void binaryTree::transplant(binaryNode* u, binaryNode* v)
{
    if (u->getParent() == NULL)
    {
        root_ = v;
    }
    else if (u == u->getParent()->getLeft())
    {
        u->getParent()->setLeft(v);
    }
    else
    {
        u->getParent()->setRight(v);
    }
    
    if (v != NULL)
    {
        v->setParent(u->getParent());
    }

}

/*
 * depth: evaluate the depth of a subtree
 * input:
 *     - subTree: binary node
 *     
 */
int binaryTree::depth(binaryNode *subTree) {
    if(subTree != NULL)
        return 1+std::max(depth(subTree->getLeft()),depth(subTree->getRight()));
    else
        return 0;
}

/*
 * clear: clear the whole tree
 * input:
 *     - NO INPUT
 *     
 */
void binaryTree::clear() { 
    deleteSubTree(root_);
    root_ = new binaryNode;
    size_ = 0;
    progIndex_ = 0;
}

/*
 * deleteSubTree: delete a certain subTree
 * input:
 *     - subTree: binary node
 *     
 */
void binaryTree::deleteSubTree(binaryNode* subTree) {
    if(subTree != NULL) {
        deleteData(subTree->getElLeft());
        deleteData(subTree->getElRight());
        deleteSubTree(subTree->getLeft());
        deleteSubTree(subTree->getRight());
        deleteData(subTree);
    }
}

/*
 * *treeMin: return pointer to the most left leaf
 * input:
 *      - node: tree/subtree of searching
 */
chemComp* binaryTree::treeMin(binaryNode *subTree) { 
    if(subTree !=NULL) {
            while(subTree->getLeft() != NULL) {
                subTree = subTree->getLeft();
            }
            return subTree->getElLeft();
        } else {
            return NULL;
        }
}

/*
 * treeNextLeaf: return the next leaf to a given one
 * input:
 *      - leaf: chemComp (tree leaf)
 */
chemComp* binaryTree::treeNextLeaf(chemComp* leaf) {
    if(size_ > 1) {
        // if leaf is the left -> look in the rigth side
        if(leaf == leaf->getNode()->getElLeft()) {
            binaryNode* pn = leaf->getNode();
            if(pn->getRight() == NULL) { // rigth side is only rigth element
                return pn->getElRight();
            } else { // rigth side is a subtree -> return the most left element of
                //this sub tree
                return treeMin(pn->getRight());
            }
        } else { // leaf is the rigth -> climb up the tree
            binaryNode *pny = leaf->getNode();
            while(pny->getParent() != NULL) {
                // look if the parent of the leaf is on the left
                if(pny == pny->getParent()->getLeft()) {
                    if(pny->getParent()->getRight() == NULL) {
                        return pny->getParent()->getElRight();
                    } else {
                        return treeMin(pny->getParent()->getRight());
                    }
                }
                pny = pny->getParent();
                
            }
            // if reach this point it means that the tree has been climbed up
            // until the root from the rigth side so there isn't next
            return NULL;
        }
    } else {
        return NULL;
    }
}

/*
 * balance: balance routine (by F. Contino)
 * input:
 *      - NO INPUT
 */
void binaryTree::balance() 
{
    if(size_ > 1) {
        // evaluation of the leaf mean value
        std::vector<chemComp*> leafList;
        leafList.reserve(size_); // memory allocation
        chemComp *x = getTreeMin();
        int nSpec_ = x->getSpec();
        VectorXd mean(nSpec_);
        VectorXd variance(nSpec_);

	for(unsigned int i=0; i<nSpec_; i++)
	{
		mean(i)=0.;
		variance(i)=0.;
	}

        while(x != NULL) {
            mean += x->getPhi();
            leafList.push_back(x);
            x = treeNextLeaf(x);
        }

        mean /= size_;

        // evaluation of the variance in each direction
        for(std::vector<chemComp*>::iterator iter = leafList.begin(), end=leafList.end(); iter != end; ++iter) {
            for(int i = 0; i < nSpec_; i++)
                variance(i) += std::pow((*iter)->getPhi()(i)-mean(i),2);
        }
        std::sort(variance.data(),variance.data()+variance.size());
        
	// find the maximun variance direction
        int maxDir = -1;
        int hmTested = 1;
        double lfLeft = 0.;
        double bestBalance = double(size_);

        while(
                ((lfLeft < 0.35*double(size_)) || (lfLeft > 0.65*double(size_))) && (hmTested < variance.rows()-1)
             ) 
        {
            lfLeft = 0.;

            int curDir = variance.rows()-hmTested;
            
            for(std::vector<chemComp*>::iterator i = leafList.begin(), end = leafList.end(); i != end; ++i) {
                if((*i)->getPhi()(curDir) < mean(curDir)) {
                    lfLeft++;
                }
            }

            if(fabs(lfLeft-double(size_)*0.5) < bestBalance) {
                maxDir = curDir;
                bestBalance = fabs(lfLeft-double(size_)*0.5);
            }

            hmTested++;
        }

        std::vector<chemComp*> tmp;
        tmp.reserve(size_);
        srand(time(NULL));
        
	for(unsigned int k = 0; k < size_; k++) {
            // random generation number
            unsigned int num = (rand() % (leafList.size()));
            unsigned int cont = 0;
            for(std::vector<chemComp*>::iterator i = leafList.begin(), end = leafList.end(); i != end; ++i) {
                if(cont == num) {
                    tmp.push_back(*i);
                    leafList.erase(i);
                    break;
                }
                cont++;
            }
        }
        leafList = tmp;
        tmp.clear();
        tmp.resize(0);

        // find most left and most rigth leaf in respect to the selected direction
        chemComp *min = NULL;
        chemComp *max = NULL;
        double maxPhi = 0.0;
        double minPhi = 1.e32;
        
	for(std::vector<chemComp*>::iterator iterMM = leafList.begin(), end = leafList.end(); iterMM != end; ++iterMM) {
            double phiDir = (*iterMM)->getPhi()(maxDir);
            if(phiDir>maxPhi)
            {
                max = *iterMM;
                maxPhi=phiDir;
            }
            if(phiDir<minPhi)
            {
                min = *iterMM;
                minPhi = phiDir;
            }
        }
       
 	
        // cancel the node of the tree
        deleteAllNode(root_);
        root_ = NULL;

        //creo un nodo per min e max
        binaryNode *newNode = new binaryNode(min, max, NULL);             
        root_ = newNode;

        min->setNode(newNode);
        max->setNode(newNode);

        // rebuild the tree
        
        for(std::vector<chemComp*>::iterator it = leafList.begin(), end = leafList.end(); it != end; ++it) {
            if( *it != min && *it != max) {
                //search tree for position
                chemComp* phi0;
                searchTreeLeaf((*it)->getPhi(),root_,phi0);
                chemComp* leaf = dynamic_cast<chemComp*>(phi0);
                //add the chemPoint
                binaryNode* nd = new binaryNode(leaf,*it, leaf->getNode());
                insertNode(nd, leaf);//make the parent of phi0 point to the newly created node
                leaf->setNode(nd);
                (*it)->setNode(nd);
            }
        }
        leafList.clear();
        leafList.resize(0);
	
	
    }
    return;
}


/*
 * deleteAllNode: delete all the nodes of the tree
 * input:
 *      - x: binary node
 */
void binaryTree::deleteAllNode(binaryNode* x) {
    if(x != NULL) {
        deleteAllNode(x->getLeft());
        deleteAllNode(x->getRight());
        deleteData(x);
    }
}
    
/*
 * chemCompSibling: given a chemComp phi0 retrieve, if exist, the other chemComp of the same root
 * input:
 *      - x: leaf
 */
chemComp* binaryTree::chemCompSibling(chemComp* x) {
    
    if(size_>1)
    {
        binaryNode* y = x->getNode();
        if(x==y->getElLeft()) //x is on the left, return right side
        {
            return y->getElRight();
        }
        else//x is on the right, return left side
        {
            return y->getElLeft();
        }
    }
    else
    {
        return NULL;
    }
}

/*
 * chemCompSibling: given a chemComp return the other side of the node parent
 * input:
 *      - x: leaf
 */
binaryNode* binaryTree::nodeSibling(chemComp* x) {
    
    if(size_>1)
    {
        binaryNode* y = x->getNode();
        if(x==y->getElLeft()) //x is on the left, return right side
        {
            return y->getRight();
        }
        else//x is on the right, return left side
        {
            return y->getLeft();
        }
    }
    else
    {
        return NULL;
    }
}

/*
 * countNode: count the number of node of a certain subTree
 * input:
 *      - tree: node
 */
int binaryTree::countNode(binaryNode *tree) { 
        if(tree == NULL) {
        return 0;
    }    
    else {
        int count = 1;
        count += countNode(tree->getLeft());
        count += countNode(tree->getRight());
        return count;
    }
}    

/*
 * ~binaryTree: destructor
 */
binaryTree::~binaryTree() {
    delete root_;
}

