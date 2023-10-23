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
|   ISAT class:                                                           |
|                                                                         |
|   ISAT algorithm operation (retrieve, growth, addition)                 |
|   linked list search (MRU, MFU)                                         |
|   cleaning and balance procedure                                        |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "ISAT.h"
#include <iostream>

/*
 * ISAT: constructor
 */
ISAT::ISAT(VectorXd scaleFactor, double epsTol, unsigned int nCol) 
{
    // create a pointer to the three
    chemistryTree_ = new binaryTree();
    
    // inizialize parameters
    scaleFactor_ = scaleFactor;
    epsTol_ = epsTol;
    nSpec_ = nCol;   
    
    // set max dimension of MRU e BT
    maxSizeMRU_ = 10;
    maxSizeMFU_ = 30;
    maxSearchMRU_ = 10;
    maxSearchMFU_ = 30;
    maxSizeBT_ = 50000;

    // set param BTSR;
    maxLevelBTSR_ = 5;
    
    // flag :
    //     - search: BTSR, MRU, MFU, BruteForce
    //     - clear if full
    flagSearchBTSR_ = true; //true;
    flagSearchMRU_ = true;//true;
    flagSearchMFU_ = true;//true;
    flagSearchBruteForce_ = false; //false;
    flagSearchCL_ = true; //true;
    flagClearingIfFull_ = false;
    flagCleanAndBalance_ = true;
    
    // balance param
    balanceFactorAdd_ = 0.01;
    balanceFactorRet_ = 2;
    maxHeightCoeff_ = 5.;
    maxGrowCoeff_ = 1.;
    minUsedCoeff_ = 0.005;
    maxTimeOldCoeff_ = 0.7;
    
    // initialize counter
    nAdd_=0;
    nGrow_=0;
    nUse_=0;
    nBTS_=0;
    nBTSR_=0;
    nMRU_=0;
    nMFU_=0;
    nCL_=0;
    nFailBTinEOA_=0;
    
    // balance counter
    nBalance_ = 0;
    nAddFromLastBalance_ = 0;
    nUseFromLastBalance_ = 0;
    nRemovedLeaves_ = 0;
    
    // initialize MRU, MFU and toRemove list
    mfu_.clear();
    mru_.clear();
    toRemove_.clear();
    
    // memorize start time
    timeStart_ = double(std::clock())/CLOCKS_PER_SEC;
    // inizialize time counter
    tStartOp_ = 0.;
    tEndOp_ = 0.;
    cpuRet_ = 0.;
    cpuInt_ = 0.;
    cpuGrw_ = 0.;
    cpuAdd_ = 0.;

    // QR type
    qrType_ = 0;
}

/*
 * retrieve: find the closest leaf to a given query point
 * input:
 *      - phiQ: vector of composition query
 *      - closest: leaf closest to phiQ (output)
 * 
 * The function call binaryTree::searchTreeLeaf in order to find the closest 
 * leaf to a given query point.
 * The output of the function is a boolean where true means that it was found 
 * the leaf and that leaf is saved in chemComp *closest.
 * If the search fails, the MRU search and the MFU search are performed and if
 * on also the BruteForce
 */
bool ISAT::retrieve(const VectorXd &phiQ, chemComp *&closest)
{
    // BinaryTree Search
    tStartOp_ = double(std::clock())/CLOCKS_PER_SEC;

    chemistryTree_->searchTreeLeaf(phiQ, chemistryTree_->getRoot(), closest);

    if(!closest) 
    { 
        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
        cpuRet_ += tEndOp_-tStartOp_;
        return false; 
    }
    else 
    {
        chemComp *phi0 = dynamic_cast<chemComp*>(closest);
        
	// Binary Tree Search
        if(phi0->inEOA(phiQ)) 
        {
            // retrieve leaf -> add to MRU list
            nBTS_++;
            phi0->nUsedPP();
            phi0->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);
            
            addToMRU(phi0);
            addToMFU(phi0);

            tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
            cpuRet_ += tEndOp_-tStartOp_;
            return true;
        }
	chemComp* prova = chemistryTree_->getTreeMin();
	Eigen::VectorXd testp = prova->getPhi()-phiQ;
	double normp = testp.norm();
	double normpmin = normp;

	chemComp* minLeaf = prova;

	//Reverse Binary Tree Search
	if(flagSearchBTSR_)
	{
		binaryNode* currentNode = phi0->getNode();
		binaryNode* previousNode = NULL;
		unsigned int level=0;
		bool foundLeaf = false;

		while(currentNode!=NULL && level<=maxLevelBTSR_ && !foundLeaf)
		{
			if(level == 0)
			{
				if(currentNode->getElRight() == phi0)
				{
					chemComp* testLeaf = currentNode->getElLeft();

					if(testLeaf!=NULL && testLeaf->inEOA(phiQ))
					{
						nBTSR_++;
						closest = testLeaf;
						testLeaf->nUsedPP();
            					testLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

						addToMRU(testLeaf);
            					addToMFU(testLeaf);

						tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
            					cpuRet_ += tEndOp_-tStartOp_;
					
						foundLeaf = true;
					}
				}
				else
				{
					chemComp* testLeaf = currentNode->getElRight();

                                        if(testLeaf!=NULL && testLeaf->inEOA(phiQ))
                                        {
                                                closest = testLeaf;

                                                nBTSR_++;
                                                testLeaf->nUsedPP();
                                                testLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

                                                addToMRU(testLeaf);
                                                addToMFU(testLeaf);

                                                tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                                                cpuRet_ += tEndOp_-tStartOp_;

                                                foundLeaf = true;
                                        }
				}
			}
			else
			{
				if(currentNode->getRight() == previousNode)
				{
					chemComp* testLeaf = currentNode->getElLeft();

					if(testLeaf!=NULL && testLeaf->inEOA(phiQ))
                                        {
                                                closest = testLeaf;

                                                nBTSR_++;
                                                testLeaf->nUsedPP();
                                                testLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

                                                addToMRU(testLeaf);
                                                addToMFU(testLeaf);

                                                tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                                                cpuRet_ += tEndOp_-tStartOp_;

                                                foundLeaf = true;
                                        }

					if(!foundLeaf)
					{
						testLeaf =currentNode->getElRight();
						
						if(testLeaf!=NULL && testLeaf->inEOA(phiQ))
						{
							closest=testLeaf;

							//nBTSR_++;
							testLeaf->nUsedPP();
							testLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

							addToMRU(testLeaf);
							addToMFU(testLeaf);

							tEndOp_=double(std::clock())/CLOCKS_PER_SEC;
							cpuRet_+=tEndOp_-tStartOp_;

							foundLeaf = true;
						}
					}
				}
				else
				{
					chemComp* testLeaf = currentNode->getElRight();

					if(testLeaf!=NULL && testLeaf->inEOA(phiQ))
                                        {
                                                closest = testLeaf;

                                                nBTSR_++;
                                                testLeaf->nUsedPP();
                                                testLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

                                                addToMRU(testLeaf);
                                                addToMFU(testLeaf);

                                                tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                                                cpuRet_ += tEndOp_-tStartOp_;

                                                foundLeaf = true;
                                        }

					if(!foundLeaf)
					{
						testLeaf=currentNode->getElLeft();

						if(testLeaf!=NULL && testLeaf->inEOA(phiQ))
						{
							closest=testLeaf;

							nBTSR_++;
							testLeaf->nUsedPP();
							testLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

							addToMRU(testLeaf);
							addToMFU(testLeaf);

							tEndOp_=double(std::clock())/CLOCKS_PER_SEC;
							cpuRet_+=tEndOp_-tStartOp_;

							foundLeaf=true;
						}
					}

				}
			}

			previousNode = currentNode;
			currentNode = currentNode->getParent();

			level++;
		}

		if(foundLeaf)
			return true;	
	
	}	
        
        // MRU search
        std::list<chemComp*> tempMFU = mfu_;
        if(flagSearchMRU_)
        { // MRU search
            unsigned int cont = 0;
            std::list<chemComp*>::iterator iter=mru_.begin();
            while(iter != mru_.end() &&  cont <= maxSearchMRU_) 
            {
                if((*iter)->getIndex() != phi0->getIndex()) 
                { // skip previous leaf  
                    if((*iter)->inEOA(phiQ)) 
                    { // if the leaf cover the query return
                        closest = (*iter);
                        (*iter)->nUsedPP();
                        (*iter)->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);
                        
                        // add to MRU and MFU
                        addToMFU((*iter));
                        addToMRU((*iter));
                        nMRU_++;
                        
                        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                        cpuRet_ += tEndOp_-tStartOp_;
                        tempMFU.clear();
                        return true;   
                    }
                }
                tempMFU.remove((*iter)); // remove from temporary MFU the leaf already used
                cont++;
                ++iter;
            }  
        }
 
        // MFU search
        if(flagSearchMFU_)
        { 
            unsigned int cont = 0;
            std::list<chemComp*>::iterator iter=tempMFU.begin();
            while(iter != tempMFU.end() &&  cont <= maxSearchMRU_) 
            {
                if((*iter)->getIndex() != phi0->getIndex()) // skip previous leaf 
                {  
                    if((*iter)->inEOA(phiQ)) // if the leaf cover the query return
                    { 
                    
                    closest = (*iter);
                    (*iter)->nUsedPP();
                    (*iter)->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);
                    
                    // add to MRU and MFU
                    addToMFU((*iter));
                    addToMRU((*iter));
                    nMFU_++;
                    
                    tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                    cpuRet_ += tEndOp_-tStartOp_;
                    tempMFU.clear();
                    return true;   
                    }
                }
                cont++;
                ++iter;
            }  
        }

	// Closest Leaf Search        
	if(flagSearchCL_)
	{
		chemComp* prova = chemistryTree_->getTreeMin();
        	Eigen::VectorXd testp = prova->getPhi()-phiQ;
        	double normp = testp.norm();
        	double normpmin = normp;
	
        	chemComp* minLeaf = prova;

		while(prova!=NULL)
        	{
                	prova = chemistryTree_->treeNextLeaf(prova);
                	if(prova!=NULL)
                	{
                        	testp = prova->getPhi()-phiQ;
                        	normp = testp.norm();

                        	if(normp<normpmin)
                        	{
                                	minLeaf=prova;
                                	normpmin = normp;
                        	}
                	}
		}

		if(minLeaf!=phi0 && minLeaf->inEOA(phiQ))
		{
			closest=minLeaf;
			nCL_++;

			minLeaf->nUsedPP();
			minLeaf->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);

			addToMRU(minLeaf);
			addToMFU(minLeaf);

			tEndOp_=double(std::clock())/CLOCKS_PER_SEC;
			cpuRet_+=tEndOp_-tStartOp_;

			return true;
		}	
	}

        // BruteForce search
        if(flagSearchBruteForce_ && (chemistryTree_->getSize() < 2500) )
        { 
            phi0 = chemistryTree_->getTreeMin();
            while(phi0 != NULL) 
            {
                if(phi0->inEOA(phiQ)) 
                {
                    nFailBTinEOA_++;
                    phi0->nUsedPP();
                    phi0->setLastTimeUse(double(std::clock())/CLOCKS_PER_SEC);
                    
                    // add to MRU and MFU
                    addToMRU(phi0);
                    addToMFU(phi0);

                    tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                    cpuRet_ += tEndOp_-tStartOp_;
                    return true;
                }
                phi0 = chemistryTree_->treeNextLeaf(phi0);               
            }
            // if Brute Force failed doesn't exist any leaf so return false
            tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
            cpuRet_ += tEndOp_-tStartOp_;
            return false;
        } 
        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
        cpuRet_ += tEndOp_-tStartOp_;
        return false;
    }
    
}

/*
 * add: add a new leaf in the ISAT tree
 * input:
 *      - phi: vector of composition
 *      - Rphi: reaction mapping 
 *      - A: mapping gradient
 *      - phi0: closest leaf
 * 
 *     First, check if the tree is full. If it is check which full table treatment
 *  and perform it adding at the end the last leaf. If it is not full add the leaf
 *     
 */
bool ISAT::add(const VectorXd &phi, const VectorXd &Rphi, const MatrixXd &A, chemComp *&phi0) {

    tStartOp_ = double(std::clock())/CLOCKS_PER_SEC;
    
    unsigned int checkSize = chemistryTree_->getSize(); // tree size before add operation
      
    if(checkSize >= maxSizeBT_)  // full table
    {
        if(mru_.size() > 0) 
        {
            if(flagClearingIfFull_) // clear the whole tree and add the MRU and MFU
            {
		std::list<chemComp*> tempLeaves;
                std::list<chemComp*> tempMFU = mfu_;
		std::vector<std::pair<bool,bool>> isInLists;
                std::list<chemComp*>::iterator it;
                for(it=mru_.begin(); it != mru_.end(); ++it) 
		{
			std::pair<bool,bool> coupleLists;
                        tempLeaves.push_front(new chemComp(**it));
			coupleLists.first = true;
			unsigned int sizetempmfu = tempMFU.size();
                        tempMFU.remove(*it);
			if(sizetempmfu-tempMFU.size()!=0)
			{
				coupleLists.second = true;
			}
			else
			{
				coupleLists.second = false;
			}
			isInLists.push_back(coupleLists);
                }

                // aggiungo le foglie della MFU che non erano presenti nella MRU
                for(it=tempMFU.begin(); it != tempMFU.end(); ++it) 
		{
                        tempLeaves.push_front(new chemComp(**it));
			std::pair<bool,bool> coupleLists;
			coupleLists.first = false;
			coupleLists.second = true;
			isInLists.push_back(coupleLists);
                }

                // pulisci albero
                chemistryTree_->clear();
                nAddFromLastBalance_ = 0;

                // Clear the lists
                mru_.clear();
                mfu_.clear();
		toRemove_.clear();

                // Insert new point
                chemComp *nPhi = NULL;
                chemistryTree_->insertNewLeaf(phi,Rphi,A,scaleFactor_,epsTol_,nSpec_,qrType_,nPhi);

                // memorizzo questa foglia in modo da aggiungerla alla MRU in testa dopo averla ripopolata
		chemComp *phiTemp = chemistryTree_->getTreeMin();
                
		// aggiungo gli elementi delle tempMRU
                std::list<chemComp*>::iterator iter;

		unsigned int itBoolList = isInLists.size()-1;

                // faccio scorrere la lista, se trovo la foglia ottenuta dalla BT la salto
                for(iter=tempLeaves.begin(); iter != tempLeaves.end(); ++iter) 
		{
		    chemComp *newPhi = NULL;

                    chemistryTree_->insertNewLeaf((*iter)->getPhi(),(*iter)->getRphi(),(*iter)->getMapGrad(),scaleFactor_,epsTol_,nSpec_, qrType_,newPhi);

                    // inserisco (se non NULL) nella MRU e nella MFU l'ultima foglia inserita nell'albero
                    if(newPhi->getNode()->getElRight()) 
		    {	     
			     if(isInLists[itBoolList].first)
                             	mru_.push_front(newPhi->getNode()->getElRight());

			     if(isInLists[itBoolList].second)
                             	mfu_.push_front(newPhi->getNode()->getElRight());
                    }

                    delete *iter;
		    itBoolList--;
                }

                // inserisci primo elemento della lista (ovvero il nuovo punto appena aggiunto)
                if(mru_.size()<maxSizeMRU_)
                	mru_.push_front(phiTemp);

		if(mfu_.size()<maxSizeMFU_)
                	mfu_.push_back(phiTemp);

                tempMFU.clear();
                tempLeaves.clear();
		nAdd_++;

                return true;
            } 
            else 
            { // delete the leaf at the end of the MFU and replace it using this new leaf  

                unsigned int toDel = mfu_.back()->getIndex();
                bool sameLeaf = false;
                if( toDel == phi0->getIndex()) {
                    sameLeaf = true;
                    // remove dalla toRemove list
                    toRemove_.remove(phi0);
                }

                cleanLeaf(mfu_.back());

                chemComp *nPhi = NULL;
                if(!sameLeaf)
                    nPhi = dynamic_cast<chemComp*>(phi0);
                else
                    chemistryTree_->searchTreeLeaf(phi, chemistryTree_->getRoot(), nPhi);
               
		//Blocco Test 4
		//test4_ -= double(std::std::clock())/CLOCKS_PER_SEC;
                chemistryTree_->insertNewLeaf(phi,Rphi,A,scaleFactor_,epsTol_,nSpec_,qrType_,nPhi);
                phi0 = nPhi;
		//test4_ += double(std::std::clock())/CLOCKS_PER_SEC;
		//Fine Blocco Test 4

                if(phi0->getNode()->getElRight()) {
                        addToMRU(phi0->getNode()->getElRight());
                        addToMFU(phi0->getNode()->getElRight());
                }

                nAdd_++;
                tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
                cpuAdd_ += tEndOp_-tStartOp_;
                
                return true;
            }
            
        } else { // MRU empty, add to the tree and to MRU
            chemistryTree_->clear();
            chemComp *nPhi = NULL;
            chemistryTree_->insertNewLeaf(phi,Rphi,A,scaleFactor_,epsTol_,nSpec_,qrType_,nPhi);
            addToMRU(chemistryTree_->getTreeMin());
            addToMFU(chemistryTree_->getTreeMin());
            nAdd_++;
            tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
            cpuAdd_ += tEndOp_-tStartOp_;
            return true;
        }
        
    } else { // tree not full
         
        chemComp *phi0ISAT = dynamic_cast<chemComp*>(phi0);
        chemistryTree_->insertNewLeaf(phi,Rphi,A,scaleFactor_,epsTol_,nSpec_,qrType_,phi0ISAT);
        phi0 = phi0ISAT;
        // add leaf to MRU and MFU
        if(phi0->getNode()->getElRight()) {
            addToMRU(phi0->getNode()->getElRight());
            addToMFU(phi0->getNode()->getElRight());
        } 
        nAdd_++;

        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
        cpuAdd_ += tEndOp_-tStartOp_;
        return true;
    }   
    
    tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
    cpuAdd_ += tEndOp_-tStartOp_;
    return false;
}
        
        
  
/*
 * grow: check if the leaf can grow 
 * input:
 *      - phi: vector of composition
 *      - Rphi: reaction mapping 
 *      - phi0: closest leaf
 *  
 *  Check if the leaf can grow both considering the growing procedure
 *  and the number of growth already occurred
 *     
 */
bool ISAT::grow(const VectorXd &phi, const VectorXd &Rphi, chemComp*& phi0) 
{
    
    tStartOp_ = double(std::clock())/CLOCKS_PER_SEC;
    
    if(!phi0) {
        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
        cpuGrw_ += tEndOp_-tStartOp_;
        return false;
    }
    chemComp *phi0_=dynamic_cast<chemComp*>(phi0);
    
    if(phi0_->canGrowEOA(phi, Rphi)) 
    {
        // check if nGrow greater than limit if true add to Remove list
        if (double(phi0_->getGrown()) >= double(maxSizeBT_)*maxGrowCoeff_ && flagCleanAndBalance_) 
        {
            addToRemoveList(phi0_);
            tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
            cpuGrw_ += tEndOp_-tStartOp_;
            return false;
        }
        // grow the leaf EOA
        phi0_->growEOA(phi);
        nGrow_++;
        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
        cpuGrw_ += tEndOp_-tStartOp_;
        return true;
    } 
    else 
    {
        tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
        cpuGrw_ += tEndOp_-tStartOp_;
        return false;
    }
        
}

/* interpol: evaluate the mapping gradient using interpolation
 * input:
 *      - phiQ: vector query composition
 *      - phi0: leaf used to interpolating
 * 
 * The function return tha mapping gradient calculated with ISAT interpolation
 * and check if necessary to clean and balance the tree
 */ 
void ISAT::interpol(const VectorXd &phiQ, VectorXd &Rphi,chemComp *&phi0) {

    tStartOp_ = double(std::clock())/CLOCKS_PER_SEC;
    chemComp *phi0_=dynamic_cast<chemComp*>(phi0);
    Rphi = phi0_->getRphi() + phi0_->getMapGrad()*(phiQ-phi0_->getPhi());
    nUse_++;
    tEndOp_ = double(std::clock())/CLOCKS_PER_SEC;
    cpuInt_ += tEndOp_-tStartOp_;

    return;
}

/*
 * cleanLeaf: remove a leaf from the BT and from the lists
 * input:
 *      - phi0: leaf
 */
void ISAT::cleanLeaf(chemComp *&phi0) 
{
    // rimuove la foglia dalla MFU
    mfu_.remove(phi0);
    // rimuove dalla MRU
    mru_.remove(phi0);
    // rimuove la foglia dall'albero (distrugge anche il puntatore)
    chemistryTree_->deleteLeaf(phi0);
    return;
}

/*
 * cleanAndBalance: perform the clean and balance procedure
 * input:
 *      - NO INPUT
 * 
 *  1) Remove leaves that have been grown too many times
 *  2) If the number of add or retrieve is larger BTSize/10 and BTSize/2
 *     respectively enter in the balance procedure: remove the leaf that 
 *     are both old and too leass used
 *  3) If the heigth of the tree is larger than a parameter perform tree balance
 */
bool ISAT::cleanAndBalance() 
{
    
    std::list<chemComp*>::iterator iter;
    
    // 1 - remove leaves that have been grown more times than limit
    for(iter=toRemove_.begin(); iter != toRemove_.end(); ++iter) 
    {
        cleanLeaf(*iter);
    }
    toRemove_.clear();
    
    // 2 - scan the entire tree and add the leaf to remove list if:
    //       a) leaf has been used more times than maxUsedTimes_
    //       b) leaf is too old and it is not frequently used
    unsigned int  numberOfEvents = nAdd_+ nUse_ + nGrow_;
    
    if(flagCleanAndBalance_)
    {
        if( 
          (double(numberOfEvents-nAddFromLastBalance_) > (double(maxSizeBT_)*balanceFactorAdd_))
          //|| 
          //(double(nUse_-nUseFromLastBalance_) > (double(maxSizeBT_)*balanceFactorRet_))
          )
          {
            nUseFromLastBalance_ = nUse_;
            nAddFromLastBalance_ = nAdd_;
            nRemovedLeaves_ = 0;    
            //balance procedure
            std::list<chemComp*> rem;
            double timeNow = double(std::clock())/CLOCKS_PER_SEC;
            chemComp *x = chemistryTree_->getTreeMin();
            double diff = timeNow - timeStart_;
            while(x != NULL) {
                if( 
                  (double(x->getUsed()) < double(mru_.back()->getUsed())*minUsedCoeff_)      // too less used
                  ||                                                                        // and
                  ((timeNow - x->getLastTimeUse()) > maxTimeOldCoeff_*diff)                    // too old
                ) {
                    rem.push_back(x);
                    nRemovedLeaves_++;
                } 
                x = chemistryTree_->treeNextLeaf(x);
            }

	    
            // remove selected leaves
            for(iter=rem.begin(); iter != rem.end(); ++iter) 
            {
                cleanLeaf(*iter);
            }
	    

            // 3 - check the heigth of the tree and balance if needed
            if(double(chemistryTree_->getHeight()) >= (maxHeightCoeff_*std::log(chemistryTree_->getSize())/std::log(2))) 
            {
                chemistryTree_->balance();
            }
            
        nAddFromLastBalance_ = numberOfEvents;
        //nUseFromLastBalance_ = nUse_;
        nBalance_++;
        return true;
        }
    }	
    return false;
}

/*
 * addToMRU: add a leaf to MRU
 * input:
 *      - phi: leaf
 * 
 */
 /*
 void ISAT::addToMRU(chemComp *phi) {

    // check if the leaf is alreasy inside the list
    for(std::list<chemComp*>::iterator iter=mru_.begin(), end=mru_.end(); iter != end; ++iter) {
        if( (*iter)->getIndex() == phi->getIndex()) { //it is in the leaf -> move it to the top
            mru_.erase(iter);
            mru_.push_front(phi); 
            return;
        }
    } 
    
    // the leaf is not inside the list
    if(mru_.size() < maxSizeMRU_) {
            mru_.push_front(phi);
    } else { // mru_Size == maxSizeMRU -> erase last element and add the new one
            mru_.pop_back(); //remove last element
            mru_.push_front(phi);
    }
    return;
        
}
*/
void ISAT::addToMRU(chemComp *phi) {
    std::list<chemComp*>::iterator iter;
    
    // controllo che la foglia sia già nella lista
    bool inMRU = false;
    for(iter=mru_.begin(); iter != mru_.end(); ++iter) 
    {
        if( (*iter) == phi) 
        {
            inMRU = true;
            break;
        }
    } 
    
    if(inMRU)
    { //it in in the leaf -> move it to the top
        mru_.remove(*iter);
        mru_.push_front(phi);
    } 
    else 
    {
        // if MRU size < max size -> add new item (phi)
        if(mru_.size() < maxSizeMRU_) 
        {
            mru_.push_front(phi);
        }
        else 
        { // mru_Size == maxSizeMRU ->cancella ultimo e metti dentro nuovo
            mru_.pop_back(); //remove last element
            mru_.push_front(phi);
        }
        
    }
        
}
/*
 * addToMFU: add a leaf to MFU
 * input:
 *      - phi: leaf
 * 
 */
 /*
void ISAT::addToMFU(chemComp *phi) {
   if(mfu_.size() > 1) {
        for(std::list<chemComp*>::iterator iter=mfu_.begin(), end=mfu_.end(); iter != end; ++iter) {
            if( (*iter)->getIndex() == phi->getIndex()) {
                mfu_.erase(iter);
                std::list<chemComp*>::reverse_iterator k;
                for(k=mfu_.rbegin(); k != mfu_.rend(); ++k) {
                    if((*k)->getUsed() > phi->getUsed()) {
                        break;
                    } 
                }
                mfu_.insert(k.base(),phi);
                return;      
            }
        }
        // leaf is not in the list
        if(mfu_.size() < maxSizeMFU_) {
            std::list<chemComp*>::reverse_iterator k;
            for(k=mfu_.rbegin(); k != mfu_.rend(); ++k) {
                if((*k)->getUsed() > phi->getUsed()) {
                    break;
                } 
            }
            mfu_.insert(k.base(),phi);
            return;
        } else { // list full check if last leaf is lass used, remove and add in the correct position
            if( mfu_.back()->getUsed() < phi->getUsed() ) {
                mfu_.pop_back();
                std::list<chemComp*>::reverse_iterator k;
                for(k=mfu_.rbegin(); k != mfu_.rend(); ++k) {
                    if((*k)->getUsed() > phi->getUsed()) {
                        break;
                    } 
                }
                mfu_.insert(k.base(),phi);
                return;
            }
        }
        return;
    } else if (mfu_.size()==1) {
        if(mfu_.front()->getIndex() != phi->getIndex()) {
            if( mfu_.front()->getUsed() < phi->getUsed()) {
                mfu_.push_front(phi);
                return;
            } else {
                mfu_.push_back(phi);
                return;
            }
        }
        return;
    } else  {
        mfu_.push_front(phi);
        return;
    }

 }
 */
 void ISAT::addToMFU(chemComp* phi) 
 {
    std::list<chemComp*>::iterator iter;
    // controllo che la foglia sia già nella lista
    bool inMFU = false;
    for(iter=mfu_.begin(); iter != mfu_.end(); ++iter) 
    {
        if( (*iter) == phi) 
        {
            inMFU = true;
            break;
        }
    } 
    
    if(inMFU) //it in in the leaf -> move to rigth position
    { 
        mfu_.remove(*iter); // rimuovo dalla lista
 
        if(mfu_.size() > 1) 
        {
            std::list<chemComp*>::iterator i;
            for(i=mfu_.end(); i != mfu_.begin(); ) 
            {
	       i--; //RICCARDO UGLIETTI: ho spostato i-- fuori dalla parentesi del for

               if((*i)->getUsed()<=phi->getUsed()) 
               {
                  break;
               }
            }
            mfu_.insert(i,phi);
         } 
         else 
         {
            mfu_.push_front(phi);
         }
    } 
    else 
    {
        // if MFU size < max size -> add new item (phi)
        if(mfu_.size() < maxSizeMFU_) 
        { // fai scorrere la MFU e inserisci al posto giusto della frequenza
            if(mfu_.size()>0) 
            {
                // trova la posizione in cui inserirlo
                std::list<chemComp*>::iterator i;
                for(i=mfu_.end(); i != mfu_.begin();) 
                {
		   i--; //RICCARDO UGLIETTI: ho spostato i-- fuori dalla parantesi del for
                   
		   if((*i)->getUsed()<=phi->getUsed()) 
                   {
                       break;
                   }
                }
                mfu_.insert(i,phi);
            }
            else 
            {
                mfu_.push_front(phi);
            }
        } 
        else 
        { // mfu_Size == maxSizeMFU -> se ultima foglie usata meno di questa add 
            if(mfu_.back()->getUsed() < phi->getUsed()) 
            {
                mfu_.pop_back(); //remove last element
                mfu_.push_back(phi);
            }
        }
    }
}
/*
 * addToRemoveList: add a leaf to the list of the leaf that will be removed
 * input:
 *      - phi: leaf
 * 
 */
void ISAT::addToRemoveList(chemComp *phi) 
{
    std::list<chemComp*>::iterator iter;
    bool inRemoveList = false;
    for(iter=toRemove_.begin(); iter != toRemove_.end(); ++iter) 
    {
        if( (*iter) == phi) 
        {
            inRemoveList = true;
            break;
        }
    } 
    
    if(!inRemoveList)
    toRemove_.push_back(phi);
        
    return;
}

/*
 * ~ISAT: destructor
 */
ISAT::~ISAT() 
{
    chemistryTree_->clear();
    mfu_.clear();
    mru_.clear();
    delete chemistryTree_;
}
