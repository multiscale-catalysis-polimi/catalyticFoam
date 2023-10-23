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


#ifndef ISAT_H
#define    ISAT_H

#if ISAT_USE_MKL == 1
    #define EIGEN_USE_MKL_ALL 
#endif

#include "binaryTree.h"
#include "binaryNode.h"
#include "Eigen/Core"
#include <ctime>
#include <list>

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;

class ISAT {
private:
    // ISAT tree
    binaryTree *chemistryTree_;
    
    // vector scaleFactor
    VectorXd scaleFactor_;
   
    // epsTol
    double epsTol_ ;
    
    // number of species
    unsigned int nSpec_;
    
    // MRU list
    std::list<chemComp*> mru_;
    // MFU list
    std::list<chemComp*> mfu_;
    // toRemove list
    std::list<chemComp*> toRemove_;
    
    //----------flag------------------
    bool flagSearchBTSR_;
    bool flagSearchMRU_;
    bool flagSearchBruteForce_;
    bool flagSearchMFU_;
    bool flagSearchCL_;
    bool flagClearingIfFull_;
    bool flagCleanAndBalance_;

    //----------balance param --------
    double balanceFactorAdd_;
    double balanceFactorRet_;
    double minUsedCoeff_;
    double maxGrowCoeff_;     
    double maxTimeOldCoeff_;
    double maxHeightCoeff_;    


    //----------maxSize----------------
    // MRU max items
    unsigned int maxSizeMRU_;
    // MFU max items
    unsigned int maxSizeMFU_;
    // MRU max search items;
    unsigned int maxSearchMRU_;
    // MFU max search items;
    unsigned int maxSearchMFU_;
    // BTSR max level search
    unsigned int maxLevelBTSR_;
    // BT max size
    unsigned int maxSizeBT_;
    //---------------------------------
    
    //----------counter----------------
    // number of leaf addition
    unsigned int nAdd_;
    // number of leaf growing
    unsigned int nGrow_;
    // number of leaf use
    unsigned int nUse_;
    
    // number of successfully MRU search
    unsigned int nMRU_;
    // number of successfully MFU search
    unsigned int nMFU_;
    // number of successfully BTS search
    unsigned int nBTS_;
    // number of successfully BTSR search
    unsigned int nBTSR_;
    // number of successfully CL search
    unsigned int nCL_;
    // number of fail binary tree search 
    unsigned int nFailBTinEOA_;

    // number of balance procedure
    unsigned int nBalance_;
    unsigned int nRemovedLeaves_;

    // counter of events between balance
    unsigned int nAddFromLastBalance_;
    unsigned int nUseFromLastBalance_;
    //----------------------------------

    //----------time info---------------
    // time info
    double tStartOp_;
    double tEndOp_;
    // start time
    double timeStart_;
    // time spent in retrieving
    double cpuRet_;
     // time spent in interpolation
    double cpuInt_;
     // time spent in growing
    double cpuGrw_;
     // time spent in adding
    double cpuAdd_;
    //----------------------------------

    //------------qr type---------------
    unsigned int qrType_;

public:
    // costructor
    ISAT(VectorXd scaleFactor, double epsTol, unsigned int nCol);

    // find the closest leaf to a given query point
    bool retrieve(const VectorXd &phiQ, chemComp *&closest);
    
    // evaluate the mapping gradient using interpolation
    void interpol(const VectorXd &phiQ, VectorXd &Rphi,chemComp *&phi0);
    
    // add a new leaf in the ISAT tree
    bool add(const VectorXd &phi, const VectorXd &Rphi, const MatrixXd &A, chemComp *&phi0);
    
    // check if the leaf has an EOA which can grow
    bool grow(const VectorXd &phi, const VectorXd &Rphi, chemComp *&phi0);
    
    // check if it is necessary to balance and do it
    bool cleanAndBalance();
   
    // destructor
    virtual ~ISAT();
    
    //----------------Access Function-------------------------------
    // return nUse
    inline int nUse() const {
        return nUse_;
    }

     // return nAdd
    inline int nAdd() const {
        return nAdd_;
    }
    
     // retrun nGrow
    inline int nGrow() const {
        return nGrow_;
    }
    
    // retrun nMRU
    inline int nMRU() const {
        return nMRU_;
    }
    
     // retrun nMFU
    inline int nMFU() const {
        return nMFU_;
    }
    
    // retrun nBTS
    inline int nBTS() const {
        return nBTS_;
    }

    // return nBTSR
    inline int nBTSR() const {
	return nBTSR_;
    }

    // return nCL
    inline int nCL() const {
	return nCL_;
    }

    // retrun nBalance
    inline int nBalance() const {
        return nBalance_;
    }
    
    // retrun maxSizeMRU
    inline int getMaxSizeMRU() const {
        return maxSizeMRU_;
    }
    
    // retrun maxSizeMRU
    inline int getMaxSizeMFU() const {
        return maxSizeMFU_;
    }
    
    // retrun maxSizeBT
    inline int getMaxSizeBT() const {
        return maxSizeBT_;
    }
    
        // retrun maxSizeMRU
    inline int getMaxSearchMRU() const {
        return maxSearchMRU_;
    }
    
    // retrun maxSizeMRU
    inline int getMaxSearchMFU() const {
        return maxSearchMFU_;
    }

    // return maxLevelBTSR
    inline int getMaxLevelBTSR() const {
	return maxLevelBTSR_;
    }
    
    // retrun nFailBTinEOA
    inline int nFailBTinEOA() const {
        return nFailBTinEOA_;
    }

    // return flagBTSRSearch
    inline bool getFlagSearchBTSR() const {
	return flagSearchBTSR_;
    }
   
    // return flagMRUSearch
    inline bool getFlagSearchMRU() const {
        return flagSearchMRU_;
    }
    
    // return flagMFUSearch
    inline bool getFlagSearchMFU() const {
        return flagSearchMFU_;
    }
    
    // return flagSearchBruteForce
    inline bool getFlagSearchBruteForce() const {
        return flagSearchBruteForce_;
    }

    // return flagSearchCL
    inline bool getFlagSearchCL() const {
	return flagSearchCL_;
    }

    // return clearingIfFull
    inline bool getFlagClearingIfFull() const {
        return flagClearingIfFull_;
    }

    // return flagCleanAndBalance
    inline bool getFlagCleaningAndBalance() const {
        return flagCleanAndBalance_;
    }

    // return MaxTimeOld (balance param)
    inline double getMaxTimeOldCoeff() const {
        return maxTimeOldCoeff_;
    }

    // return maxGrowCoeff (balance param)
    inline double getMaxGrowCoeff() const {
        return maxGrowCoeff_;
    }

    // return maxUsedCoeff (balance param)
    inline double getMinUsedCoeff() const {
        return minUsedCoeff_;
    }

    // return maxHeightCoeff (balance param)
    inline double getMaxHeightCoeff() const {
        return maxHeightCoeff_;
    }
    
    // return balanceFactorAdd (balance param)
    inline double getBalanceFactorAdd() const {
        return balanceFactorAdd_;
    }

    // return balanceFactorRet (balance param)
    inline double getBalanceFactorRet() const {
        return balanceFactorRet_;
    }
    
    // return scaling factors
    inline VectorXd getScaleFactor() const {
        return scaleFactor_;
    }
    
    // retrun the qrType
    inline unsigned int getQRType() const {
        return qrType_;
    }
    
    // return average retrieve time
    inline double getCpuRetrieveTime() const {
        if(nUse() > 0)
            return cpuRet_/double(nUse());
        else
         return 0.;
            
    }
    
    // return average interpolation time
    inline double getCpuInterpolationTime() const {
        if(nUse() > 0)
            return cpuInt_/double(nUse());
        else
            return 0.;
    }
    
    // return average grow time
    inline double getCpuGrowTime() const {
        if(nGrow() > 0)
            return cpuGrw_/double(nGrow());
        else
            return 0.;
    }
    
    // return average add time
    inline double getCpuAddTime() const {
        if(nAdd() > 0)
            return cpuAdd_/double(nAdd());
        else
            return 0.;
    }

    // return pointer to binary tree 
    inline binaryTree *getBTree() const {
        return chemistryTree_;
    }

    // return MRU list
    inline std::list<chemComp*> getMRUList() const {
        return mru_;
    }
   
    // return MFU list
    inline std::list<chemComp*> getMFUList() const {
        return mfu_;
    }
    
    // return ISAT tolerance
    inline double getEpsTol() const {
        return epsTol_;
    }
    
    // removed leaves in balance operation
    inline double getNumberRemovedLeaves() const {
        return nRemovedLeaves_;
    }

    //----------------Set ISAT param-------------------------------
    // set maxSizeBT
    inline void setMaxSizeBT(int maxSizeBT) {
        maxSizeBT_ = maxSizeBT;
    }
    
    // set maxSizeMRU
    inline void setMaxSizeMRU(int maxSizeMRU) {
        maxSizeMRU_ = maxSizeMRU;
    }

    // set maxSizeMFU
    inline void setMaxSizeMFU(int maxSizeMFU) {
    maxSizeMFU_ = maxSizeMFU;
    }

    // set maxSearchMFU
    inline void setMaxSearchMFU(int maxSearchMFU) {
    maxSearchMFU_ = maxSearchMFU;
    }

    // set maxSearchMRU
    inline void setMaxSearchMRU(int maxSearchMRU) {
    maxSearchMRU_ = maxSearchMRU;
    }

    // set maxLevelBTSR
    inline void setMaxLevelBTSR(int maxLevel) {
    maxLevelBTSR_ = maxLevel;
    }

    // set flagBTSRSearch
    inline void setFlagSearchBTSR(bool flagBTSR) {
    	flagSearchBTSR_ = flagBTSR;
    }

    // set flagSearchMRU
    inline void setFlagSearchMRU(bool flagSearchMRU) {
        flagSearchMRU_ = flagSearchMRU;
    }
    
    // set flagMFUSearch
    inline void setFlagSearchMFU(bool flagSearchMFU) {
        flagSearchMFU_ = flagSearchMFU;
    }
    
    // set flagSearchBruteForce
    inline void setFlagSearchBruteForce(bool flagSearchBruteForce) {
        flagSearchBruteForce_ = flagSearchBruteForce;
    }

    // set flagSearchCL
    inline void setFlagSearchCL(bool flagCL) {
	flagSearchCL_ = flagCL;
    }

    // set clearingIfFull
    inline void setFlagClearingIfFull(bool flagClearingIfFull) {
        flagClearingIfFull_ = flagClearingIfFull;
    }

    // set flagCleanAndBalance
    inline void setFlagCleanAndBalance(bool flagCleanAndBalance) {
        flagCleanAndBalance_ = flagCleanAndBalance;
    }

    // set qrType
    inline void setQRType(unsigned int qrType) {
        qrType_ = qrType;
    }

    // set maxTimeOld (balance param)
    inline void setMaxTimeOldCoeff(double maxTimeOldCoeff) {
        maxTimeOldCoeff_ = maxTimeOldCoeff;
    }

    // set maxHeightCoeff (balance param)
    inline void setMaxHeightCoeff(double maxHeightCoeff) {
        maxHeightCoeff_ = maxHeightCoeff;
    }

    // set minUsedCoeff
    inline void setMinUsedCoeff(double minUsedCoeff) {
        minUsedCoeff_ = minUsedCoeff;
    }

    // set maxGrowCoeff
    inline void setMaxGrowCoeff(double maxGrowCoeff) {
        maxGrowCoeff_ = maxGrowCoeff;
    }
    
    // set balanceFactorAdd
    inline void setBalanceFactorAdd(double balanceFactorAdd) {
        balanceFactorAdd_ = balanceFactorAdd;
    }
    
    // set maxGrowCoeff
    inline void setBalanceFactorRet(double balanceFactorRet) {
        balanceFactorRet_ = balanceFactorRet;
    }

    //--------------------------Reset Stats -------------------------    
    inline void resetStats() {
        nUse_ = 0;
        nAdd_ = 0;
        nGrow_ = 0;
        nBTS_ = 0;
	nBTSR_ = 0;
        nMRU_ = 0;
        nMFU_ = 0;
        nFailBTinEOA_ = 0;
    }
    
private:
    // add leaf to MRU list
    void addToMRU(chemComp *phi);
    // add leaf to MFU list
    void addToMFU(chemComp *phi);
    // add leaf to remove list
    void addToRemoveList(chemComp *phi);
    // remove leaf from tree and from list
    void cleanLeaf(chemComp *&phi0);
    
};

#include "ISAT.C"

#endif    /* ISAT_H */

