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
|   Evaluation of the mapping gradient for catalyticFOAM                  |
|   Compulsory                                                            |                
|   - Eigen Library                                                       |
|   - function to evaluate the Jacobian of the ODE system                 |
|                                                                         |
\*-----------------------------------------------------------------------*/


#include "Eigen/Dense"
 #include <iostream>

/*
 * calcMappingGradient: compute mapping gradient using an implicit method
 * input:
 *      - phiQ: vector of composition 
 *      - RphiQ: vector of composition oafter integration
 *      - A: mapping gradient (output)
 *      - sF: scale factors vector
 *      - typeSolver: type of solver used 
 *              - 0: LU Full Pivoting factorization
 *              - 1: LU Partial Pivoting factorization
 *      - nsubsteps: number of substeps
 *              - 0: evaluation of mapping gradient on initial point (phiQ)
 *              - 1: evaluation of mapping gradient on final point (RphiQ)
 *              - >1: evaluation of mapping gradient using nsubstesp step
 *      - dt: integration time step
 *      - ode: ODE system
 * 
 * The function evaluate the mappingGradient using an implicit integration method using a user fixed number
 * of steps
 */ 
template<typename ODESystem>
void calcMappingGradient(
                        const Eigen::VectorXd &phi0, 
                        const Eigen::VectorXd &phiF, 
                        Eigen::MatrixXd &A, 
                        const Eigen::VectorXd &sF, 
                        const unsigned int typeSolver, 
                        const unsigned int nsubsteps, 
                        const double dt, 
                        ODESystem *ode
                        )
{
    const unsigned int NEQ  = phi0.rows();
    Eigen::VectorXd omega0  = phi0.cwiseQuotient(sF);
    Eigen::VectorXd omegaF  = phiF.cwiseQuotient(sF);

    // Jacobian matrix
    Eigen::MatrixXd J(NEQ,NEQ);

    // Initial condition (Identity matrix, Bii=0)
    Eigen::MatrixXd B(NEQ,NEQ);
    B.setIdentity(NEQ,NEQ);

    // Linear solvers (only declaration)
    Eigen::FullPivLU<Eigen::MatrixXd>    *luFull;
    Eigen::PartialPivLU<Eigen::MatrixXd>    *luPartial;

    // Linear system matrix
    Eigen::MatrixXd M(NEQ,NEQ);

    if (nsubsteps==0) // evaluation on initial
    {
        // Update the Jacobian
        numericalJacobian(omega0, J ,sF, ode);
          // Assembling of linear system matrix
        M = -J*dt;
        for(unsigned int i=0;i<NEQ;i++)
            M(i,i) += 1.;
        
        // LU factorization
        if(typeSolver == 0)
        {
            luFull = new Eigen::FullPivLU<Eigen::MatrixXd>();
            luFull->compute(M);
        }
        else if (typeSolver == 1)
        {
            luPartial = new Eigen::PartialPivLU<Eigen::MatrixXd>();
            luPartial->compute(M);
        }

        // solve with LU factorization
        if(typeSolver == 0)     
            M = luFull->solve(B);
        else            
            M = luPartial->solve(B);
            
        B = M;
    }
    else //final o substeps
    {
        for (unsigned int k=1;k<=nsubsteps;k++)
        {
            // Update the Jacobian
            Eigen::VectorXd omega = omega0+(omegaF-omega0)*double(k)/double(nsubsteps);
            
            numericalJacobian(omega, J ,sF, ode);
            
            // Assembling of linear system matrix
            M = -J*dt/double(nsubsteps);
            for(unsigned int i=0;i<NEQ;i++)
                M(i,i) += 1.;
            
            // LU factorization
            if(typeSolver == 0)
            {
                luFull = new Eigen::FullPivLU<Eigen::MatrixXd>();
                luFull->compute(M);
            }
            else if (typeSolver == 1)
            {
                luPartial = new Eigen::PartialPivLU<Eigen::MatrixXd>();
                luPartial->compute(M);
            }

            // solve with LU factorization
            if(typeSolver == 0)     
                M = luFull->solve(B);
            else            
                M = luPartial->solve(B);
            
            B = M;
        
            }
    }
    if(typeSolver == 0) delete luFull;
    if(typeSolver == 1) delete luPartial;
    A = B;
};

