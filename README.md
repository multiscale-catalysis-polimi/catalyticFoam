catalyticFoam
============
CFD solver for heterogeneous reacting flows with detailed kinetic mechanisms based on OpenFOAM 

If you are using this software, please cite:
> Maestri, M. and Cuoci, A. (2013) Coupling CFD with detailed microkinetic modeling
> in heterogeneous catalysis. Chemical Engineering Science. Volume 96. Pages 106-117.
> doi.org/10.1016/j.ces.2013.03.048.
                                                                      
If you are using ISAT within this software, please cite:
> Bracconi, M., Maestri, M. and Cuoci, A. (2017) In situ adaptive tabulation for the CFD simulation of
> heterogeneous reactors based on operator-splitting algorithm.AIChE J., Volume 63. Pages 95–104.
> doi:10.1002/aic.15441 

## Authors:
**catalyticFoam** has been developed in the Multiscale Catalysis Group of the Laboratory of Catalysis and Catalytic Processes of Politecnico di Milano.

## Compulsory libraries:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [RapidXML](http://rapidxml.sourceforge.net/)
- [Boost C++](http://www.boost.org/)
- [OpenSMOKE++][1] (provided with the current version of catalyticFoam)

## Compilation
Instructions:
1. Open the `mybashrc`, choose the version of OpenFOAM you are using (4.x = 40, 5.x = 50) and adjust the paths to the compulsory external libraries
2. Type: `source mybashrc`
3. Type: `./Allwmake`

## Run your first case
The folder `example/case` contains a simple test case (2D honeycomb channel).

1. Build the mesh using the `blockMesh` utility, 
2. Run the case using the `catalyticPimpleFoam` solver. 

The folder `example/case_isat` contains a simple test case (2D honeycomb channel) simulated with the ISAT library.

1. Build the mesh using the `blockMesh` utility, 
2. Run the case using the `catalyticPimpleFoam` solver. 

The folder `example/sphere` contains a 4-spheres string reactor.
1. Build the mesh with `./makeMesh`
2. Run the case using the `catalyticPimpleFoam` solver. 

[1]: https://www.opensmokepp.polimi.it
