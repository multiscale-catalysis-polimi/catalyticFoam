catalyticFoam
============
CFD solver for heterogeneous reacting flows with detailed kinetic mechanisms

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

## Information:
**catalyticFoam** has been developed on top of the OpenFOAM framework and it is currently compatible with OpenFOAM version from 4.x to 9. 
**catalyticFoam** is not part of the official OpenFOAM release nor endorsed by The OpenFOAM Foundation.

## Compulsory libraries:
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) 
- [RapidXML](http://rapidxml.sourceforge.net/)
- [Boost C++](http://www.boost.org/)
- [OpenSMOKE++][1] (provided with the current version of catalyticFoam)

## Compilation on Linux:
Instructions:
1. Open a terminal and move to the desidered position of the catalyticFoam folder
2. Type: `git clone https://github.com/multiscale-catalysis-polimi/catalyticFoam.git`
3. Open the `mybashrc`and adjust the paths to the compulsory external libraries
4. Type: `source mybashrc`
5. Type: `./Allwmake`

## Compilation on Windows 10:
1. Install Windows Subsystem for Linux (WSL) following the instructions available [here]( https://docs.microsoft.com/en-gb/windows/wsl/install-win10)
2. Install OpenFOAM on WSL following the instructions available [here](https://openfoam.org/download/)
3. Open a terminal and move to the desidered position of the catalyticFoam folder
4. Type: `git clone https://github.com/multiscale-catalysis-polimi/catalyticFoam.git`
5. Open the `mybashrc` and adjust the paths to the compulsory external libraries
6. Type: `source mybashrc`
7. Type: `./Allwmake`

## Run your first case:
The folder `example/case` contains a simple test case (2D honeycomb channel).
1. Build the mesh using the `blockMesh` utility, 
2. Run the case using the `catalyticPimpleFoam` solver. 

The folder `example/case_isat` contains a simple test case (2D honeycomb channel) simulated with the ISAT library.
1. Build the mesh using the `blockMesh` utility, 
2. Run the case using the `catalyticPimpleFoam` solver. 

The folder `example/sphere` contains a 4-spheres string reactor.
1. Build the mesh with `./makeMesh`
2. Run the case using the `catalyticPimpleFoam` solver. 

## Preprocessing of a kinetic scheme:
The preprocessing of a kinetic scheme requires:
1. Create a new folder with the kinetic scheme
2. Create the kinetic.kin, surface.sur according to the kinetic mechanism
3. Provide the termodynamic (thermo.tdc) and transport (transport.tra) databases
4. Create the input.dic input file (see example)
5. Compile the mechanism with `catalyticFoam_CHEMKINPreProcessor`


[1]: https://www.opensmokepp.polimi.it
