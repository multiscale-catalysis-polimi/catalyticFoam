/*------------------------------------------------------------------------\
|   catalyticFOAM                                                         |
|   http://www.catalyticfoam.polimi.it/                                   |
|                                                                         |
|   Authors:                                                              |
|                                                                         |
|   Alberto Cuoci <alberto.cuoci@polimi.it>                               |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|   Matteo Maestri <matteo.maestri@polimi.it>                             |
|   Department of Energy                                                  |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|   Mauro Bracconi <mauro.bracconi@polimi.it>                             |
|   Department of Energy                                                  |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of catalyticFOAM framework.                         |
|                                                                         |
|   Copyright(C) 2020-2011, A.Cuoci, M.Maestri,                           |
|                2020-2014, M. Bracconi                                   |
|                2015-2013, S.Rebughini                                   |
|                     2013, T.Maffei                                      |
|                     2013, G.Gentile, F.Manelli                          |
|                     2012, M.Calonaci, F.Furnari                         |
|                     2011, S.Goisis, A.Osio                              |
|                                                                         |
|   catalyticFOAM is distributed in the hope that it will be useful,      |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with catalyticFOAM. If not, see <http://www.gnu.org/licenses/>. |
|                                                                         |
\*-----------------------------------------------------------------------*/

{
	Info << endl;
	Info << endl;
	Info << "/-------------------------------------------------------------------------\\" << endl;
	Info << "|   catalyticFOAM                                                         |" << endl;
	Info << "|   http://www.catalyticfoam.polimi.it/                                   |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Authors:                                                              |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Alberto Cuoci <alberto.cuoci@polimi.it>                               |" << endl;
	Info << "|   Department of Chemistry, Materials and Chemical Engineering           |" << endl;
	Info << "|   Politecnico di Milano                                                 |" << endl;
	Info << "|   P.zza Leonardo da Vinci 32, 20133 Milano                              |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Matteo Maestri <matteo.maestri@polimi.it>                             |" << endl;
	Info << "|   Department of Energy                                                  |" << endl;
	Info << "|   Politecnico di Milano                                                 |" << endl;
	Info << "|   P.zza Leonardo da Vinci 32, 20133 Milano                              |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Mauro Bracconi <mauro.bracconi@polimi.it>                             |" << endl;
	Info << "|   Department of Energy                                                  |" << endl;
	Info << "|   Politecnico di Milano                                                 |" << endl;
	Info << "|   P.zza Leonardo da Vinci 32, 20133 Milano                              |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Copyright(C) 2020-2011, A.Cuoci, M.Maestri,                           |" << endl;
	Info << "|                2020-2014, M.Bracconi                                    |" << endl;
	Info << "|                2018-2013, S.Rebughini                                   |" << endl;
	Info << "|                     2013, T.Maffei                                      |" << endl;
	Info << "|                     2013, G.Gentile                                     |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|-------------------------------------------------------------------------|" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Some parts of catalyticFOAM are based on OpenFOAM® and Copyright on   |" << endl;
	Info << "|   these parts is held by the OpenFOAM® Foundation (www.openFoam.org)    |" << endl;
	Info << "|   and potentially other parties.                                        |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   Note: this code is not part of OpenFOAM®.                             |" << endl;
	Info << "|   This offering is not approved nor endorsed by OpenCFD Limited,        |" << endl;
	Info << "|   the producer of the OpenFOAM® software and owner of the OPENFOAM®     |" << endl;
	Info << "|   trade mark.                                                           |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|-------------------------------------------------------------------------|" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   If you are using this software, please cite:                          |" << endl;
	Info << "|   *   Maestri, M. and Cuoci, A. (2013)                                  |" << endl;
	Info << "|   *   Coupling CFD with detailed microkinetic modeling                  |" << endl;
	Info << "|       in heterogeneous catalysis.                                       |" << endl;
	Info << "|   *   Chemical Engineering Science. Volume 96. Pages 106-117.           |" << endl;
	Info << "|       doi.org/10.1016/j.ces.2013.03.048.                                |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "|   If you are using ISAT within this software, please cite:              |" << endl;
	Info << "|   *   Bracconi, M., Maestri, M. and Cuoci, A. (2017)                    |" << endl;
	Info << "|   *   In situ adaptive tabulation for the CFD simulation of             |" << endl;
	Info << "|       heterogeneous reactors based on operator-splitting algorithm.     |" << endl;
	Info << "|   *   AIChE J., Volume 63. Pages 95–104.                                |" << endl;
	Info << "|       doi:10.1002/aic.15441                                             |" << endl;
	Info << "|                                                                         |" << endl;
	Info << "\\-------------------------------------------------------------------------/" << endl;
	Info << endl;
	Info << endl;
}


