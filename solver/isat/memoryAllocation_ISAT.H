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

//- Homogeneous
Eigen::VectorXd phiISAT_HOM;
Eigen::VectorXd RphiISAT_HOM;
Eigen::MatrixXd mapGrad_HOM;

//- Heterogeneous
Eigen::VectorXd phiISAT_HET;
Eigen::VectorXd RphiISAT_HET;
Eigen::MatrixXd mapGrad_HET;

if(isatHomoCheck == true)
{
	const unsigned int GAS_NC   = thermodynamicsMapXML->NumberOfSpecies();
	const unsigned int NEQ = GAS_NC+1;

	phiISAT_HOM.resize(NEQ);
	RphiISAT_HOM.resize(NEQ);
	mapGrad_HOM.resize(NEQ,NEQ);
}

if(isatHeteroCheck == true)
{
	const unsigned int GAS_NC   = thermodynamicsMapXML->NumberOfSpecies();
	const unsigned int SURF_NC  = surfaceThermodynamicsMapXML->number_of_site_species();
	const unsigned int SURF_NP  = surfaceThermodynamicsMapXML->number_of_site_phases(0);
	const unsigned int NEQ = GAS_NC+1+SURF_NP+SURF_NC+1;

	phiISAT_HET.resize(NEQ);
	RphiISAT_HET.resize(NEQ);
	mapGrad_HET.resize(NEQ,NEQ);
}
