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

	catalyticReactorClass  *catalyticReactors;
	labelList       		isNeighbourCell;
	labelList				indexCatalyst;

	label totalNumberOfNonCatalyticReactors = 0;
	label totalNumberOfCatalyticFaces    	= 0;
	label totalNumberOfCatalyticReactors 	= 0;
	label totalNumberOfReactors				= 0;
	label maxNumberOfCatalyticFaces 		= 0;
	label minNumberOfCatalyticFaces 		= 0;
	
	isNeighbourCell.resize(Y[0].size());
	for (int j=0;j<Y[0].size();j++)
		isNeighbourCell[j] = 0;

	if(catalyst == true)
	{
		Info << " Catalyst topology Analysis" << endl;
		
		labelList       catalyticWallSize;
		labelList       catalystNeighbour;
		labelList       catalyticPatchIndex;
		labelList       catalyticFaceIndex;
		scalarField     catalyticCellArea;
		scalarList      catalyticNeighbourCellVolume;

		//- Set gradients to zero & find catalytic cells
		forAll(Y[0].boundaryField(), patchi)
		{
			/*forAll(Y[0].boundaryFieldRef()[patchi], facei)
			{
				for (label i=0; i<Ygrad.size(); i++)
				{
					Ygrad[i].boundaryFieldRef()[patchi][facei] = 0.0;
				}
				
				grad_T.boundaryFieldRef()[patchi][facei] = 0.0;
			}*/
    
				forAll(catalyticWalls,patchj)
				{
					if (patchi == mesh.boundaryMesh().findPatchID(catalyticWalls[patchj]))
					{
						indexCatalyst.append(patchi);
						catalyticWallSize.append(Y[0].boundaryFieldRef()[patchi].size());
						totalNumberOfCatalyticFaces += Y[0].boundaryFieldRef()[patchi].size();
					}
				}
		}
		

		//- Look for cells near catalytic wall
		if(totalNumberOfCatalyticFaces != 0)
		{
			Info << endl;
			Info << "-------------------------------------------------------------------------" << endl;
			Info << "                        Catalytic patches                                " << endl;
			Info << "-------------------------------------------------------------------------" << endl;
			Info << " Patch\tFaces\tmaxA\t\tminA\t\tmeanA\t\ttotA                           " << endl;
			Info << "-------------------------------------------------------------------------" << endl;
			Info<<totalNumberOfCatalyticFaces << indexCatalyst.size() << catalyticWallSize.size() << endl;
		
			for(int i=0; i<indexCatalyst.size(); i++)
			{
				Info<<totalNumberOfCatalyticFaces << indexCatalyst[i] << catalyticWallSize[i] << endl;
			}
			
			{
				OpenSMOKE::OpenSMOKEVectorInt    		catalystNeighbourTemp(totalNumberOfCatalyticFaces);
				OpenSMOKE::OpenSMOKEVectorInt    		catalyticPatchIndexTemp(totalNumberOfCatalyticFaces);
				OpenSMOKE::OpenSMOKEVectorInt    		catalyticFaceIndexTemp(totalNumberOfCatalyticFaces);
				OpenSMOKE::OpenSMOKEVectorDouble		catalyticCellAreaTemp(totalNumberOfCatalyticFaces);
				OpenSMOKE::OpenSMOKEVectorDouble       	catalyticNeighbourCellVolumeTemp(totalNumberOfCatalyticFaces);

				label counter = 1;
				for(int i=0; i<indexCatalyst.size(); i++)
				{
					double totalArea = 0.;
					double maxArea   = -1.e16;
					double minArea   =  1e16;

					const labelList& cellOwner = Y[0].boundaryFieldRef()[indexCatalyst[i]].patch().faceCells();
					forAll(Y[0].boundaryFieldRef()[indexCatalyst[i]] , facei)
					{

						double area = mesh.magSf().boundaryField()[indexCatalyst[i]][facei];
						
						catalyticPatchIndexTemp[counter]		= indexCatalyst[i];
						catalyticFaceIndexTemp[counter]		= facei;   
						catalystNeighbourTemp[counter]		= cellOwner[facei];
						catalyticNeighbourCellVolumeTemp[counter]	= mesh.V()[cellOwner[facei]];
						catalyticCellAreaTemp[counter]		= area;

						totalArea += area;
						if (area>maxArea) 		maxArea = area;
						else if (area < minArea) 	minArea = area;

						counter++;
					}
					
					if(catalyticWallSize[i]!=0)
					{
						Info << " " 	<< indexCatalyst[i] << "\t" << catalyticWallSize[i] << "\t" << maxArea << "\t" 
							<< minArea << "\t" << totalArea/double(catalyticWallSize[i]) << "\t" <<  totalArea << endl;
					}
					else
					{
						Info << " "     << indexCatalyst[i] << "\tno catalytic area for this patch in master processor" << endl;
					}
					
				}

				Info << endl << endl;

				// Sort provisional vectors
				OpenSMOKE::OpenSMOKEVectorInt Indices(totalNumberOfCatalyticFaces);
				Sort(&catalystNeighbourTemp, &Indices);
				Reorder(&catalyticFaceIndexTemp, Indices);
				Reorder(&catalyticPatchIndexTemp, Indices);
				Reorder(&catalyticCellAreaTemp, Indices);
				Reorder(&catalyticNeighbourCellVolumeTemp, Indices);

				// Indices
				catalyticFaceIndex.resize(totalNumberOfCatalyticFaces);
				catalyticPatchIndex.resize(totalNumberOfCatalyticFaces);
				catalyticCellArea.resize(totalNumberOfCatalyticFaces);
				catalystNeighbour.resize(totalNumberOfCatalyticFaces);
				catalyticNeighbourCellVolume.resize(totalNumberOfCatalyticFaces);

				for (int j=1;j<=Indices.Size();j++)
				{
					catalyticFaceIndex[j-1]  		= catalyticFaceIndexTemp[j];
					catalyticPatchIndex[j-1] 		= catalyticPatchIndexTemp[j];
					catalyticCellArea[j-1]   		= catalyticCellAreaTemp[j];
					catalystNeighbour[j-1]   		= catalystNeighbourTemp[j];
					catalyticNeighbourCellVolume[j-1]   	= catalyticNeighbourCellVolumeTemp[j];
				}



				// NeighbourCell flag
				int k=0;	    		
				while(k<catalystNeighbour.size())
					{
					int counter = 1;
					for(int j=k+1; j<catalystNeighbour.size(); j++)
					{
						if (catalystNeighbour[k] == catalystNeighbour[j])
							counter++;
						else
							break;
					}

					isNeighbourCell[catalystNeighbour[k]] = counter;
					
					k+=counter;
				}

				// Count total number of catalytic reactors
				scalarField& TCells = T.ref();
				forAll(TCells, celli)
				{
					if (isNeighbourCell[celli] > 0)
						totalNumberOfCatalyticReactors++;
				}
			}

			Info << "Total number of catalytic reactors: " << totalNumberOfCatalyticReactors << endl;
			catalyticReactors = new catalyticReactorClass[totalNumberOfCatalyticReactors];
		

			int k=0;	    		
			int index = 0;
			while(k<catalystNeighbour.size())
			{
				int counter = 1;
				for(int j=k+1; j<catalystNeighbour.size(); j++)
				{
					if (catalystNeighbour[k] == catalystNeighbour[j])
						counter++;
					else
						break;
				}

				for(int j=k; j<k+counter; j++)
				{
					catalyticReactors[index].globalIndex = catalystNeighbour[j];
					catalyticReactors[index].volume = catalyticNeighbourCellVolume[j];
					catalyticReactors[index].numberOfFaces++;
					catalyticReactors[index].patchIndex.Append(catalyticPatchIndex[j]);
					catalyticReactors[index].faceIndex.Append(catalyticFaceIndex[j]);
					catalyticReactors[index].faceArea.Append(catalyticCellArea[j]);	
				}
				
				k+=counter;
				index++;
			}

			totalNumberOfReactors = T.size();
			totalNumberOfNonCatalyticReactors = totalNumberOfReactors - totalNumberOfCatalyticReactors;

			// Summary on video
			{
				double meanNumberOfCatalyticFaces = 0;
				maxNumberOfCatalyticFaces 	= 0;
				minNumberOfCatalyticFaces 	= 1e5;

				for(int j=0;j<totalNumberOfCatalyticReactors;j++)
				{
					if      (catalyticReactors[j].numberOfFaces > maxNumberOfCatalyticFaces) maxNumberOfCatalyticFaces = catalyticReactors[j].numberOfFaces;
					else if (catalyticReactors[j].numberOfFaces < minNumberOfCatalyticFaces) minNumberOfCatalyticFaces = catalyticReactors[j].numberOfFaces;
					meanNumberOfCatalyticFaces += catalyticReactors[j].numberOfFaces;
				}
				meanNumberOfCatalyticFaces /= double(totalNumberOfCatalyticReactors);


				Info << " * Total number of reactors:              " << totalNumberOfReactors << endl;
				Info << " * Total number of catalytic reactors:    " << totalNumberOfCatalyticReactors << " (" << double(totalNumberOfCatalyticReactors)/double(totalNumberOfReactors)*100. << " %)" << endl;
				Info << " * Total number of noncatalytic reactors: " << totalNumberOfNonCatalyticReactors << " (" << double(totalNumberOfNonCatalyticReactors)/double(totalNumberOfReactors)*100. << " %)" << endl;
				Info << endl;
				Info << " * Cat. faces belonging to cat. reactors: Max " << maxNumberOfCatalyticFaces << " Min " << minNumberOfCatalyticFaces << " Mean " << meanNumberOfCatalyticFaces << endl;
			}

			// Write on file
			{
				std::ofstream fCatalyticReactors;
				fCatalyticReactors.open("CatalyticReactors", std::ios::out);
				fCatalyticReactors << "#\tpatch\tcellIndex\tedges\tvolume\tfaceIndex\tfaceArea" << endl;
				for(int i=0;i<totalNumberOfCatalyticReactors;i++)
				{
						fCatalyticReactors	<< i 						<< "\t"
											<< catalyticReactors[i].globalIndex 		<< "\t"
											<< catalyticReactors[i].volume 		<< "\t" 
											<< catalyticReactors[i].numberOfFaces	<< "\t"
											<< endl;
				}
				fCatalyticReactors.close();
			}
		}
		else
		{
			totalNumberOfNonCatalyticReactors = T.size();
		}
	}
	else
	{
		Info << "No catalytic reactors: " << endl;
		catalyticReactors = new catalyticReactorClass;
	}
