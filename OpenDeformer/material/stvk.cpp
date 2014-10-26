#include "stdafx.h"
#include "stvk.h"
#include "nodeIndexer.h"

namespace ODER{
	StVKMaterial::StVKMaterial(double rho, double lameFirst, double lameSecond, const Reference<Mesh> &mesh, int orderNum)
		:MecMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear)),
		lambda(lameFirst), mu(lameSecond), orders(orderNum){

		D[0] = lambda + 2.0 * mu;
		D[1] = lambda;
		D[2] = mu;

		intergration[0] = NULL;
		intergration[1] = NULL;


		loadParameter = allocAligned<double>(orders);
	}

	void StVKMaterial::generateStiffnessMatrix(const Reference<Mesh> &mesh, const NodeIndexer &indexer, SparseMatrixAssembler& matrix) const{
		const int numNodesPerElement = mesh->numNodesPerElement;
		double subStiffness[3 * 3];

		Element *element = mesh->getEmptyMaterialElement(type);

		for (int elementIndex = 0; elementIndex < mesh->numElements; elementIndex++){
			mesh->updateElementInfo(elementIndex, element);

			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++){
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++){
					element->generateSubStiffnessMatrix(aNodeIndex, bNodeIndex, D, subStiffness);
					//assemble to stiffness matrix
					for (int subRow = 0; subRow < 3; subRow++){
						if (aNodeIndex != bNodeIndex){
							for (int subColumn = 0; subColumn < 3; subColumn++){
								if (subStiffness[3 * subRow + subColumn] != 0.0){
									int i_index = indexer.getGlobalIndex(*element, aNodeIndex, subRow);
									int j_index = indexer.getGlobalIndex(*element, bNodeIndex, subColumn);
									if (i_index >= 0 && j_index >= 0){
										int rowIndex = max(i_index, j_index);
										int columnIndex = min(i_index, j_index);
										matrix.addEntry(rowIndex, columnIndex, subStiffness[3 * subRow + subColumn]);
									}
								}
							}
						}
						else{
							for (int subColumn = 0; subColumn <= subRow; subColumn++){
								if (subStiffness[3 * subRow + subColumn] != 0.0){
									int rowIndex = indexer.getGlobalIndex(*element, aNodeIndex, subRow);
									int columnIndex = indexer.getGlobalIndex(*element, bNodeIndex, subColumn);
									if (rowIndex >= 0 && columnIndex >= 0)
										matrix.addEntry(rowIndex, columnIndex, subStiffness[3 * subRow + subColumn]);
								}
							}
						}
					}
				}
			}
		}

		delete element;
	}

	void StVKMaterial::preprocessWithReduction(const Reference<Mesh> &mesh, const NodeIndexer &indexer, int dofs, const double *basises){
		const int numNodes = mesh->numNodes;
		const int numNodePerElement = mesh->numNodesPerElement;
		double *mem = allocAligned<double>(2 * dofs);
		intergration[0] = mem;
		intergration[1] = mem + dofs;

	}

	void StVKMaterial::getNodeForces(double *ds, double *forces) const{
		for (int i = 0; i < orders; i++){

		}
	}

	StVKMaterial::~StVKMaterial(){
		if (intergration[0])
		    freeAligned(intergration[0]);
		freeAligned(loadParameter);
	}

}


