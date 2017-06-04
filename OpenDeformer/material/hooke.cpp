#include "stdafx.h"
#include "hooke.h"
#include "nodeIndexer.h"

namespace ODER{
	HookeMaterial::HookeMaterial(Scalar rho, Scalar YOUNGS, Scalar POSSION)
		:MechMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_Linear)),
		youngs(YOUNGS), possion(POSSION){
		//Lame prameters
		Scalar lambda = possion*youngs / ((Scalar(1.0) + possion) * (Scalar(1.0) - Scalar(2.0) * possion));
		Scalar mu = youngs / (Scalar(2.0) * (Scalar(1.0) + possion));

		D[0] = lambda + Scalar(2.0) * mu;
		D[1] = lambda;
		D[2] = mu;
	}

	void HookeMaterial::generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const{
		const int numNodesPerElement = mesh->getNodePerElementCount();
		Scalar subStiffness[3 * 3];

		LinearIsotropicElement *element = dynamic_cast<LinearIsotropicElement *>(mesh->getMaterialElement(type));

		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++){
			element->setNodeIndices(elementIndex);
			element->setBMatrixs();

			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++){
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++){
					element->generateSubStiffnessMatrix(aNodeIndex, bNodeIndex, D, subStiffness);
					//assemble to stiffness matrix
					for (int subRow = 0; subRow < 3; subRow++){
						if (aNodeIndex != bNodeIndex){
							for (int subColumn = 0; subColumn < 3; subColumn++){
								if (subStiffness[3 * subRow + subColumn] != Scalar(0)){
									int i_index = indexer->getGlobalIndex(*element, aNodeIndex, subRow);
									int j_index = indexer->getGlobalIndex(*element, bNodeIndex, subColumn);
									if (i_index >= 0 && j_index >= 0){
										int rowIndex = std::max(i_index, j_index);
										int columnIndex = std::min(i_index, j_index);
										matrix.addEntry(rowIndex, columnIndex, subStiffness[3 * subRow + subColumn]);
									}
								}
							}
						}
						else{
							for (int subColumn = 0; subColumn <= subRow; subColumn++){
								if (subStiffness[3 * subRow + subColumn] != Scalar(0)){
									int rowIndex = indexer->getGlobalIndex(*element, aNodeIndex, subRow);
									int columnIndex = indexer->getGlobalIndex(*element, bNodeIndex, subColumn);
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
}