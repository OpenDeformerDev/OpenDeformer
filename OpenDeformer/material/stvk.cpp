#include "stdafx.h"
#include "stvk.h"
#include "nodeIndexer.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	StVKMaterial::StVKMaterial(double rho, double lameFirst, double lameSecond, int orderNum)
		: HyperelasticMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear), orderNum),
		lambda(lameFirst), mu(lameSecond){

		D[0] = lambda + 2.0 * mu;
		D[1] = lambda;
		D[2] = mu;

		stressNonlinear = NULL;
		intergration[0] = NULL;
		intergration[1] = NULL;
	}

	void StVKMaterial::generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const{
		const int numNodesPerElement = mesh->getNodePerElementCount();
		double subStiffness[3 * 3];

		Element *element = mesh->getEmptyMaterialElement(type);

		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++){
			element->setNodeIndexs(elementIndex);
			element->setBMatrixs();

			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++){
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++){
					element->generateSubStiffnessMatrix(aNodeIndex, bNodeIndex, D, subStiffness);
					//assemble to stiffness matrix
					for (int subRow = 0; subRow < 3; subRow++){
						if (aNodeIndex != bNodeIndex){
							for (int subColumn = 0; subColumn < 3; subColumn++){
								if (subStiffness[3 * subRow + subColumn] != 0.0){
									int i_index = indexer->getGlobalIndex(*element, aNodeIndex, subRow);
									int j_index = indexer->getGlobalIndex(*element, bNodeIndex, subColumn);
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

	void StVKMaterial::preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer){
		const int numNodes = mesh->getNodeCount();
		const int numNodePerElement = mesh->getNodePerElementCount();
		const int numElements = mesh->getElementCount();

		int commonEntryNum = numNodePerElement*numNodePerElement*numNodePerElement;
		int nlEntries = commonEntryNum * 3;
		int nnEntries = commonEntryNum * numNodePerElement;

		const int entrys = (getNonlinearAsymptoticOrder() - 2)*numElements*numElements;
		if (entrys > 0){
			stressNonlinear = new double[entrys];
			memset(stressNonlinear, 0, sizeof(double) * entrys);
		}

		double *mem = allocAligned<double>((nlEntries + nnEntries)*numElements);
		memset(mem, 0, sizeof(double) * (nlEntries + nnEntries)*numElements);
		intergration[0] = mem;
		intergration[1] = mem + nlEntries*numElements;

		double *nlpart = allocAligned<double>(nlEntries);
		double *nnpart = allocAligned<double>(nnEntries);

		Element *element = mesh->getEmptyMaterialElement(type);
		for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
			//set new element info
			element->setNodeIndexs(elementIndex);
			element->setBMatrixs();

			//get nl and nn
			element->Intergration(&D[1], nlpart, nnpart);
			memcpy(intergration[0] + elementIndex * nlEntries, nlpart, sizeof(double)*nlEntries);
			memcpy(intergration[1] + elementIndex * nnEntries, nnpart, sizeof(double)*nnEntries);
		}

		delete element;
		freeAligned(nlpart);
		freeAligned(nnpart);
	}

	void StVKMaterial::getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, int totalDofs, const double *ds, double *forces){
		const int numElements = mesh->getElementCount();
		const int numNodesPerElement = mesh->getNodePerElementCount();
		const int numElements2 = numElements*numElements;

		Element *element = mesh->getEmptyElement();
		int commonEntryNum = numNodesPerElement*numNodesPerElement*numNodesPerElement;
		int nlEntries = commonEntryNum * 3;
		int nnEntries = commonEntryNum * numNodesPerElement;

		int *elementNodeIndices = new int[3 * numNodesPerElement];

		bool lowerOrder = (getNonlinearAsymptoticOrder() - order) > 1;
		memset(forces, 0, totalDofs*sizeof(double));
		VectorBase<double> da, db;
		for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
			//set new element info
			element->setNodeIndexs(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			const double *nlpart = intergration[0] + elementIndex*nlEntries;
			const double *nnpart = intergration[1] + elementIndex*nnEntries;
			for (int i = 0; i < order; i++){
				const double *di = ds + totalDofs*i;
				const double *dj = ds + totalDofs*(order - i - 1);

				int orderOffset = (order - 1)*numElements2;
				for (int a = 0; a < numNodesPerElement; a++){
					//get node a displacement
					getNodeDisplacements(di, &elementNodeIndices[numNodesPerElement*a], da);

					for (int b = 0; b < numNodesPerElement; b++){
						//get node b displacement
						getNodeDisplacements(dj, &elementNodeIndices[numNodesPerElement*b], db);

						double factor = da*db;
						for (int c = 0; c < numNodesPerElement; c++){
							double factor2 = 0.0;
							for (int axis = 0; axis < 3; axis++){
								//assemble part without stresses
								int index = elementNodeIndices[numNodesPerElement*c + axis];
								if (index >= 0){
									forces[index] -= factor * nlpart[a * 48 + b * 12 + c * 3 + axis] * 0.5;
									factor2 += da[axis] * nlpart[b * 48 + c * 12 + a * 3 + axis];
								}
							}
							for (int axis = 0; axis < 3; axis++){
								int index = elementNodeIndices[numNodesPerElement*c + axis];
								if (index >= 0){
									forces[index] -= factor2*db[axis];
								}
							}

							if (lowerOrder){
								int cNodeIndex = element->getNodeIndex(c);
								for (int d = 0; d < numNodesPerElement; d++){
									//assemble part with stresses
									int dNodeIndex = element->getNodeIndex(d);
									stressNonlinear[orderOffset + cNodeIndex * numElements + dNodeIndex]
										+= factor*nnpart[a * 64 + b * 16 + c * 4 + d];
								}
							}
						}
					}
				}
			}
		}

		//assmble lower order nonlinear part
		for (int i = 0; i < order - 1; i++){
			const double *d = ds + totalDofs*(order - i - 2);
			int orderOffset = i*numElements;
			for (int aNodeIndex = 0; aNodeIndex < mesh->getNodeCount(); aNodeIndex++){
				int offset = aNodeIndex*numElements + orderOffset;
				for (int bNodeIndex = 0; bNodeIndex < mesh->getNodeCount(); bNodeIndex++){
					for (int axis = 0; axis < 3; axis++){
						int index = indexer->getGlobalIndex(bNodeIndex, axis);
						if (index >= 0){
							forces[index] -= stressNonlinear[offset + bNodeIndex] * d[index];
						}
					}
				}
			}
		}

		delete[] elementNodeIndices;
		delete element;
	}

	void StVKMaterial::getNodeDisplacements(const double *ds, const int *nodeIndices, VectorBase<double>& d) const{
		for (int axis = 0; axis < 3; axis++){
			int index = nodeIndices[axis];
			if (index >= 0) d[axis] = ds[index];
			else d[axis] = 0.0;
		}
	}

	StVKMaterial::~StVKMaterial(){
		if (intergration[0]) freeAligned(intergration[0]);
		if (stressNonlinear) delete[] stressNonlinear;
	}
}


