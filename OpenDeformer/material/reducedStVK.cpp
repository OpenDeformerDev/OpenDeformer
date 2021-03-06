#include "stdafx.h"
#include "reducedStvk.h"
#include "nodeIndexer.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	ReducedStVKMaterial::ReducedStVKMaterial(Scalar rho, Scalar lameFirst, Scalar lameSecond, int orderNum)
		:ReducedHyperelasticMaterial(rho, orderNum),
		lambda(lameFirst), mu(lameSecond){

		D[0] = lambda + Scalar(2.0) * mu;
		D[1] = lambda;
		D[2] = mu;

		stressNonlinear = NULL;
		intergration[0] = NULL;
		intergration[1] = NULL;

	}

	void ReducedStVKMaterial::generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const{
		const int numNodesPerElement = mesh->getNodePerElementCount();
		Scalar subStiffness[3 * 3];

		ReducedIsotropicElement *element = dynamic_cast<ReducedIsotropicElement *>(mesh->getMaterialElement(type));

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
								if (subStiffness[3 * subRow + subColumn] != 0.0){
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

	void ReducedStVKMaterial::preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer){
		const int numNodes = mesh->getNodeCount();
		const int numNodePerElement = mesh->getNodePerElementCount();
		const int numElements = mesh->getElementCount();

		int commonEntryNum = numNodePerElement*numNodePerElement*numNodePerElement;
		int nlEntries = commonEntryNum * 3;
		int nnEntries = commonEntryNum * numNodePerElement;

		const int entrys = (getNonlinearAsymptoticOrder() - 2)*numNodes*numNodes;
		if (entrys > 0){
			stressNonlinear = new Scalar[entrys];
			memset(stressNonlinear, 0, sizeof(Scalar) * entrys);
		}

		Scalar *mem = allocAligned<Scalar>((nlEntries + nnEntries)*numElements);
		intergration[0] = mem;
		intergration[1] = mem + nlEntries*numElements;

		Scalar *nlpart = allocAligned<Scalar>(nlEntries);
		Scalar *nnpart = allocAligned<Scalar>(nnEntries);

		ReducedIsotropicElement *element = dynamic_cast<ReducedIsotropicElement *>(mesh->getMaterialElement(type));
		for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
			//set new element info
			element->setNodeIndices(elementIndex);
			element->setBMatrixs();
			//get nl and nn
			element->Intergration(&D[1], nlpart, nnpart);
			memcpy(intergration[0] + elementIndex * nlEntries, nlpart, sizeof(Scalar)*nlEntries);
			memcpy(intergration[1] + elementIndex * nnEntries, nnpart, sizeof(Scalar)*nnEntries);
		}

		delete element;
		freeAligned(nlpart);
		freeAligned(nnpart);
	}

	void ReducedStVKMaterial::getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, int totalDofs, const Scalar *ds, Scalar *forces){
		const int numElements = mesh->getElementCount();
		const int numNodesPerElement = mesh->getNodePerElementCount();
		const int numNodes = mesh->getNodeCount();
		const int numNodes2 = numNodes*numNodes;
		memset(forces, 0, totalDofs*sizeof(Scalar));

		GeometricElement *element = mesh->getGeometricElement();
		int commonEntryNum = numNodesPerElement*numNodesPerElement*numNodesPerElement;
		int nlEntries = commonEntryNum * 3;
		int nnEntries = commonEntryNum * numNodesPerElement;

		int *elementNodeIndices = (int *)alloca(3 * numNodesPerElement*sizeof(Scalar));

		int orderCount = getNonlinearAsymptoticOrder();
		bool lowerOrder = (orderCount - order) > 1;

		VectorBase<Scalar> da, db;
		int orderOffset = (order - 1)*numNodes2;
		if (lowerOrder)
			memset(stressNonlinear + orderOffset, 0, sizeof(Scalar) * numNodes2);

		for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
			//set new element info
			element->setNodeIndices(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			const Scalar *nlpart = intergration[0] + elementIndex*nlEntries;
			const Scalar *nnpart = intergration[1] + elementIndex*nnEntries;
			for (int i = 0; i < order; i++){
				const Scalar *di = ds + totalDofs*i;
				const Scalar *dj = ds + totalDofs*(order - i - 1);

				for (int a = 0; a < numNodesPerElement; a++){
					//get node a displacement
					getNodeDisplacements(di, &elementNodeIndices[3 * a], da);

					for (int b = 0; b < numNodesPerElement; b++){
						//get node b displacement
						getNodeDisplacements(dj, &elementNodeIndices[3 * b], db);

						Scalar factor = da * db;
						for (int c = 0; c < numNodesPerElement; c++){
							Scalar factor2 = 0.0;
							for (int axis = 0; axis < 3; axis++){
								//assemble part without stresses
								int index = elementNodeIndices[3 * c + axis];
								if (index >= 0)
									forces[index] -= factor * nlpart[a * 48 + b * 12 + c * 3 + axis] * Scalar(0.5);
								factor2 += da[axis] * nlpart[b * 48 + c * 12 + a * 3 + axis];
							}
							for (int axis = 0; axis < 3; axis++){
								int index = elementNodeIndices[3 * c + axis];
								if (index >= 0)
									forces[index] -= factor2 * db[axis];
							}

							if (lowerOrder){
								int cNodeIndex = element->getNodeIndex(c);
								for (int d = 0; d < numNodesPerElement; d++){
									//assemble part with stresses
									int dNodeIndex = element->getNodeIndex(d);
									stressNonlinear[orderOffset + cNodeIndex * numNodes + dNodeIndex]
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
			const Scalar *d = ds + totalDofs*(order - i - 2);
			int orderOffset = i*numNodes2;
			for (int aNodeIndex = 0; aNodeIndex < numNodes; aNodeIndex++){
				int offset = aNodeIndex*numNodes + orderOffset;
				//init node a indices
				int aNodeGlobalIndices[3];
				bool valid[3];
				aNodeGlobalIndices[0] = indexer->getGlobalIndex(aNodeIndex, 0);
				aNodeGlobalIndices[1] = indexer->getGlobalIndex(aNodeIndex, 1);
				aNodeGlobalIndices[2] = indexer->getGlobalIndex(aNodeIndex, 2);
				valid[0] = aNodeGlobalIndices[0] >= 0;
				valid[1] = aNodeGlobalIndices[1] >= 0;
				valid[2] = aNodeGlobalIndices[2] >= 0;
				if (valid[0]  || valid[1] || valid[2]){
					Scalar entry = 0.0;
					for (int bNodeIndex = 0; bNodeIndex < numNodes; bNodeIndex++){
						for (int bNodeAxis = 0; bNodeAxis < 3; bNodeAxis++){
							int bNodeGlobalIndex = indexer->getGlobalIndex(bNodeIndex, bNodeAxis);
							if (bNodeGlobalIndex >= 0){
								entry += stressNonlinear[offset + bNodeIndex] * d[bNodeGlobalIndex];
							}
						}
					}
					for (int axis = 0; axis < 3; axis++){
						if (valid[axis])
							forces[aNodeGlobalIndices[axis]] -= entry;
					}
				}
			}
		}

		delete element;
	}

	void ReducedStVKMaterial::getNodeDisplacements(const Scalar *ds, const int *nodeIndices, VectorBase<Scalar>& d) const{
		for (int axis = 0; axis < 3; axis++){
			int index = nodeIndices[axis];
			d[axis] = index >= 0 ? ds[index] : 0;
		}
	}

	ReducedStVKMaterial::~ReducedStVKMaterial(){
		if (intergration[0]) freeAligned(intergration[0]);
		delete[] stressNonlinear;
	}
}


