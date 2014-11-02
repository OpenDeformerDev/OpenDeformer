#include "stdafx.h"
#include "stvk.h"
#include "nodeIndexer.h"

namespace ODER{
	StVKMaterial::StVKMaterial(double rho, double lameFirst, double lameSecond, int numElements, int orderNum)
		:MecMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear)),
		lambda(lameFirst), mu(lameSecond), orders(orderNum){

		D[0] = lambda + 2.0 * mu;
		D[1] = lambda;
		D[2] = mu;

		stressNonlinear = allocAligned<double **>(orders);
		double **pointerMem = allocAligned<double *>(orders*numElements);
		double *mem = allocAligned<double>(orders*numElements*numElements);
		memset(mem, 0, orders*numElements*numElements);
		for (int i = 0; i < orders; i++){
			stressNonlinear[i] = pointerMem;
			pointerMem += numElements;
			for (int j = 0; j < numElements; j++){
				stressNonlinear[i][j] = mem;
				mem += numElements;
			}
		}

		intergration[0] = NULL;
		intergration[1] = NULL;
	}

	void StVKMaterial::generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const{
		const int numNodesPerElement = mesh->numNodesPerElement;
		double subStiffness[3 * 3];

		Element *element = mesh->getEmptyMaterialElement(type);

		for (int elementIndex = 0; elementIndex < mesh->numElements; elementIndex++){
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
		const int numNodes = mesh->numNodes;
		const int numNodePerElement = mesh->numNodesPerElement;
		const int numElements = mesh->numElements;
		double **pointerMem = allocAligned<double *>(2 * numElements);
		intergration[0] = pointerMem;
		intergration[1] = pointerMem + numElements;

		int commonEntryNum = numNodePerElement*numNodePerElement*numNodePerElement;
		int nlEntries = commonEntryNum * 3;
		int nnEntries = commonEntryNum * numNodePerElement;

		double *mem = allocAligned<double>((nlEntries + nnEntries)*numElements);
		for (int i = 0; i < numElements; i++){
			intergration[0][i] = mem;
			mem += nlEntries;
		}
		for (int i = 0; i < numElements; i++){
			intergration[1][i] = mem;
			mem += nnEntries;
		}


		double *nlpart = allocAligned<double>(nlEntries);
		double *nnpart = allocAligned<double>(nnEntries);

		int totalDofs = 3 * indexer->getMatrixOrder(mesh);
		double *nl = new double[totalDofs];
		double *nn = new double[totalDofs];
		memset(nl, 0, totalDofs*sizeof(double));
		memset(nn, 0, totalDofs*sizeof(double));

		Element *element = mesh->getEmptyMaterialElement(type);
		for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
			//set new element info
			element->setNodeIndexs(elementIndex);
			element->setBMatrixs();

			//get nl and nn
			element->Intergration(&D[1], nlpart, nnpart);
			memcpy(intergration[0][elementIndex], nlpart, nlEntries);
			memcpy(intergration[1][elementIndex], nnpart, nnEntries);
		}

		delete element;
		delete[] nl;
		delete[] nn;
		freeAligned(nlpart);
		freeAligned(nnpart);
	}

	void StVKMaterial::getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, const double *ds, double *forces) {
		const int numElements = mesh->numElements;
		const int numNodesPerElement = mesh->numNodesPerElement;

		Element *element = mesh->getEmptyMaterialElement(type);
		int commonEntryNum = numNodesPerElement*numNodesPerElement*numNodesPerElement;
		int nlEntries = commonEntryNum * 3;
		int nnEntries = commonEntryNum * numNodesPerElement;

		double *nlpart = NULL;
		double *nnpart = NULL;

		int *elementNodeIndices = (int *)alloca(numNodesPerElement*sizeof(int));

		memset(forces, 0, indexer->getMatrixOrder(mesh)*sizeof(double));
		VectorBase<double> da, db, dc, dd;
		for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
			//set new element info
			element->setNodeIndexs(elementIndex);
			element->setBMatrixs();
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			nlpart = intergration[0][elementIndex*nlEntries];
			nnpart = intergration[1][elementIndex*nnEntries];
			for (int i = 0; i < order; i++){
				const double *di = ds + numElements*i;
				const double *dj = ds + numElements*(order - i - 1);

				for (int a = 0; a < numNodesPerElement; a++){
					//get node a displacement
					getNodeDisplacements(di, &elementNodeIndices[numNodesPerElement*a], da);

					for (int b = 0; b < numNodesPerElement; b++){
						//get node b displacement
						getNodeDisplacements(dj, &elementNodeIndices[numNodesPerElement*b], db);

						double factor = da*db;
						for (int c = 0; c < numNodesPerElement; c++){
							//get c node displacement
							getNodeDisplacements(dj, &elementNodeIndices[numNodesPerElement*c], dc);

							double factor2 = 0.0;
							for (int axis = 0; axis < 3; axis++){
								//assemble part without stresses
								int index = elementNodeIndices[numNodesPerElement*c + axis];
								if (index >= 0){
									forces[index] += factor*nlpart[a * 48 + b * 12 + c * 3 + axis];
									factor2 += da[axis] * nlpart[a * 48 + b * 12 + c * 3 + axis];
								}
							}
							for (int axis = 0; axis < 3; axis++){
								int index = elementNodeIndices[numNodesPerElement*c + axis];
								if (index >= 0){
									forces[index] += factor2*db[axis];
								}
							}

							int cNodeIndex = element->nodeIndexs[c];
							for (int d = 0; d < numNodesPerElement; d++){
								//assemble part with stresses
								int dNodeIndex = element->nodeIndexs[d];
								stressNonlinear[order - 1][dNodeIndex][cNodeIndex] += factor*nnpart[a * 64 + b * 16 * c * 4 + d];
							}
						}
					}
				}
				//assmble lower order nonlinear part
				if (i < order - 1){
					for (int aNodeIndex = 0; aNodeIndex < mesh->numNodes; aNodeIndex++){
						for (int bNodeIndex = 0; bNodeIndex < mesh->numNodes; bNodeIndex++){
							for (int axis = 0; axis < 3; axis++){
								int index = indexer->getGlobalIndex(bNodeIndex, axis);
								if (index >= 0){
									forces[index] += stressNonlinear[i][aNodeIndex][bNodeIndex] * dj[index];
								}
							}
						}
					}
				}
			}
		}
		delete element;
	}

	void StVKMaterial::getNodeDisplacements(const double *ds, const int *nodeIndices, VectorBase<double>& d) const{
		for (int axis = 0; axis < 3; axis++){
			int index = nodeIndices[axis];
			if (index >= 0)
				d[index] = ds[index];
			else
				d[index] = 0.0;
		}
	}

	StVKMaterial::~StVKMaterial(){
		freeAligned(stressNonlinear[0][0]);
		freeAligned(stressNonlinear[0]);
		freeAligned(stressNonlinear);
		if (intergration[0])
			freeAligned(intergration[0][0]);
		    freeAligned(intergration[0]);
	}
}


