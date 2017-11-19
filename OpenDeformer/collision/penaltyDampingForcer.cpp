#include "stdafx.h"
#include "penaltyDampingForcer.h"
#include "intersection.h"

namespace ODER {
	namespace Collision {
		void MeshPlaneCollisionDampingForcer::getDampingMatrix(Scalar *subMat, int *indices) const {
			Assert(elementIndex >= 0);
			const int nodePerElementCount = mesh->getNodePerElementCount();
			nodeIndexer->getElementNodesGlobalIndices(mesh, elementIndex, nodePerElementCount, indices);

			if (nodePerElementCount == 4) {
				constexpr int diagIndices[12] = { 0, 12, 23, 33, 42, 50, 57, 63, 68, 72, 75, 77 };

				DynamicTetrahedronShape tet;
				const int *tetPointIndices = mesh->getElementNodeReference(elementIndex);
				for (int i = 0; i < 4; i++)
					tet.points[i] = mesh->getVertex(tetPointIndices[i]) + mesh->getVertexDisplacementConst(tetPointIndices[i]);

				Vector3 centroid, dir;
				Intersection(plane, tet, centroid, dir);

				Scalar coords[4];
				computeBarycentricCoordinates(tet, centroid, coords);

				Tensor2<Scalar> t = dir ^ dir;

				for (int aNode = 0; aNode < 4; aNode++) {
					for (int i = 0; i < 3; i++) 
						for (int j = 0; j < 3 - i; j++)
							subMat[diagIndices[aNode * 3 + i] + j] = coords[aNode] * coords[aNode] * t(i, j);

					for (int bNode = aNode + 1; bNode < 4; bNode++)
						for (int i = 0; i < 3; i++)
							for (int j = 0; j < 3; j++)
								subMat[diagIndices[aNode * 3 + i] + ((bNode - aNode) * 3 - i) + j] = coords[aNode] * coords[bNode] * t(i, j);
				}
			}
			else 
				Severe("Unimplemented features in MeshPlaneCollisionDampingForcer::getDampingMatrix");

		}
	}
}