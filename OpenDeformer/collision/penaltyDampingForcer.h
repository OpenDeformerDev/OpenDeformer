#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_COLLISION_PENALTYDAMPINGFORCER
#define ODER_COLLISION_PENALTYDAMPINGFORCER

#include "oder.h"
#include "forcer.h"
#include "nodeIndexer.h"
#include "collisionShapes.h"

namespace ODER {
	namespace Collision {
		class PenaltyDampingForcer : public DampingForcer {
		public:
			PenaltyDampingForcer(Scalar dampingCoefficient, const Reference<Mesh>& mesh, const Reference<NodeIndexer>& indexer) 
				: dampingCoeff(dampingCoefficient), mesh(mesh), nodeIndexer(indexer) {}
		protected:
			Reference<Mesh> mesh;
			Reference<NodeIndexer> nodeIndexer;
			Scalar dampingCoeff;
		};

		class MeshPlaneCollisionDampingForcer : public PenaltyDampingForcer {
			MeshPlaneCollisionDampingForcer(Scalar dampingCoefficient, const Reference<Mesh>& mesh, const Reference<NodeIndexer>& indexer)
				: PenaltyDampingForcer(dampingCoefficient, mesh, indexer), elementIndex(-1) {}
			void setCollisionPair(const StaticPlaneShape& plane, int elementIndex);
			void addDampingMatrix(const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& mat) const;
		private:
			int elementIndex;
			StaticPlaneShape plane;
		};

		inline void MeshPlaneCollisionDampingForcer::setCollisionPair(const StaticPlaneShape& plane, int elementIndex) {
			this->plane = plane;
			this->elementIndex = elementIndex;
		}

	}
}


#endif // ! ODER_COLLISION_PENALTYDAMPINGFORCER

