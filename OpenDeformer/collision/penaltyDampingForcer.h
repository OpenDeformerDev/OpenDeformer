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
		class PenaltyDampingForcer : public StaDampingForcer {
		public:
			PenaltyDampingForcer(Scalar dampingCoefficient, const Reference<Mesh>& mesh, const Reference<NodeIndexer>& indexer) 
				: dampingCoeff(dampingCoefficient), mesh(mesh), nodeIndexer(indexer) {}
		protected:
			Reference<Mesh> mesh;
			Reference<NodeIndexer> nodeIndexer;
			Scalar dampingCoeff;
		};

		class MeshPlaneCollisionDampingForcer : public PenaltyDampingForcer {
		public:
			MeshPlaneCollisionDampingForcer(Scalar dampingCoefficient, const Reference<Mesh>& mesh, const Reference<NodeIndexer>& indexer,
				const Vector3 *velocities): PenaltyDampingForcer(dampingCoefficient, mesh, indexer), velocities(velocities), elementIndex(-1) {}
			void setCollisionPair(const StaticPlaneShape& plane, int elementIndex);
			void addDampingMatrix(const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& mat) const;
		private:
			int elementIndex;
			StaticPlaneShape plane;
			const Vector3 *velocities;
		};

		inline void MeshPlaneCollisionDampingForcer::setCollisionPair(const StaticPlaneShape& plane, int elementIndex) {
			this->plane = plane;
			this->elementIndex = elementIndex;
		}

		class MeshMeshCollisionDampingForcer : public DynDampingForcer<BlockedSymSpMatrix> {
		public:
			MeshMeshCollisionDampingForcer(Scalar dampingCoefficient, const Reference<Mesh>& aColliderMesh, const Reference<Mesh>& bColliderMesh,
				const Reference<NodeIndexer>& aNodeIndexer, const Reference<NodeIndexer>& bNodeIndexer,
				const Vector3 *velocities) : dampingCoeff(dampingCoefficient), aColliderMesh(aColliderMesh), bColliderMesh(bColliderMesh),
				aNodeIndexer(aNodeIndexer), bNodeIndexer(bNodeIndexer), velocities(velocities),
				aElementIndex(-1), bElementIndex(-1) {}
			void setCollisionPair(int aElementIndex, int bElementIndex);
			void addDampingMatrix(AggSpMatrix<BlockedSymSpMatrix>& mat) const;
		private:
			Scalar dampingCoeff;

			int aElementIndex, bElementIndex;
			Reference<Mesh> aColliderMesh, bColliderMesh;
			Reference<NodeIndexer> aNodeIndexer, bNodeIndexer;

			const Vector3 *velocities;
		};

		inline void MeshMeshCollisionDampingForcer::setCollisionPair(int aElementIndex, int bElementIndex) {
			this->aElementIndex = aElementIndex;
			this->bElementIndex = bElementIndex;
		}

	}
}


#endif // ! ODER_COLLISION_PENALTYDAMPINGFORCER

