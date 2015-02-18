#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MESHRELABELER_H
#define ODER_CORE_MESHRELABELER_H

#include "oder.h"
#include "latool.h"
#include "memory.h"

namespace ODER{
	class MeshGraphArcNode{
	public:
		MeshGraphArcNode() :lable(-1), next(NULL){};
		int getLable() const { return lable; }
		void setLable(int LABLE){ lable = LABLE; }
		MeshGraphArcNode* getNextNode() const{ return next; }
		void setNextNode(MeshGraphArcNode *node){ next = node; }
	private:
		int lable;
		MeshGraphArcNode* next;
	};

	class MeshGraphVertexNode{
	public:
		MeshGraphVertexNode() :newLable(-1), oldLable(-1), degree(0), adjacent(NULL){};
		void insertArc(int lable, MemoryPool<MeshGraphArcNode>& pool);
		void deleteArc(int lable, MemoryPool<MeshGraphArcNode>& pool);
		int getDegree() const { return degree; }
		void setNewLable(int lable){ newLable = lable; }
		int getNewLable() const{ return newLable; }
		void setOldLable(int lable){ oldLable = lable; }
		int getOldLable() const{ return oldLable; }
		MeshGraphArcNode* getAdjacentNodes() const{ return adjacent; }
	private:
		int newLable;
		int oldLable;
		int degree;
		MeshGraphArcNode* adjacent;
	};

	struct MeshGraphEdge{
		MeshGraphEdge(int v0, int v1){
			if (v0 > v1) std::swap(v0, v1);
			head = v0; tail = v1;
		}
		bool operator<(const MeshGraphEdge& edge) const{
			if (head == edge.head)
				return tail < edge.tail;
			return head < edge.head;
		}
		int head;
		int tail;
	};

	struct DegreeComparer{
		bool operator()(const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){ 
			return left->getDegree() < right->getDegree(); 
		}
	};

	class RootedLevelStructure{
	public:
		RootedLevelStructure() :root(NULL){}
		RootedLevelStructure(MeshGraphVertexNode *rootNode) :root(rootNode){}
		RootedLevelStructure(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer);
		RootedLevelStructure(const RootedLevelStructure&) = default;
		RootedLevelStructure& operator=(const RootedLevelStructure&) = default;
		RootedLevelStructure(RootedLevelStructure&& levelstructure) :levels(std::move(levelstructure.levels)){
			root = levelstructure.root;
			levelstructure.root = NULL;
		}
		RootedLevelStructure& operator=(RootedLevelStructure&& levelstructure){
			std::swap(root, levelstructure.root);
			levels = std::move(levelstructure.levels);
			return *this;
		}
		void setRoot(MeshGraphVertexNode *node){ root = node; }
		MeshGraphVertexNode *getRoot(){ return root; }
		void generateLevels(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer);
		void Clear();
		void insertVertex(int levelNum, MeshGraphVertexNode *node){
			Assert(levelNum > 0);
			levels[levelNum - 1].push_back(node);
		}
		vector<MeshGraphVertexNode *>::iterator getLevelNodeIterBeign(int levelNum){
			Assert(levelNum > 0);
			return levels[levelNum - 1].begin();
		}
		vector<MeshGraphVertexNode *>::iterator getLevelNodeIterEnd(int levelNum){
			Assert(levelNum > 0);
			return levels[levelNum - 1].end();
		}
		void sortLevel(int level){
			Assert(level > 0);
			std::sort(levels[level - 1].begin(), levels[level - 1].end(),
				[](const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){return left->getDegree() < right->getDegree(); });
		}
		int getLevelSize() const { return levels.size() + 1;}
		int getWidth() const;
		~RootedLevelStructure() = default;
	private:
		MeshGraphVertexNode *root;
		vector<vector<MeshGraphVertexNode *>> levels;
	};

	class MeshRelabeler{
	public:
		MeshRelabeler(int nodeCount);
		MeshRelabeler(const Mesh& mesh);
		MeshRelabeler() = delete;
		MeshRelabeler& operator=(const MeshRelabeler&) = delete;

		void getNewLables(int elementCount, int nodePerElementCount, int *newNodeLables, int *elements);
		void getNewLables(int *newNodeLables, Mesh& mesh);
		~MeshRelabeler(){ 
			freeAligned(graphVerts[0]);  
		}
	private:
		void setNewNodeLables(int *newNodeLables);
		void processSigleVert(MeshGraphVertexNode * vert, std::list<MeshGraphVertexNode *>& levelNodes,
			std::queue<MeshGraphVertexNode *>& working, std::priority_queue<MeshGraphVertexNode *, vector<MeshGraphVertexNode *>, DegreeComparer> &toBeAssigned);
		void generateGraphFromElementIndices(int elementCount, int nodePerElementCount, int *elements);
		void generateGraphFromMesh(const Mesh& mesh);
		struct OrderedElement{
			OrderedElement() :score(-1), elementNodeIndices(NULL){}
			OrderedElement(int nodePerElementCount, const int *indices, int *mem);
			bool operator<(const OrderedElement& element) const{ return score < element.score; }
			int score;
			int *elementNodeIndices;
		};

		int nextFreeLable;
		MemoryPool<MeshGraphArcNode> graphArcPool;
		vector<MeshGraphVertexNode *> graphVerts;
	};
}

#endif