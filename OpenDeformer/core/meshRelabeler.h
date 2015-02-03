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
		void deleteArc(int lable);
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

	class LevelStructure{
	public:
		LevelStructure() :root(NULL){};
		LevelStructure(MeshGraphVertexNode *rootNode) :root(rootNode){}
		LevelStructure(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer);
		void setRoot(MeshGraphVertexNode *node){ root = node; }
		MeshGraphVertexNode *getRoot(){ return root; }
		void generateLevels(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer);
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
			std::sort(levels[level - 1].begin(), levels[level - 1].end(),
				[](const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){return left->getDegree() < right->getDegree(); });
		}
		int getLevelSize() const { return levels.size();}
		int getWidth() const;
	private:
		MeshGraphVertexNode *root;
		vector<vector<MeshGraphVertexNode *>> levels;
	};

	class MeshRelabeler{
	public:
		MeshRelabeler(int nodeCount, int edgeCount, const int *edgeNodeIndices);
		MeshRelabeler() = delete;
		MeshRelabeler& operator=(const MeshRelabeler&) = delete;

		void getNewLables(int *newNodeLables) const;
		~MeshRelabeler(){ delete[] graphVerts[0]; }
	private:
		int nextFreeLable;
		MemoryPool<MeshGraphArcNode> graphArcPool;
		vector<MeshGraphVertexNode *> graphVerts;
	};

}

#endif