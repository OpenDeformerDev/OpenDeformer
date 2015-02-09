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

	class RootedLevelStructure{
	public:
		RootedLevelStructure() :root(NULL){};
		RootedLevelStructure(MeshGraphVertexNode *rootNode) :root(rootNode){}
		RootedLevelStructure(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer);
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
		MeshRelabeler(int nodeCount, int elementCount, int nodePerElement, const int *eles);
		MeshRelabeler() = delete;
		MeshRelabeler& operator=(const MeshRelabeler&) = delete;

		void getNewLables(int *newNodeLables, int *newElements);
		~MeshRelabeler(){ 
			freeAligned(elements);
			freeAligned(graphVerts[0]);  
		}
	private:
		template<class Queue, class PriorityQueue> 
		void processSigleVert(MeshGraphVertexNode * vert, std::list<MeshGraphVertexNode *>& levelNodes,
			Queue &working, PriorityQueue &toBeAssigned){
			vert->setNewLable(nextFreeLable++);
			working.push(vert);
			while (!levelNodes.empty()){
				MeshGraphArcNode *node = levelNodes.front()->getAdjacentNodes();
				levelNodes.pop_front();
				while (node){
					MeshGraphVertexNode *vert = graphVerts[node->getLable()];
					if (std::find(levelNodes.begin(), levelNodes.end(), vert) != levelNodes.end()){
						toBeAssigned.push(vert);
						levelNodes.remove(vert);
					}
					node = node->getNextNode();
				}

				while (!toBeAssigned.empty()){
					MeshGraphVertexNode *vert = toBeAssigned.top();
					vert->setNewLable(nextFreeLable++);
					working.push(vert);
					toBeAssigned.pop();
				}
			}
		}

		struct OrderedElement{
			OrderedElement(int nodePerElementCount, const int *unorderedElement, int *mem);
			bool operator<(const OrderedElement& element) const{ return score < element.score; }
			int score;
			int *elementNodeIndices;
		};
		
		int nextFreeLable;
		int nodePerElementCount;
		int numElements;
		int *elements;
		MemoryPool<MeshGraphArcNode> graphArcPool;
		vector<MeshGraphVertexNode *> graphVerts;
	};

}

#endif