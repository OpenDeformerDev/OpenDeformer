#include "stdafx.h"
#include "meshRelabeler.h"

namespace ODER{
	void MeshGraphVertexNode::insertArc(int lable, MemoryPool<MeshGraphArcNode>& pool){
		degree++;
		MeshGraphArcNode *node = pool.Alloc();
		node->setLable(lable);
		node->setNextNode(adjacent);
		adjacent = node;
	}

	void MeshGraphVertexNode::deleteArc(int lable){
		if (adjacent->getLable() != lable){
			MeshGraphArcNode *pre = adjacent;
			MeshGraphArcNode *now = pre->getNextNode();
			while (now && now->getLable() != lable){
				pre = now;
				now = now->getNextNode();
			}
			if (now){
				pre->setNextNode(now->getNextNode());
				degree--;
			}
		}
		else{
			adjacent = adjacent->getNextNode();
			degree--;
		}
	}

	void MeshGraphVertexNode::deleteArc(int lable, MemoryPool<MeshGraphArcNode>& pool){
		if (adjacent->getLable() != lable){
			MeshGraphArcNode *pre = adjacent;
			MeshGraphArcNode *now = pre->getNextNode();
			while (now && now->getLable() != lable){
				pre = now;
				now = now->getNextNode();
				degree--;
			}
			if (now){
				pre->setNextNode(now->getNextNode());
				pool.Dealloc(now);
			}
		}
		else{
			MeshGraphArcNode *deleted = adjacent;
			adjacent = adjacent->getNextNode();
			degree--;
			pool.Dealloc(deleted);
		}
	}

	MeshRelabeler::MeshRelabeler(int nodeCount, int edgeCount, const int *edgeNodeIndices){
		nextFreeLable = 0;
		graphVerts.reserve(nodeCount);
		MeshGraphVertexNode *mem = new MeshGraphVertexNode[nodeCount];
		for (int i = 0; i < nodeCount; i++){
			graphVerts[i] = mem + i;
			graphVerts[i]->setOldLable(i);
		}

		for (int i = 0; i < edgeCount; i++){
			int a = edgeNodeIndices[2 * i];
			int b = edgeNodeIndices[2 * i + 1];
			graphVerts[a]->insertArc(b, graphArcPool);
			graphVerts[b]->insertArc(a, graphArcPool);
		}
	}
	LevelStructure::LevelStructure(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer){
		root = rootNode;
		visitedBuffer[root->getOldLable()] = true;
		MeshGraphArcNode *node = root->getAdjacentNodes();
		std::queue<MeshGraphVertexNode *> lableStarters;
		//level 1 construction
		while (node){
			int lable = node->getLable();
			MeshGraphVertexNode *vertNode = graphVerts[lable];
			levels[0].push_back(vertNode);
			lableStarters.push(vertNode);
			visitedBuffer[lable] = true;
			node = node->getNextNode();
		}
		//other levels
		int levelNum = 1;
		while (!lableStarters.empty()){
			int unLebledSize = lableStarters.size();
			for (int i = 0; i < unLebledSize; i++){
				node = lableStarters.front()->getAdjacentNodes();
				lableStarters.pop();
				while (node){
					int lable = node->getLable();
					if (!visitedBuffer[lable]){
						MeshGraphVertexNode *vertNode = graphVerts[lable];
						levels[levelNum].push_back(vertNode);
						lableStarters.push(vertNode);
						visitedBuffer[lable] = true;
					}
					node = node->getNextNode();
				}
			}
			levelNum++;
		}
	}

	void MeshRelabeler::getNewLables(int *newNodeLables) const{
		//find endpoints of pseudo-diameter
		int size = graphVerts.size();
		bool *visited = new bool[size];
		memset(visited, 0, sizeof(bool)*size);
		


	}

}