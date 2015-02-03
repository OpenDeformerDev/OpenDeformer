#include "stdafx.h"
#include "meshRelabeler.h"
#include "datastructure.h"

namespace ODER{
	void MeshGraphVertexNode::insertArc(int lable, MemoryPool<MeshGraphArcNode>& pool){
		MeshGraphArcNode *node = pool.Alloc();
		node->setLable(lable);
		node->setNextNode(adjacent);
		adjacent = node;
		degree++;
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
			}
			if (now){
				pre->setNextNode(now->getNextNode());
				degree--;
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
		generateLevels(rootNode, graphVerts, visitedBuffer);
	}

	int LevelStructure::getWidth() const{
		int levelSize = getLevelSize();
		int width = INT_MIN;
		for (int i = 0; i < levelSize; i++){
			int levelWidth = levels[i].size();
			if (levelWidth > width)
				width = levelWidth;
		}

		return width;
	}

	void LevelStructure::generateLevels(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer){
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
		int vertSize = graphVerts.size();
		bool *visited = new bool[vertSize];
		memset(visited, 0, sizeof(bool)*vertSize);
		MeshGraphVertexNode *start = *std::min_element(graphVerts.begin(), graphVerts.end(),
			[](const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){return left->getDegree() < right->getDegree(); });

		bool find = false;
		LevelStructure startLevels(start, graphVerts, visited);
		LevelStructure endLevels, mayEndLevels;
		while (!find){
			find = true;
			int endWidth = INT_MAX;
			int lastLevel = startLevels.getLevelSize() - 1;
			startLevels.sortLevel(lastLevel);
			auto lastLevelNode = startLevels.getLevelNodeIterBeign(lastLevel);
			auto lastLevelEnd = startLevels.getLevelNodeIterEnd(lastLevel);
			while (lastLevelNode != lastLevelEnd){
				memset(visited, 0, sizeof(bool)*vertSize);
				mayEndLevels.generateLevels(*lastLevelNode, graphVerts, visited);
				int levelSize = mayEndLevels.getLevelSize();
				if (levelSize > startLevels.getLevelSize()){
					find = false;
					startLevels = mayEndLevels;
					break;
				}
				int width = mayEndLevels.getWidth();
				if (width < endWidth){
					endLevels = mayEndLevels;
					endWidth = width;
				}
				lastLevelNode++;
			}
		}

		//minimizing level width
		int startLevelSize = startLevels.getLevelSize(), endLevelSize = endLevels.getLevelSize();
		int size = max(startLevelSize, endLevelSize);

		struct OrderedPair{
			OrderedPair(int start = INT_MIN, int end = INT_MIN) :startLevel(start), endLevel(end){}
			int startLevel;
			int endLevel;
		};

		struct Comparer{
			bool operator()(const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){
				return size_t(left) < size_t(right);
			}
		};

		map<MeshGraphVertexNode *, OrderedPair, Comparer> pairs;

		pairs.insert(std::pair<MeshGraphVertexNode *, OrderedPair>(startLevels.getRoot(), OrderedPair(0)));
		for (int i = 1; i < startLevelSize + 1; i++){
			auto levelNode = startLevels.getLevelNodeIterBeign(i);
			auto levelEnd = startLevels.getLevelNodeIterEnd(i);
			while (levelNode != levelEnd)
				pairs.insert(std::pair<MeshGraphVertexNode *, OrderedPair>(*levelNode++, OrderedPair(i + 1)));
		}

		pairs[endLevels.getRoot()].endLevel = endLevelSize;
		for (int i = 1; i < endLevelSize + 1; i++){
			auto levelNode = endLevels.getLevelNodeIterBeign(i);
			auto levelEnd = endLevels.getLevelNodeIterEnd(i);
			while (levelNode != levelEnd)
				pairs[*levelNode++].endLevel = endLevelSize - i;
		}

		//delete arcs
		vector<MeshGraphVertexNode *> verts = graphVerts;
		for (int i = 0; i < verts.size(); i++){
			MeshGraphVertexNode *vert = graphVerts[i];
			auto pair = pairs[vert];
			if (pair.startLevel == pair.endLevel){
				verts[i] = NULL;
				int label = vert->getOldLable();
				MeshGraphArcNode* arc = vert->getAdjacentNodes();
				while (arc){
					graphVerts[arc->getLable()]->deleteArc(label);
					arc = arc->getNextNode();
				}
			}
		}

		DisjointSets<MeshGraphVertexNode *, Comparer> connetComponents;
		for (auto vert : verts){
			if (vert){
				connetComponents.makeSet(vert);
			}
		}

	}

}