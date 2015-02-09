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

	RootedLevelStructure::RootedLevelStructure(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer){
		generateLevels(rootNode, graphVerts, visitedBuffer);
	}

	int RootedLevelStructure::getWidth() const{
		int levelSize = getLevelSize();
		int width = INT_MIN;
		for (int i = 0; i < levelSize; i++){
			int levelWidth = levels[i].size();
			if (levelWidth > width)
				width = levelWidth;
		}

		return width;
	}

	void RootedLevelStructure::generateLevels(MeshGraphVertexNode *rootNode, const vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer){
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

	MeshRelabeler::OrderedElement::OrderedElement(int nodePerElementCount, const int *unorderedElement, int *mem){
		memcpy(elementNodeIndices, unorderedElement, sizeof(int)*nodePerElementCount);
		score = 0;
		for (int i = 0; i < nodePerElementCount; i++)
			score += unorderedElement[i];
	}

	MeshRelabeler::MeshRelabeler(int nodeCount, int elementCount, int nodePerElement, const int *eles){
		nextFreeLable = 0;
		nodePerElementCount = nodePerElement;
		numElements = elementCount;
		graphVerts.reserve(nodeCount);
		elements = allocAligned<int>(elementCount*nodePerElement);

		MeshGraphVertexNode *mem = allocAligned<MeshGraphVertexNode>(nodeCount);
		for (int i = 0; i < nodeCount; i++){
			graphVerts[i] = mem + i;
			graphVerts[i]->setOldLable(i);
		}
		memcpy(elements, eles, sizeof(int)*elementCount*nodePerElement);

		struct Edge{
			Edge(int v0, int v1){
				if (v0 > v1) std::swap(v0, v1);
				head = v0; tail = v1;
			}
			bool operator<(const Edge& edge) const{
				if (head == edge.head)
					return tail < edge.tail;
				return head < edge.head;
			}
			int head;
			int tail;
		};

		set<Edge> edges;
		const int* element = elements;
		for (int i = 0; i < elementCount; i++){
			for (int j = 0; j < nodePerElement; j++){
				for (int k = j + 1; k < nodePerElement; k++)
					edges.insert(Edge(element[j], element[k]));
			}
			element += elementCount;
		}

		for (auto edge : edges){
			graphVerts[edge.head]->insertArc(edge.tail, graphArcPool);
			graphVerts[edge.tail]->insertArc(edge.head, graphArcPool);
		}
	}

	void MeshRelabeler::getNewLables(int *newNodeLables, int *newElements){
		//find endpoints of pseudo-diameter
		int vertSize = graphVerts.size();
		bool *visited = new bool[vertSize];
		memset(visited, 0, sizeof(bool)*vertSize);
		MeshGraphVertexNode *start = *std::min_element(graphVerts.begin(), graphVerts.end(),
			[](const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){return left->getDegree() < right->getDegree(); });

		bool find = false;
		RootedLevelStructure startLevels(start, graphVerts, visited);
		RootedLevelStructure endLevels, mayEndLevels;
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
		int size = max(startLevelSize, endLevelSize) + 1;
		int startSmaller = endLevels.getWidth() > startLevels.getWidth();

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
		vector<MeshGraphVertexNode *> deletedVerts;
		vector<std::list<MeshGraphVertexNode *>> newLevels;
		newLevels.reserve(size);

		pairs.insert(std::pair<MeshGraphVertexNode *, OrderedPair>(startLevels.getRoot(), OrderedPair(0)));
		for (int i = 1; i < startLevelSize + 1; i++){
			auto levelNode = startLevels.getLevelNodeIterBeign(i);
			auto levelEnd = startLevels.getLevelNodeIterEnd(i);
			while (levelNode != levelEnd)
				pairs.insert(std::pair<MeshGraphVertexNode *, OrderedPair>(*levelNode++, OrderedPair(i)));
		}

		MeshGraphVertexNode *vert = endLevels.getRoot();
		auto pair = pairs[vert];
		int end = endLevelSize;
		pair.endLevel = end;
		if (pair.startLevel == end){
			deletedVerts.push_back(vert);
			newLevels[end].push_back(vert);
		}
		for (int i = 1; i < endLevelSize + 1; i++){
			end = endLevelSize - i;
			auto levelNode = endLevels.getLevelNodeIterBeign(i);
			auto levelEnd = endLevels.getLevelNodeIterEnd(i);
			while (levelNode != levelEnd){
				vert = *levelNode++;
				pair = pairs[vert];
				pair.endLevel = end;
				if (pair.startLevel == end){
					deletedVerts.push_back(vert);
					newLevels[end].push_back(vert);
				}
			}
		}

		//delete arcs
		for (auto vert : deletedVerts){
			int label = vert->getOldLable();
			graphVerts[label] = NULL;
			MeshGraphArcNode *arc = vert->getAdjacentNodes();
			while (arc){
				graphVerts[arc->getLable()]->deleteArc(label, graphArcPool);
				arc = arc->getNextNode();
			}
		}

		std::queue<MeshGraphVertexNode *> working;
		vector<set<MeshGraphVertexNode *, Comparer>> connectComponents;
		memset(visited, 0, sizeof(bool)*vertSize);

		//depth-first search for connected components
		int componentIndex = -1;
		for (auto vert : graphVerts){
			if (vert && !visited[vert->getOldLable()]){
				componentIndex++;
				visited[vert->getOldLable()] = true;
				connectComponents[componentIndex].insert(vert);
				MeshGraphArcNode *node = vert->getAdjacentNodes();
				while (node){
					working.push(graphVerts[node->getLable()]);
					node = node->getNextNode();
				}

				while (!working.empty()){
					MeshGraphVertexNode *adjant = working.front();
					if (!visited[adjant->getOldLable()]){
						visited[adjant->getOldLable()] = true;
						connectComponents[componentIndex].insert(vert);
						node = adjant->getAdjacentNodes();
						while (node){
							working.push(graphVerts[node->getLable()]);
							node = node->getNextNode();
						}
					}
					working.pop();
				}
			}
		}
		delete[] visited;

		std::sort(connectComponents.begin(), connectComponents.end(),
			[](const set<MeshGraphVertexNode *, Comparer> &left, const set<MeshGraphVertexNode *, Comparer> &right)
		{ return left.size() < right.size(); });

		int *h = new int[size], *l = new int[size];
		bool selectSecond = false, interchange = false;

		for (int i = 0; i < connectComponents.size(); i++){
			memset(h, 0, size*sizeof(int));
			memset(l, 0, size*sizeof(int));
			int componentSize = connectComponents[i].size();
			int levelWidth = newLevels[0].size();
			if (connectComponents[i].find(startLevels.getRoot()) != connectComponents[i].end())
				h[0] = levelWidth + 1;
			if (connectComponents[i].find(endLevels.getRoot()) != connectComponents[i].end())
				l[0] = levelWidth + 1;
			for (int m = 1; m < size; m++){
				h[m] = newLevels[m].size(); l[m] = newLevels[m].size();
				auto startIter = startLevels.getLevelNodeIterBeign(m);
				auto startEnd = startLevels.getLevelNodeIterEnd(m);
				auto endIter = endLevels.getLevelNodeIterBeign(m);
				auto endEnd = endLevels.getLevelNodeIterEnd(m);
				while (startIter != startEnd){
					if (connectComponents[i].find(*startIter) != connectComponents[i].end())
						h[m]++;
				}
				while (endIter != endEnd){
					if (connectComponents[i].find(*endIter) != connectComponents[i].end())
						l[m]++;
				}
			}
			int maxh = *std::max_element(h, h + size);
			int maxl = *std::max_element(l, l + size);

			if (i == 0)
				if (maxh > maxl || (maxh == maxl && !startSmaller)) selectSecond = true;

			if (maxh < maxl){
				for (auto vert : connectComponents[i])
					newLevels[pairs[vert].startLevel].push_back(vert);
			}
			else if (maxh > maxl){
				for (auto vert : connectComponents[i])
					newLevels[pairs[vert].endLevel].push_back(vert);
			}
			else{
				if (startSmaller){
					for (auto vert : connectComponents[i])
						newLevels[pairs[vert].startLevel].push_back(vert);
				}
				else{
					for (auto vert : connectComponents[i])
						newLevels[pairs[vert].endLevel].push_back(vert);
				}
			}
		}

		//restore original graph
		for(auto vert : deletedVerts){
			int label = vert->getOldLable();
			graphVerts[label] = vert;
			MeshGraphArcNode *node = vert->getAdjacentNodes();
			while (node){
				graphVerts[node->getLable()]->insertArc(label, graphArcPool);
				node = node->getNextNode();
			}
		}

		delete[] h; delete[] l;
		auto degreeCompare = [](const MeshGraphVertexNode *left, const MeshGraphVertexNode *right){return left->getDegree() < right->getDegree(); };
		for (int i = 0; i < size; i++)
			newLevels[i].sort(degreeCompare);

		//numbering
		MeshGraphVertexNode *root = startLevels.getRoot();
		int init = 0, pace = 1;
		std::function<bool(int)> compare{ [size](int i){ return i < size; } };

		if (startLevels.getRoot()->getDegree() > endLevels.getRoot()->getDegree()){
			interchange = true;
			root = endLevels.getRoot();
			init = size - 1; pace = -1;
			compare = [](int i){return i >= 0; };
		}

		auto pair = pairs[root];
		if (pair.startLevel == pair.endLevel)
			newLevels[init].remove(root);
		else{
			newLevels[pair.startLevel].remove(root);
			newLevels[pair.endLevel].remove(root);
		}

		int level = init;
		std::priority_queue<MeshGraphVertexNode *, vector<MeshGraphVertexNode *>, decltype(degreeCompare)> toBeAssigned(degreeCompare);

		//dealing with root
		processSigleVert(root, newLevels[level], working, toBeAssigned);

		//dealing with other verts in level 0
		while (!newLevels[level].empty()){
			MeshGraphVertexNode *vert = newLevels[level].front();
			newLevels[level].pop_front();
			processSigleVert(vert, newLevels[level], working, toBeAssigned);
		}

		for (level = init + pace; compare(level); level += pace){
			bool unfound = true;
			MeshGraphVertexNode *process = NULL;
			while (unfound && !working.empty()){
				MeshGraphVertexNode *vert = working.front();
				MeshGraphArcNode *node = vert->getAdjacentNodes();
				while (node){
					MeshGraphVertexNode *mayProcess = graphVerts[node->getLable()];
					if (std::find(newLevels[level].begin(), newLevels[level].end(), mayProcess)
						!= newLevels[level].end()){
						unfound = false;
						if (!process || process->getDegree() > mayProcess->getDegree())
							mayProcess = process;
					}
				}
				working.pop();
			}
			//clear queue
			while (!working.empty()) working.pop();

			Assert(!process);

			newLevels[level].remove(process);
			processSigleVert(process, newLevels[level], working, toBeAssigned);

			while (!newLevels[level].empty()){
				MeshGraphVertexNode *vert = newLevels[level].front();
				newLevels[level].pop_front();
				processSigleVert(vert, newLevels[level], working, toBeAssigned);
			}

		}

		int index = 0;
		if (interchange != selectSecond){
			for (auto vert : graphVerts)
				newNodeLables[index++] = vert->getNewLable();
		}
		else{
			for (auto vert : graphVerts)
				newNodeLables[index++] = (vertSize - 1) - vert->getNewLable();
		}

		//relabel elements
		std::priority_queue<OrderedElement> elementQueue;
		int *mem = allocAligned<int>(nodePerElementCount*numElements);
		int *unordered = new int[nodePerElementCount];

		int bound = numElements * nodePerElementCount;
		for (int i = 0; i < bound; i += nodePerElementCount){
			const int *iter = elements + i;
			for (int j = 0; j < nodePerElementCount; j++)
				unordered[j] = newNodeLables[iter[j]];
			elementQueue.push(OrderedElement(nodePerElementCount, unordered, mem));
		}

		for (int i = 0; i < bound; i += nodePerElementCount){
			int *iter = newElements + i;
			OrderedElement element = elementQueue.top();
			memcpy(iter, element.elementNodeIndices, sizeof(int)*nodePerElementCount);
			elementQueue.pop();
		}

		freeAligned(mem);
	}

}