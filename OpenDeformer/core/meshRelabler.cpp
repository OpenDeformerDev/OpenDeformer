#include "stdafx.h"
#include "meshRelabeler.h"
#include "datastructure.h"
#include "mesh.h"
#include <unordered_map>
#include <unordered_set>

namespace ODER{
	void MeshGraphVertexNode::insertArc(int lable, MemoryPool<MeshGraphArcNode>& pool){
		MeshGraphArcNode *node = pool.Alloc();
		node->setLable(lable);
		node->setNextNode(adjacent);
		adjacent = node;
		degree++;
	}

	void MeshGraphVertexNode::deleteArc(int lable, MemoryPool<MeshGraphArcNode>& pool){
		if (adjacent->getLable() != lable){
			MeshGraphArcNode *pre = adjacent;
			MeshGraphArcNode *now = pre->getNextNode();
			while (now && now->getLable() != lable){
				pre = now;
				now = now->getNextNode();
			}
			Assert(now != NULL);
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

	RootedLevelStructure::RootedLevelStructure(MeshGraphVertexNode *rootNode, const std::vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer){
		generateLevels(rootNode, graphVerts, visitedBuffer);
	}

	int RootedLevelStructure::getWidth() const{
		int levelSize = levels.size();
		int width = INT_MIN;
		for (int i = 0; i < levelSize; i++){
			int levelWidth = levels[i].size();
			if (levelWidth > width)
				width = levelWidth;
		}

		return width;
	}

	void RootedLevelStructure::generateLevels(MeshGraphVertexNode *rootNode, const std::vector<MeshGraphVertexNode *> &graphVerts, bool *visitedBuffer){
		root = rootNode;
		visitedBuffer[root->getOldLable()] = true;
		MeshGraphArcNode *node = root->getAdjacentNodes();
		std::queue<MeshGraphVertexNode *> lableStarters;
		//level 1 construction
		levels.emplace_back();
		while (node){
			int lable = node->getLable();
			MeshGraphVertexNode *vertNode = graphVerts[lable];
			levels[0].push_back(vertNode);
			lableStarters.push(vertNode);
			visitedBuffer[lable] = true;
			node = node->getNextNode();
		}
		//other levels
		int levelNum = 0;
		while (!lableStarters.empty()){
			levelNum++;
			levels.emplace_back();
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
		}
		levels.pop_back();
	}

	MeshRelabeler::OrderedElement::OrderedElement(int nodePerElementCount, const int *indices, int *mem){
		elementNodeIndices = mem;
		memcpy(elementNodeIndices, indices, sizeof(int)*nodePerElementCount);
		score = 0;
		for (int i = 0; i < nodePerElementCount; i++)
			score += elementNodeIndices[i];
	}

	MeshRelabeler::MeshRelabeler(int nodeCount) :nextFreeLable(0), graphVerts(nodeCount){
		MeshGraphVertexNode *mem = allocAligned<MeshGraphVertexNode>(nodeCount);
		for (int i = 0; i < nodeCount; i++){
			MeshGraphVertexNode *pointer = mem + i;
			new (pointer)MeshGraphVertexNode();
			graphVerts[i] = pointer;
			graphVerts[i]->setOldLable(i);
		}
	}

	MeshRelabeler::MeshRelabeler(const Mesh& mesh) :nextFreeLable(0), graphVerts(mesh.getNodeCount()){
		int nodeCount = mesh.getNodeCount();
		MeshGraphVertexNode *mem = allocAligned<MeshGraphVertexNode>(nodeCount);
		for (int i = 0; i < nodeCount; i++){
			MeshGraphVertexNode *pointer = mem + i;
			new (pointer)MeshGraphVertexNode();
			graphVerts[i] = pointer;
			graphVerts[i]->setOldLable(i);
		}
	}

	void MeshRelabeler::generateGraphFromElementIndices(int elementCount, int nodePerElementCount, int *elements){
		std::unordered_set<MeshGraphEdge, std::function<size_t(const MeshGraphEdge&)>>
			edges(2 * elementCount,
			[](const MeshGraphEdge& edge)->size_t{ return ((edge.tail * (edge.tail + 1)) >> 1) + edge.head; });
		const int* element = elements;
		for (int i = 0; i < elementCount; i++){
			for (int j = 0; j < nodePerElementCount; j++){
				for (int k = j + 1; k < nodePerElementCount; k++)
					edges.insert(MeshGraphEdge(element[j], element[k]));
			}
			element += nodePerElementCount;
		}

		for (auto edge : edges){
			graphVerts[edge.head]->insertArc(edge.tail, graphArcPool);
			graphVerts[edge.tail]->insertArc(edge.head, graphArcPool);
		}
	}

	void MeshRelabeler::generateGraphFromMesh(const Mesh& mesh){
		std::unordered_set<MeshGraphEdge, std::function<size_t(const MeshGraphEdge&)>> 
			edges(2 * mesh.getElementCount(), 
			[](const MeshGraphEdge& edge)->size_t{ return ((edge.tail * (edge.tail + 1)) >> 1) + edge.head; });
		int elementCount = mesh.getElementCount();
		int nodePerElementCount = mesh.getNodePerElementCount();
		for (int i = 0; i < elementCount; i++){
			const int* element = mesh.getElementNodeReference(i);
			for (int j = 0; j < nodePerElementCount; j++){
				for (int k = j + 1; k < nodePerElementCount; k++){
					edges.insert(MeshGraphEdge(element[j], element[k]));
				}
			}
		}

		for (auto edge : edges){
			graphVerts[edge.head]->insertArc(edge.tail, graphArcPool);
			graphVerts[edge.tail]->insertArc(edge.head, graphArcPool);
		}
	}

	void MeshRelabeler::processSigleVert(MeshGraphVertexNode * vert, std::list<MeshGraphVertexNode *>& levelNodes,
		std::queue<MeshGraphVertexNode *>& working, 
		std::priority_queue<MeshGraphVertexNode *, std::vector<MeshGraphVertexNode *>, DegreeComparer> &toBeAssigned){

		vert->setNewLable(nextFreeLable++);
		working.push(vert);
		MeshGraphArcNode *node = vert->getAdjacentNodes();
		while (node){
			MeshGraphVertexNode *vertex = graphVerts[node->getLable()];
			auto found = std::find(levelNodes.begin(), levelNodes.end(), vertex);
			if (found != levelNodes.end()){
				toBeAssigned.push(vertex);
				levelNodes.erase(found);
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

	void MeshRelabeler::setNewNodeLables(int *newNodeLables){
		//find endpoints of pseudo-diameter
		int vertSize = graphVerts.size();
		bool *visited = new bool[vertSize];
		memset(visited, 0, sizeof(bool)*vertSize);

		DegreeComparer degreeCompare;
		MeshGraphVertexNode *start = *std::min_element(graphVerts.begin(), graphVerts.end(), degreeCompare);

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
				mayEndLevels.Clear();
				mayEndLevels.generateLevels(*lastLevelNode++, graphVerts, visited);
				int levelSize = mayEndLevels.getLevelSize();
				if (levelSize > startLevels.getLevelSize()){
					find = false;
					std::swap(startLevels, mayEndLevels);
					break;
				}
				int width = mayEndLevels.getWidth();
				if (width < endWidth){
					std::swap(endLevels, mayEndLevels);
					endWidth = width;
				}
			}
		}

		//minimizing level width
		int startLevelSize = startLevels.getLevelSize(), endLevelSize = endLevels.getLevelSize();
		int size = std::max(startLevelSize, endLevelSize);
		int startLevelWidth = startLevels.getWidth();
		bool startSmaller = endLevels.getWidth() > startLevelWidth;

		struct OrderedPair{
			int startLevel;
			int endLevel;
		};
		
		OrderedPair *pairs = new OrderedPair[vertSize];
		std::deque<MeshGraphVertexNode *> deletedVerts;
		std::queue<MeshGraphVertexNode *> remainedVertices;
		std::vector<std::list<MeshGraphVertexNode *>> newLevels(size, std::list<MeshGraphVertexNode *>());

		pairs[startLevels.getRoot()->getOldLable()] = OrderedPair{ 0, INT_MIN };
		for (int i = 1; i < startLevelSize; i++){
			auto levelNode = startLevels.getLevelNodeIterBeign(i);
			auto levelEnd = startLevels.getLevelNodeIterEnd(i);
			while (levelNode != levelEnd){
				MeshGraphVertexNode *v = *levelNode++;
				pairs[v->getOldLable()] = OrderedPair{ i, INT_MIN };
			}
		}

		MeshGraphVertexNode *vert = endLevels.getRoot();
		OrderedPair& endLevelRootPair = pairs[vert->getOldLable()];
		int end = endLevelSize - 1;
		endLevelRootPair.endLevel = end;
		if (endLevelRootPair.startLevel == end){
			deletedVerts.push_back(vert);
			newLevels[end].push_back(vert);
		}
		else
			remainedVertices.push(vert);
		for (int i = 1; i < endLevelSize; i++){
			end--;
			auto levelNode = endLevels.getLevelNodeIterBeign(i);
			auto levelEnd = endLevels.getLevelNodeIterEnd(i);
			while (levelNode != levelEnd){
				vert = *levelNode++;
				OrderedPair& pair = pairs[vert->getOldLable()];
				pair.endLevel = end;
				if (pair.startLevel == end){
					deletedVerts.push_back(vert);
					newLevels[end].push_back(vert);
				}
				else
					remainedVertices.push(vert);
			}
		}

		//delete arcs
		for (auto vert : deletedVerts){
			int label = vert->getOldLable();
			MeshGraphArcNode *arc = vert->getAdjacentNodes();
			while (arc){
				graphVerts[arc->getLable()]->deleteArc(label, graphArcPool);
				arc = arc->getNextNode();
			}
		}

		std::queue<MeshGraphVertexNode *> working;
		std::vector<std::vector<MeshGraphVertexNode *>> connectComponents;
		memset(visited, 0, sizeof(bool)*vertSize);

		//depth-first search for connected components
		int componentIndex = -1;
		while (!remainedVertices.empty()){
			MeshGraphVertexNode *vert = remainedVertices.front();
			if (!visited[vert->getOldLable()]){
				componentIndex++;
				connectComponents.emplace_back();
				visited[vert->getOldLable()] = true;
				connectComponents[componentIndex].push_back(vert);
				MeshGraphArcNode *node = vert->getAdjacentNodes();
				while (node){
					working.push(graphVerts[node->getLable()]);
					node = node->getNextNode();
				}

				while (!working.empty()){
					MeshGraphVertexNode *adjant = working.front();
					if (!visited[adjant->getOldLable()]){
						visited[adjant->getOldLable()] = true;
						connectComponents[componentIndex].push_back(adjant);
						node = adjant->getAdjacentNodes();
						while (node){
							working.push(graphVerts[node->getLable()]);
							node = node->getNextNode();
						}
					}
					working.pop();
				}
			}
			remainedVertices.pop();
		}
		delete[] visited;

		std::sort(connectComponents.begin(), connectComponents.end(),
			[](const std::vector<MeshGraphVertexNode *> &left, const std::vector<MeshGraphVertexNode *> &right)
		{ return left.size() < right.size(); });

		int *mem =new int[size * 2];
		int *firstLevelCount = mem, *secondLevelCount = mem + size;
		bool selectSecond = false, interchange = false;

		std::queue<OrderedPair> toBeInsert;
		for (int i = 0; i < connectComponents.size(); i++){
			memset(firstLevelCount, 0, size*sizeof(int));
			memset(secondLevelCount, 0, size*sizeof(int));
			for (auto vertex : connectComponents[i]){
				OrderedPair pair = pairs[vertex->getOldLable()];
				firstLevelCount[pair.startLevel]++;
				secondLevelCount[pair.endLevel]++;
				toBeInsert.push(pair);
			}
			int maxFirstCount = *std::max_element(firstLevelCount, firstLevelCount + size);
			int maxSecondCount = *std::max_element(secondLevelCount, secondLevelCount + size);

			if (i == 0)
				if (maxFirstCount > maxSecondCount || (maxFirstCount == maxSecondCount && !startSmaller)) selectSecond = true;

			if (maxFirstCount < maxSecondCount){
				for (auto vertex : connectComponents[i]){
					newLevels[toBeInsert.front().startLevel].push_back(vertex);
					toBeInsert.pop();
				}
			}
			else if (maxFirstCount > maxSecondCount){
				for (auto vertex : connectComponents[i]){
					newLevels[toBeInsert.front().endLevel].push_back(vertex);
					toBeInsert.pop();
				}
			}
			else{
				if (startSmaller){
					for (auto vertex : connectComponents[i]){
						newLevels[toBeInsert.front().startLevel].push_back(vertex);
						toBeInsert.pop();
					}
				}
				else{
					for (auto vertex : connectComponents[i]){
						newLevels[toBeInsert.front().endLevel].push_back(vertex);
						toBeInsert.pop();
					}
				}
			}
		}

		//restore original graph
		for (auto vert : deletedVerts){
			int label = vert->getOldLable();
			MeshGraphArcNode *node = vert->getAdjacentNodes();
			while (node){
				graphVerts[node->getLable()]->insertArc(label, graphArcPool);
				node = node->getNextNode();
			}
		}

		delete[] mem;
		//DegreeComparer degreeCompare;
		for (int i = 0; i < size; i++)
			newLevels[i].sort(degreeCompare);

		//numbering
		MeshGraphVertexNode *root = startLevels.getRoot();
		int init = 0, pace = 1;
		std::function<bool(int)> compare{ [size](int i) { return i < size; } };

		if (startLevels.getRoot()->getDegree() > endLevels.getRoot()->getDegree()){
			interchange = true;
			root = endLevels.getRoot();
			init = size - 1; pace = -1;
			compare = [](int i) { return i >= 0; };
		}

		OrderedPair rootPair = pairs[root->getOldLable()];
		if (rootPair.startLevel == rootPair.endLevel)
			newLevels[init].remove(root);
		else{
			newLevels[rootPair.startLevel].remove(root);
			newLevels[rootPair.endLevel].remove(root);
		}

		delete[] pairs;

		int level = init;
		std::vector<MeshGraphVertexNode *> container;
		container.reserve(startLevelWidth);
		std::priority_queue<MeshGraphVertexNode *, std::vector<MeshGraphVertexNode *>, DegreeComparer> 
			toBeAssigned(degreeCompare, std::move(container));

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
							process = mayProcess;
					}
					node = node->getNextNode();
				}
				working.pop();
			}
			//clear queue
			std::queue<MeshGraphVertexNode *>().swap(working);

			Assert(process != NULL);

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
	}

	void MeshRelabeler::getNewLables(int elementCount, int nodePerElementCount, int *newNodeLables, int *elements){
		generateGraphFromElementIndices(elementCount, nodePerElementCount, elements);
		setNewNodeLables(newNodeLables);
		//relabel elements
		OrderedElement *elementQueue = allocAligned<OrderedElement>(elementCount);
		int *mem = allocAligned<int>(elementCount*nodePerElementCount);
		int *elementIndices = new int[nodePerElementCount];

		const int *eleIter = elements;
		for (int i = 0; i < elementCount; i++){
			for (int j = 0; j < nodePerElementCount; j++)
				elementIndices[j] = newNodeLables[eleIter[j]];
			elementQueue[i] = OrderedElement(nodePerElementCount, elementIndices, mem + i*nodePerElementCount);
			eleIter += nodePerElementCount;
		}
		
		std::sort(elementQueue, elementQueue + elementCount);
		int *newElementsIter = elements;
		for (int i = 0; i < elementCount; i++){
			memcpy(newElementsIter, elementQueue[i].elementNodeIndices, sizeof(int)*nodePerElementCount);
			newElementsIter += nodePerElementCount;
		}

		delete[] elementIndices;
		freeAligned(mem);
		freeAligned(elementQueue);
	}

	void MeshRelabeler::getNewLables(int *newNodeLables, Mesh& mesh){
		generateGraphFromMesh(mesh);
		setNewNodeLables(newNodeLables);
		//relabel elements
		int elementCount = mesh.getElementCount();
		int nodePerElementCount = mesh.getNodePerElementCount();
		OrderedElement *elementQueue = allocAligned<OrderedElement>(elementCount);
		int *mem = allocAligned<int>(elementCount*nodePerElementCount);
		int *elementIndices = new int[nodePerElementCount];

		for (int i = 0; i < elementCount; i++){
			const int *element = mesh.getElementNodeReference(i);
			for (int j = 0; j < nodePerElementCount; j++)
				elementIndices[j] = newNodeLables[element[j]];
			elementQueue[i] = OrderedElement(nodePerElementCount, elementIndices, mem + i*nodePerElementCount);
		}

		std::sort(elementQueue, elementQueue + elementCount);
		for (int i = 0; i < elementCount; i++)
			mesh.setElement(i, elementQueue[i].elementNodeIndices);

		delete[] elementIndices;
		freeAligned(mem);
		freeAligned(elementQueue);
	}
}