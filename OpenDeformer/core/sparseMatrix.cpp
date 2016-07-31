#include "stdafx.h"
#include "sparseMatrix.h"

namespace ODER{
	void SparseMatrixAssembler::addEntry(int row, int column, double data){
		auto found = rowEntries[column].find(row);
		if (found != rowEntries[column].end())
			found->second += data;
		else{
			std::pair<int, double> entry(row, data);
			rowEntries[column].insert(entry);
		}
	}

	void SparseMatrixAssembler::removeRow(int row){
		for (int i = 0; i < numColumns; i++){
			auto iter = rowEntries[i].lower_bound(row);
			if (iter == rowEntries[i].end()) continue;
			if (iter->first == row)
				rowEntries[i].erase(iter++);
			do{
				int newRow = iter->first - 1;
				double data = iter->second;
				rowEntries[i].erase(iter++);
				rowEntries[i].insert(iter, std::pair<int, double>(newRow, data));
			} while (iter != rowEntries[i].end());
		}
	}

	void SparseMatrixAssembler::removeColumn(int column){
		numColumns--;
		rowEntries.erase(rowEntries.begin() + column);
	}

	void SparseMatrixAssembler::removeColumnAndRow(int num){
		numColumns--;
		rowEntries.erase(rowEntries.begin() + num);

		for (int i = 0; i < numColumns; i++){
			auto iter = rowEntries[i].lower_bound(num);
			if (iter == rowEntries[i].end()) continue;
			if (iter->first == num)
				rowEntries[i].erase(iter++);
			do{
				int newRow = iter->first - 1;
				double data = iter->second;
				rowEntries[i].erase(iter++);
				rowEntries[i].insert(iter, std::pair<int, double>(newRow, data));
			} while (iter != rowEntries[i].end());
		}
	}

	SparseMatrix::SparseMatrix(const Assembler &assembler){
		numColumns = assembler.numColumns;
		pcol = (int*)allocAligned<int>(numColumns + 1);

		pcol[0] = 0;
		for (int i = 1; i < numColumns + 1; i++){
			pcol[i] = pcol[i - 1] + assembler.rowEntries[i - 1].size();
		}

		value = (double *)allocAligned<double>(pcol[numColumns]);
		rows = (int *)allocAligned<int>(pcol[numColumns]);

		int j = 0;
		for (int i = 0; i < numColumns; i++){
			auto pos = assembler.rowEntries[i].begin();
			while (pos != assembler.rowEntries[i].end()){
				rows[j] = pos->first;
				value[j] = pos->second;
				j++;
				pos++;
			}
		}
	}

	std::vector<std::unordered_map<int, int>> SparseMatrix::getIndices() const {
		std::vector<std::unordered_map<int, int>> indices(numColumns);

		for (int i = 0; i < numColumns; i++) {
			int index = pcol[i];
			int end = pcol[i + 1];

			while (index < end)
				indices[i].emplace(rows[index], index++);
		}

		return indices;
	}

	SparseMatrix::~SparseMatrix(){
		freeAligned(value);
		freeAligned(rows);
		freeAligned(pcol);
	}

}