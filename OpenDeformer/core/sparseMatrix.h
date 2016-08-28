#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_SPARSEMARTIX_H
#define ODER_CORE_SPARSEMARTIX_H

#include "oder.h"
#include "memory.h"
#include "latool.h"
#include <map>
#include <vector>
#include <unordered_map>
#include <array>

namespace ODER{
	class SparseMatrixAssembler{
	public:
		SparseMatrixAssembler(int nc) :numColumns(nc){
			std::map<int, double> empty;
			for (int i = 0; i < numColumns; i++)
				rowEntries.push_back(empty);
		}
		void addEntry(int row, int column, double data);
		void removeRow(int row);
		void removeColumn(int column);
		void removeColumnAndRow(int num);

		~SparseMatrixAssembler() = default;
		int numColumns;
		std::vector<std::map<int, double>> rowEntries;
	};

	class SparseMatrix{
	public:
		using Assembler = SparseMatrixAssembler;
		SparseMatrix(const Assembler& assembler);
		SparseMatrix(const SparseMatrix&) = delete;
		SparseMatrix& operator=(const SparseMatrix&) = delete;
		SparseMatrix(SparseMatrix&& mat) noexcept{
			numColumns = mat.numColumns;
			value = mat.value;
			rows = mat.rows;
			pcol = mat.pcol;

			mat.numColumns = 0;
			mat.value = NULL; mat.rows = NULL; mat.pcol = NULL;
		}
		SparseMatrix& operator=(SparseMatrix&& mat) noexcept{
			std::swap(numColumns, mat.numColumns);
			std::swap(value, mat.value);
			std::swap(rows, mat.rows);
			std::swap(pcol, mat.pcol);
		}
		void setZeros() {
			Initiation(value, pcol[numColumns]);
		}
		void Add(double alpha, const SparseMatrix& mat) {
			Assert(numColumns == mat.numColumns);
			Assert(pcol[numColumns] == mat.pcol[numColumns]);

			int dataCount = pcol[numColumns];
			for (int i = 0; i < dataCount; i++)
				value[i] += alpha * mat.value[i];
		}
		void Scale(double alpha) {
			int dataCount = pcol[numColumns];
			for (int i = 0; i < dataCount; i++)
				value[i] *= alpha;
		}
		int getNumColumns() const { return numColumns; }
		std::vector<std::unordered_map<int, int>> getIndices() const;
		~SparseMatrix();
	private:
		int numColumns;
		//CSC format
		double* value;
		int* rows;//row number of each entry
		int* pcol;//the index of the start of each column and the end

		friend class EigenSolver;
	};

	constexpr size_t diagIndicesGen(size_t length, size_t index) noexcept { return index * length - (((index - 1)*index) >> 1); }
	template<size_t... index> constexpr std::array<size_t, sizeof...(index)> getDiagIndices(std::index_sequence<index...>&& seqs) noexcept {
		return { diagIndicesGen(sizeof...(index), index)... };
	}

	template<int blockLength, int blockWidth> class BlockedSymSparseMatrixAssembler{
	public:
		BlockedSymSparseMatrixAssembler(int columnCount) 
			: numColumn(columnCount), numBlockColumn(numColumn / blockLength), numRemainedColumn(numColumn - numBlockColumn*blockLength),
			blockEntries(numBlockColumn), remainedEntryCount(NULL){

			if (numRemainedColumn > 0){
				remainedEntryCount = new int[numRemainedColumn];
				Initiation(remainedEntryCount, numRemainedColumn);
			}
		}

		void addEntry(int row, int column, double data){
			Assert(row >= column);
			int blockCloumn = column / blockLength;
			int columnStart = blockCloumn*blockLength;
			int subColumn = column - columnStart;
			if (blockCloumn < numBlockColumn){
				if (row - columnStart < blockLength){
					int subRow = row - columnStart;
					int index = subRow - subColumn + diagIndices[subColumn];
					auto found = blockEntries[blockCloumn].find(columnStart);
					if (found != blockEntries[blockCloumn].end())
						found->second[index] += data;
					else{
						constexpr int numBlockEntries = (blockLength + 1) * blockLength / 2;
						double *vals = entryMem.Alloc(numBlockEntries);
						vals[index] = data;
						blockEntries[blockCloumn].insert(std::pair<int, double*>(columnStart, vals));
					}
				}

				else{
					int normalBlockStartRow = columnStart + blockLength;
					int blockIndex = (row - normalBlockStartRow) / blockWidth;
					int blockStartRow = normalBlockStartRow + blockIndex*blockWidth;
					int subRow = row - blockStartRow;
					
					int width = numColumn - blockStartRow;
					if (width >= blockWidth){
						auto found = blockEntries[blockCloumn].find(blockStartRow);
						if (found != blockEntries[blockCloumn].end())
							found->second[subColumn*blockWidth + subRow] += data;
						else{
							constexpr int numBlockEntries = blockLength*blockWidth;
							double *vals = entryMem.Alloc(numBlockEntries);
							vals[subColumn*blockWidth + subRow] = data;
							blockEntries[blockCloumn].insert(std::pair<int, double*>(blockStartRow, vals));
						}
					}
					else{
						auto found = blockEntries[blockCloumn].find(blockStartRow);
						if (found != blockEntries[blockCloumn].end())
							found->second[subColumn*width + subRow] += data;
						else{
							int numBlockEntries = blockLength*width;
							double *vals = entryMem.Alloc(numBlockEntries);
							vals[subColumn*width + subRow] = data;
							blockEntries[blockCloumn].insert(std::pair<int, double*>(blockStartRow, vals));
						}
					}
				}
			}
			else{
				int subRow = row - columnStart;
				int index = subColumn*numRemainedColumn + subRow;
				auto found = remainedEntries.find(index);
				if (found != remainedEntries.end())
					found->second += data;
				else{
					remainedEntryCount[subColumn]++;
					remainedEntries.insert(std::pair<int, double>(index, data));
				}
			}
		}

		~BlockedSymSparseMatrixAssembler(){ if (remainedEntryCount) delete[] remainedEntryCount; }

	private:

		int numColumn;
		int numBlockColumn;
		int numRemainedColumn;
		std::vector<std::map<int, double*>> blockEntries;
		std::map<int, double> remainedEntries;
		int *remainedEntryCount;

		static constexpr std::array<size_t, blockLength> diagIndices = getDiagIndices(std::make_index_sequence<blockLength>());
		MemoryArena<double> entryMem;

		friend class BlockedSymSparseMatrix <blockLength, blockWidth>;
	};

	template<int blockLength, int blockWidth> class BlockedSymSparseMatrix{
	public:
		using Assembler = BlockedSymSparseMatrixAssembler<blockLength, blockWidth>;
		BlockedSymSparseMatrix(const Assembler& assembler){
			numColumns = assembler.numColumn;
			numBlockColumn = assembler.numBlockColumn;
			numRemainedColumn = assembler.numRemainedColumn;

			blockPcol = allocAligned<int>(numBlockColumn + numRemainedColumn + 1);
			blockColumnOris = allocAligned<int>(numBlockColumn + numRemainedColumn + 1);

			constexpr int regularSize = blockLength * blockWidth;
			constexpr int diagSize = (blockLength + 1) * blockLength / 2;

			blockPcol[0] = 0;
			blockColumnOris[0] = 0;

			int dataCount = 0;
			for (int i = 0; i < numBlockColumn; i++){
				//counting blocks
				int blockCount = assembler.blockEntries[i].size();
				blockPcol[i + 1] = blockPcol[i] + blockCount;
				//counting datas
				if (blockCount != 0){
					int blockDataCount = blockCount * regularSize;
					auto start = assembler.blockEntries[i].cbegin();
					bool hasDiag = start->first == i * blockLength;
					if (hasDiag)
						blockDataCount -= (regularSize - diagSize);
					auto last = --(assembler.blockEntries[i].cend());
					if ((!hasDiag || start != last) && numColumns - last->first < blockWidth)
						blockDataCount -= (last->first + blockWidth - numColumns)*blockLength;
					dataCount += blockDataCount;
				}
				blockColumnOris[i + 1] = dataCount;
			}

			for (int i = 0; i < numRemainedColumn; i++){
				int entryCount = assembler.remainedEntryCount[i];
				blockPcol[i + numBlockColumn + 1] = blockPcol[i + numBlockColumn] + entryCount;
				blockColumnOris[i + numBlockColumn + 1] = blockColumnOris[i + numBlockColumn] + entryCount;
			}

			int blockRowCount = blockPcol[numBlockColumn + numRemainedColumn];
			dataCount += blockRowCount - blockPcol[numBlockColumn];

			values = allocAligned<double>(dataCount);
			blockRows = allocAligned<int>(blockRowCount);

			int index = 0;
			double *valueIter = values;
			for (int i = 0; i < numBlockColumn; i++){
				if (blockPcol[i + 1] != blockPcol[i]){
					auto pos = assembler.blockEntries[i].cbegin();
					auto last = assembler.blockEntries[i].cend();
					auto end = last--;
					if (pos->first == i*blockLength){
						blockRows[index++] = pos->first;
						memcpy(valueIter, pos->second, sizeof(double)*diagSize);
						valueIter += diagSize;
						pos++;
					}
					if (pos != end){
						while (pos != last){
							blockRows[index++] = pos->first;
							memcpy(valueIter, pos->second, sizeof(double)*regularSize);
							valueIter += regularSize;
							pos++;
						}
						if (numColumns - pos->first >= blockWidth){
							blockRows[index++] = pos->first;
							memcpy(valueIter, pos->second, sizeof(double)*regularSize);
							valueIter += regularSize;
						}
						else{
							blockRows[index++] = pos->first;
							memcpy(valueIter, pos->second, sizeof(double)*(numColumns - pos->first)*blockLength);
							valueIter += (numColumns - pos->first)*blockLength;
						}
					}
				}
			}

			int rowOffset = blockLength*numBlockColumn;
			auto pos = assembler.remainedEntries.cbegin();
			auto end = assembler.remainedEntries.cend();
			while (pos != end){
				blockRows[index++] = pos->first % numRemainedColumn + rowOffset;
				*valueIter++ = pos->second;
				pos++;
			}
		}


		void setZeros() {
			Initiation(values, blockColumnOris[numBlockColumn + numRemainedColumn]);
		}

		void Add(double alpha, const BlockedSymSparseMatrix& mat) {
			Assert(numColumns == mat.numColumns);
			Assert(numBlockColumn == mat.numBlockColumn);
			Assert(blockColumnOris[numBlockColumn + numRemainedColumn] == mat.blockColumnOris[numBlockColumn + numRemainedColumn]);

			int dataCount = blockColumnOris[numBlockColumn + numRemainedColumn];
			for (int i = 0; i < dataCount; i++)
				values[i] += alpha * mat.values[i];
		}

		void Scale(double alpha) {
			int dataCount = blockColumnOris[numBlockColumn + numRemainedColumn];
			for (int i = 0; i < dataCount; i++)
				values[i] *= alpha;
		}

		void getColumn(int column, SparseVector& vector) const{
			vector.Clear();
			constexpr int regularSize = blockLength * blockWidth;
			constexpr int diagSize = (blockLength + 1) * blockLength / 2;

			int blockIndex = column / blockLength;
			int blockStartColumn = blockIndex * blockLength;
			if (blockIndex < numBlockColumn){
				int offset = column - blockStartColumn;
				int index = blockPcol[blockIndex], end = blockPcol[blockIndex + 1];

				if (index == end) return;

				const double *vals = values + blockColumnOris[blockIndex];
				if (blockRows[index] == blockStartColumn){
					int start = diagIndices[offset];
					for (int i = 0; i < blockLength - offset; i++) {
						if (vals[start + i] != 0.0)
						    vector.emplaceBack(column + i, vals[start + i]);
					}
					index++;
					vals += diagSize;
				}

				if (index == end) return;

				while (index < end - 1){
					int row = blockRows[index++];
					int start = offset*blockWidth;
					for (int i = 0; i < blockWidth; i++) {
						if (vals[start + i] != 0.0)
						    vector.emplaceBack(row + i, vals[start + i]);
					}
					vals += regularSize;
				}

				int row = blockRows[index];
				int mayDegenWidth = numColumns - row;
				if (mayDegenWidth >= blockWidth){
					int start = offset * blockWidth;
					for (int i = 0; i < blockWidth; i++) {
						if (vals[start + i] != 0.0)
						    vector.emplaceBack(row + i, vals[start + i]);
					}
				}
				else{
					int start = offset * mayDegenWidth;
					for (int i = 0; i < mayDegenWidth; i++) {
						if (vals[start + i] != 0.0)
						    vector.emplaceBack(row + i, vals[start + i]);
					}
				}
			}
			else{
				int start = column - blockStartColumn + numBlockColumn;
				int rowIndexStart = blockPcol[start];
				int entryCount = blockPcol[start + 1] - rowIndexStart;

				const int* row = blockRows + rowIndexStart;
				const double* val = values + blockColumnOris[start];
				for (int i = 0; i < entryCount; i++) {
					if (val[i] != 0.0)
					    vector.emplaceBack(row[i], val[i]);
				}
			}
		}

		void addEntry(int row, int column, double data) {
			Assert(row >= column);
			constexpr int regularSize = blockLength * blockWidth;
			constexpr int diagSize = (blockLength + 1) * blockLength / 2;

			int blockIndex = column / blockLength;
			int startColumn = blockIndex * blockLength;
			int subColumn = column - startColumn;
			if (blockIndex < numBlockColumn) {
				bool hasDiag = blockIndex * blockLength == blockRows[blockPcol[blockIndex]];
				int blockStartRow = startColumn;
				int *pRowIndex = &blockRows[blockPcol[blockIndex]];
				if (row >= startColumn + blockLength) {
					blockStartRow += blockLength + ((row - startColumn - blockLength) / blockWidth) * blockWidth;
					pRowIndex = std::lower_bound(&blockRows[blockPcol[blockIndex]], &blockRows[blockPcol[blockIndex + 1]], blockStartRow);
				}
				Assert(*pRowIndex == blockStartRow);
				int subRow = row - blockStartRow;
				int preBlockCount = pRowIndex - &blockRows[blockPcol[blockIndex]];
				int preEntryCount = preBlockCount * regularSize;
				if (hasDiag) preEntryCount -= (regularSize - diagSize);
				if (pRowIndex != &blockRows[blockPcol[blockIndex]] || !hasDiag) {
					double *blockValueStart = values + blockColumnOris[blockIndex] + preEntryCount;
					double *val = blockValueStart + subColumn * std::min(blockWidth, numColumns - blockStartRow) + subRow;
					*val += data;
				}
				else {
					int index = subRow - subColumn + diagIndices[subColumn];
					double *val = values + blockColumnOris[blockIndex] + index;
					*val += data;
				}
			}
			else {
				int columnIndex = blockIndex + subColumn;
				int *pRowIndex = std::lower_bound(&blockRows[blockPcol[columnIndex]], &blockRows[blockPcol[columnIndex + 1]], row);
				Assert(*pRowIndex == row);
				double *val = values + blockColumnOris[columnIndex] + (pRowIndex - &blockRows[blockPcol[columnIndex]]);
				*val += data;
			}
		}

		void addEntry(int row, int column, double data, 
			const std::vector<std::unordered_map<int, int>>& indices) {
			auto pair = indices[column].find(row);
			Assert(pair != indices[column].cend());
			values[pair->second] += data;
		}

		void getDiagonal(double* diags) const{
			int i = 0;
			const double *vals = values;
			for (; i < numBlockColumn; i++){
				int columnStart = blockPcol[i];
				for (int j = 0; j < blockLength; j++)
					*diags++ = vals[diagIndices[j]];
				vals += blockColumnOris[i + 1] - blockColumnOris[i];
			}
			for (; i < numBlockColumn + numRemainedColumn; i++){
				*diags++ = *vals;
				vals += blockColumnOris[i + 1] - blockColumnOris[i];
			}
		}

		void getDiaginalWithCheck(double *diags) const{
			int column = 0, i = 0;
			const int* pcol = blockPcol;
			const double *vals = values;
			for (; column < numBlockColumn * blockLength; column += blockLength, i++){
				int columnStart = *pcol++;
				if (blockRows[columnStart] == column){
					for (int j = 0; j < blockLength; j++)
						*diags++ = vals[diagIndices[j]];
				}
				else{
					memset(diags, 0, sizeof(double)*blockLength);
					diags += blockLength;
				}
				vals += blockColumnOris[i + 1] - blockColumnOris[i];
			}
			for (; column < numColumns; column++, i++){
				int columnStart = *pcol++;
				if (blockRows[columnStart] == column)
					*diags++ = *vals;
				else
					*diags++ = 0.0;
				vals += blockColumnOris[i + 1] - blockColumnOris[i];
			}
		}

		std::vector<std::vector<std::pair<int, int>>> getFullIndices() const {
			std::vector<std::vector<std::pair<int, int>>> entryPairs(numColumns);
			const int* colPointer = blockPcol;

			const int blockBoundary = numBlockColumn * blockLength;
			int index = 0;
			for (int column = 0; column < blockBoundary; column += blockLength) {
				int blockIndex = *colPointer++;
				int blockEnd = *colPointer;

				if (blockIndex == blockEnd) continue;

				//diag block
				if (blockRows[blockIndex] == column) {
					for (int i = 0; i < blockLength; i++) {
						entryPairs[column + i].push_back(std::make_pair(column + i, index++));
						for (int j = i + 1; j < blockLength; j++) {
							entryPairs[column + i].push_back(std::make_pair(column + j, index));
							entryPairs[column + j].push_back(std::make_pair(column + i, index));
							index += 1;
						}
					}
					blockIndex += 1;
				}

				if (blockIndex == blockEnd) continue;

				while (blockIndex < blockEnd - 1) {
					int row = blockRows[blockIndex++];
					for (int i = 0; i < blockLength; i++) {
						for (int j = 0; j < blockWidth; j++) {
							entryPairs[column + i].push_back(std::make_pair(row + j, index));
							entryPairs[row + j].push_back(std::make_pair(column + i, index));
							index += 1;
						}
					}
				}

				//last block
				int row = blockRows[blockIndex];
				int lastWidth = numColumns - row;
				if (lastWidth >= blockWidth) lastWidth = blockWidth;
				for (int i = 0; i < blockLength; i++) {
					for (int j = 0; j < lastWidth; j++) {
						entryPairs[column + i].push_back(std::make_pair(row + j, index));
						entryPairs[row + j].push_back(std::make_pair(column + i, index));
						index += 1;
					}
				}
			}

			for (int column = blockBoundary; column < numColumns; column++) {
				int iter = *colPointer++;
				int end = *colPointer;

				if (iter == end) continue;

				int row = blockRows[iter];
				if (row == column) {
					entryPairs[column].push_back(std::make_pair(column, index++));
					iter += 1;
				}

				while (iter < end) {
					row = blockRows[iter++];
					entryPairs[column].push_back(std::make_pair(row, index));
					entryPairs[row].push_back(std::make_pair(column, index));
					index += 1;
				}
			}

			return entryPairs;
		}

		std::vector<std::unordered_map<int, int>> getIndices() const {
			std::vector<std::unordered_map<int, int>> indices(numColumns);

			const int* colPointer = blockPcol;

			const int blockBoundary = numBlockColumn * blockLength;
			int index = 0;
			for (int column = 0; column < blockBoundary; column += blockLength) {
				int blockIndex = *colPointer++;
				int blockEnd = *colPointer;

				if (blockIndex == blockEnd) continue;

				//diag block
				if (blockRows[blockIndex] == column) {
					for (int i = 0; i < blockLength; i++) {
						for (int j = i; j < blockLength; j++)
							indices[column + i].emplace(column + j, index++);
					}
					blockIndex += 1;
				}

				if (blockIndex == blockEnd) continue;

				while (blockIndex < blockEnd - 1) {
					int row = blockRows[blockIndex++];
					for (int i = 0; i < blockLength; i++) {
						for (int j = 0; j < blockWidth; j++)
							indices[column + i].emplace(row + j, index++);
					}
				}

				//last block
				int row = blockRows[blockIndex];
				int lastWidth = numColumns - row;
				if (lastWidth >= blockWidth) lastWidth = blockWidth;
				for (int i = 0; i < blockLength; i++) {
					for (int j = 0; j < lastWidth; j++)
						indices[column + i].emplace(row + j, index++);
				}
			}

			for (int column = blockBoundary; column < numColumns; column++) {
				int iter = *colPointer++;
				int end = *colPointer;

				while (iter < end) 
					indices[column].emplace(blockRows[iter++], index++);
			}

			return indices;
		}

		double getNonzeroRatio() const {
			int dataCount = blockColumnOris[numBlockColumn + numRemainedColumn];
			int nonzeroCount = 0;
			for (int i = 0; i < dataCount; i++)
				if (values[i] != 0.0) nonzeroCount += 1;

			return (double)nonzeroCount / (double)dataCount;
		}

		BlockedSymSparseMatrix(const BlockedSymSparseMatrix&) = delete;
		BlockedSymSparseMatrix& operator=(const BlockedSymSparseMatrix&) = delete;
		BlockedSymSparseMatrix(BlockedSymSparseMatrix&& mat) noexcept{
			numColumns = mat.numColumns;
			numBlockColumn = mat.numBlockColumn;
			numRemainedColumn = mat.numRemainedColumn;
			values = mat.values;
			blockRows = mat.blockRows;
			blockPcol = mat.blockPcol;
			blockColumnEntries = mat.blockColumnOris;

			mat.numColumns = 0;
			mat.numBlockColumn = 0;
			mat.numRemainedColumn = 0;
			mat.values = NULL; mat.blockRows = NULL; mat.blockPcol = NULL; mat.blockColumnOris = NULL;
		}
		BlockedSymSparseMatrix& operator=(BlockedSymSparseMatrix&& mat) noexcept{
			std::swap(numColumns, mat.numColumns);
			std::swap(numBlockColumn, mat.numBlockColumn);
			std::swap(numRemainedColumn, mat.numRemainedColumn);
			std::swap(values, mat.values);
			std::swap(blockRows, mat.blockRows);
			std::swap(blockPcol, mat.blockPcol);
			std::swap(blockColumnEntries, mat.blockColumnOris);
		}

		~BlockedSymSparseMatrix(){
			freeAligned(values);
			freeAligned(blockRows);
			freeAligned(blockPcol);
			freeAligned(blockColumnOris);
		}

		int getNumColumns() const{ return numColumns; }
	private:

		int numColumns;
		int numBlockColumn;
		int numRemainedColumn;
		double *values;
		int *blockRows;
		int *blockPcol;
		int *blockColumnOris;

		static constexpr std::array<size_t, blockLength> diagIndices = getDiagIndices(std::make_index_sequence<blockLength>());

		template<int blockLength, int blockWidth>
		friend void SpMDV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat, 
			const double *src, double *dest);
		template<int blockLength, int blockWidth>
		friend void SpMSV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat,
			const std::vector<std::vector<std::pair<int, int>>>& fullIndices, 
			const SparseVector& src, double *dest);
	};
}

#endif