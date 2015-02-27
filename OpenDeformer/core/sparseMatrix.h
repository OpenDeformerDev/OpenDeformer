#if defined(_MSC_VER)
#pragma once
#endif

#include "oder.h"
#include "memory.h"
#include "latool.h"

#ifndef ODER_CORE_SPARSEMARTIX_H
#define ODER_CORE_SPARSEMARTIX_H

namespace ODER{
	class SparseMatrixAssembler{
	public:
		SparseMatrixAssembler(int nc) :numColumns(nc){
			map<int, double> empty;
			for (int i = 0; i < numColumns; i++)
				rowEntries.push_back(empty);
		}
		void addEntry(int row, int column, double data);
		void removeRow(int row);
		void removeColumn(int column);
		void removeColumnAndRow(int num);

		int numColumns;
		vector<map<int, double>> rowEntries;
	};

	class SparseMatrix{
	public:
		SparseMatrix(const SparseMatrixAssembler& assembler);
		~SparseMatrix();
		int getNumColumns() const;
		void Print() const;
	private:
		int numColumns;
		//CSC format
		double* value;
		int* rows;//row number of each entry
		int* pcol;//the index of the start of each column and the end

		friend class EigenSlover;
	};

	template<int blockLength, int blockWidth> class BlockedSymSparseMatrixAssembler{
	public:
		BlockedSymSparseMatrixAssembler(int columnCount) 
			: numColumn(columnCount), numBlockColumn(numColumn / blockLength), numRemainedColumn(numColumn - numBlockColumn*blockLength),
			blockEntries(numBlockColumn), remainedEntryCount(NULL){

			if (numRemainedColumn > 0) remainedEntryCount = new int[numRemainedColumn];

			diagIndices[0] = 0;
			for (int i = 0; i < blockLength - 1; i++)
				diagIndices[i + 1] = diagIndices[i] + (blockLength - i);

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
					auto found = blockEntries[blockCloumn].find(0);
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
							vals[subColumn*blockWidth + subRow] = data;
							blockEntries[blockCloumn].insert(std::pair<int, double*>(blockStartRow, vals));
						}
					}
				}
			}
			else{
				int subRow = numColumn - row;
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
		vector<map<int, double*>> blockEntries;
		map<int, double> remainedEntries;
		int *remainedEntryCount;

		double diagIndices[blockLength];
		MemoryArena<double> entryMem;

		friend class BlockedSymSparseMatrix < blockLength, blockWidth > ;
	};

	template<int blockLength, int blockWidth> class BlockedSymSparseMatrix{
	public:
		BlockedSymSparseMatrix(const BlockedSymSparseMatrixAssembler<blockLength, blockWidth>& assembler){
			numColumns = assembler.numColumn;
			numBlockColumn = assembler.numBlockColumn;
			numRemainedColumn = assembler.numRemainedColumn;

			blockPcol = allocAligned<int>(numBlockColumn + numRemainedColumn + 1);

			constexpr int regularSize = blockLength*blockWidth;
			constexpr int diagSize = (blockLength + 1) * blockLength / 2;

			blockPcol[0] = 0;
			int dataCount = 0;
			for (int i = 0; i < numBlockColumn; i++){
				//counting blocks
				int blockCount = assembler.blockEntries[i].size();
				blockPcol[i + 1] = blockPcol[i] + blockCount;
				//counting datas
				if (blockCount != 0){
					int blockDataCount = blockCount*regularSize;
					auto start = assembler.blockEntries[i].cbegin();
					if (start->first == i*blockLength)
						blockDataCount -= (regularSize - diagSize);
					auto last = --(assembler.blockEntries[i].cend());
					if (numColumns - last->first < blockWidth)
						blockDataCount -= (last->first + blockWidth - numColumns)*blockLength;
					dataCount += blockDataCount;
				}
			}

			for (int i = 0; i < numRemainedColumn; i++)
				blockPcol[i + numBlockColumn + 1] = blockPcol[i + numBlockColumn] + assembler.remainedEntryCount[i];

			dataCount += assembler.remainedEntries.size();


			values = allocAligned<double>(dataCount);
			blockRows = allocAligned<int>(blockPcol[numBlockColumn + numRemainedColumn + 1]);

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

		~BlockedSymSparseMatrix(){
			freeAligned(values);
			freeAligned(blockRows);
			freeAligned(blockPcol);
		}

		int getNumColumns() const{ return numColumns; }
	private:
		int numColumns;
		int numBlockColumn;
		int numRemainedColumn;
		double *values;
		int *blockRows;
		int *blockPcol;

		template<int blockLength, int blockWidth> 
		friend void SpMV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat, 
			const DenseVector& src, DenseVector& dest);
	};
}

#endif