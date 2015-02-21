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
			blockEntries(numBlockColumn), remainedEntries(NULL){

			if (numRemainedColumn > 0){
				int numRemainedEntries = numRemainedColumn*numRemainedColumn;
				remainedEntries = entryMem.Alloc(numRemainedEntries);
			}

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
				int subRow = row - columnStart;
				remainedEntries[subColumn * numRemainedColumn + subRow] += data;
			}
		}

	private:
		int numColumn;
		int numBlockColumn;
		int numRemainedColumn;
		vector<map<int, double*>> blockEntries;
		double *remainedEntries;

		double diagIndices[blockLength];
		MemoryArena<double> entryMem;

		friend class BlockedSymSparseMatrix < blockLength, blockWidth > ;
	};

	template<int blockLength, int blockWidth> class BlockedSymSparseMatrix{
	public:
		BlockedSymSparseMatrix(const BlockedSymSparseMatrixAssembler<blockLength, blockWidth>& assembler);
	private:
		int numColumns;
		int numBlockColumn;
		int numRemainedColumn;
		double *values;
		double *blockRows;
		double *blockPcol;
	};
}

#endif