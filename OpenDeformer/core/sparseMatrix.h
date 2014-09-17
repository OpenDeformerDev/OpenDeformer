#if defined(_MSC_VER)
#pragma once
#endif

#include "oder.h"
#include "memory.h"

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
}

#endif