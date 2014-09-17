#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_DATASTRUCTURE_H
#define ODER_CORE_DATASTRUCTURE_H

#include "oder.h"
#include "memory.h"

namespace ODER{
	static size_t primeArray[] = { 17, 37, 79, 163, 331, 673, 1361, 2729,
		5471, 10949, 21911, 43853, 87719, 175447,
		350899, 701819, 1403641, 2807303, 5614657,
		11229331, 22458671, 44917381, 89834777, 179669557,
		359339171, 718678369, 1437356741, 2147483647 };

	enum HashEntryType{ Legitimate, Empty, Deleted };

	template<class T> struct HashEntry{
		HashEntry() :type(Empty){}
		HashEntry(const T& e) :element(e), type(Legitimate){}
		T element;
		HashEntryType type;
	};

	template<class Key, class Value> struct HashPairEntry{
		HashPairEntry() :type(Empty){}
		HashPairEntry(const Key& k, const Value& v)
			:key(k), value(v), type(Legitimate){}
		Key key;
		Value value;
		HashEntryType type;
	};

	template<class T, class Hasher> class HashSet{
	public:
		HashSet(int sizeIndex = 7, float rehashFactor = 0.7f){
			factor = rehashFactor;
			load = 0;
			primeIndex = sizeIndex;
			tableSize = primeArray[primeIndex];
			table = new HashEntry<T>[tableSize];
		}
		~HashSet(){
			delete[] table;
		}
		void insert(const T& key){
			size_t collison = 0;
			size_t pos = hasher(key) % tableSize;

			while (table[pos].type != Empty && table[pos].type != Deleted
				&& table[pos].element != key){
				pos += (2 * ++collison - 1);
				if (pos >= tableSize)
					pos -= tableSize;
			}

			if (table[pos].type != Legitimate){
				table[pos] = HashEntry<T>(key);
				load++;
			}
			if (float(load) / float(tableSize) > factor)
				rehash();
		}

		HashEntry<T>* find(const T& key){
			size_t pos = findPos(key);
			if (table[pos].type == Legitimate){
				return table + pos;
			}
			return table + tableSize;
		}

		const HashEntry<T>* find(const T& key) const{
			size_t pos = findPos(key);
			if (table[pos].type == Legitimate){
				return table + pos;
			}
			return table + tableSize;
		}

		size_t findPos(const T& key) const{
			size_t collison = 0;
			size_t pos = hasher(key) % tableSize;

			while (table[pos].type != Empty && table[pos].element != key){
				pos += (2 * ++collison - 1);
				if (pos >= tableSize)
					pos -= tableSize;
			}
			return pos;
		}
		void erase(const T& key){
			size_t pos = findPos(key);
			if (table[pos].type == Legitimate){
				table[pos].type = Deleted;
				load--;
			}
		}
		void rehash(){
			int oldSize = tableSize;
			tableSize = primeArray[++primeIndex];
			HashEntry<T> *oldTable = table;
			table = new HashEntry<T>[tableSize];
			for (int i = 0; i < oldSize; i++){
				if (oldTable[i].type == Legitimate){
					const T& key = oldTable[i].element;
					size_t collison = 0;
					size_t pos = hasher(key) % tableSize;

					while (table[pos].type != Empty && table[pos].element != key){
						pos += (2 * ++collison - 1);
						if (pos >= tableSize)
							pos -= tableSize;
					}

					table[pos] = HashEntry<T>(key);
				}
			}
			delete[] oldTable;
		}
		int size() const{
			return load;
		}
		const HashEntry<T>* begin() const{
			return table;
		}
		const HashEntry<T>* end() const{
			return table + tableSize;
		}

	private:
		HashEntry<T> *table;
		size_t tableSize;
		int primeIndex;
		size_t load;
		float factor;

		Hasher hasher;
	};

	template<class Key, class Value, class Hasher> class HashMap{
	public:
		HashMap(int sizeIndex = 8, float rehashFactor = 0.7f){
			factor = rehashFactor;
			load = 0;
			primeIndex = sizeIndex;
			tableSize = primeArray[primeIndex];
			table = new HashPairEntry<Key, Value>[tableSize];
		}
		~HashMap(){
			delete[] table;
		}
		void insert(const Key& k, const Value& v){
			size_t collsion = 0;
			size_t pos = hasher(k) % tableSize;

			while (table[pos].type != Empty && table[pos].type != Deleted
				&&table[pos].key != k){
				pos += (2 * ++collsion - 1);
				if (pos >= tableSize)
					pos -= tableSize;
			}

			if (table[pos].type != Legitimate){
				table[pos] = HashPairEntry<Key, Value>(k, v);
				load++;
			}

			if (float(load) / float(tableSize) > factor)
				rehash();
		}
		size_t findPos(const Key& k) const{
			size_t collison = 0;
			size_t pos = hasher(k) % tableSize;

			while (table[pos].type != Empty && table[pos].key != k){
				pos += (2 * ++collison - 1);
				if (pos >= tableSize)
					pos -= tableSize;
			}
			return pos;
		}
		HashPairEntry<Key, Value> *find(const Key& k){
			size_t pos = findPos(k);
			if (table[pos].type == Legitimate)
				return table + pos;
			return table + tableSize;
		}
		const HashPairEntry<Key, Value> *find(const Key& k) const{
			size_t pos = findPos(k);
			if (table[pos].type == Legitimate)
				return table + pos;
			return table + tableSize;
		}

		Value& operator[](const Key& k){
			size_t pos = findPos(k);
			if (table[pos].type != Legitimate){
				table[pos].type = Legitimate;
				table[pos].key = k;
				load++;
			}
			return table[pos].value;
		}
		const Value& operator[](const Key& k) const{
			size_t pos = findPos(k);
			return table[pos].value;
		}
		void erase(const Key& k){
			size_t pos = findPos(k);
			if (table[pos].type == Legitimate){
				table[pos].type = Deleted;
				load--;
			}
		}
		void rehash(){
			int oldSize = tableSize;
			tableSize = primeArray[++primeIndex];
			HashPairEntry<Key, Value> *oldTable = table;
			table = new HashPairEntry<Key, Value>[tableSize];
			for (int i = 0; i < oldSize; i++){
				if (oldTable[i].type == Legitimate){
					const Key& key = oldTable[i].key;
					const Value& value = oldTable[i].value;
					size_t collison = 0;
					size_t pos = hasher(key) % tableSize;

					while (table[pos].type != Empty && table[pos].key != key){
						pos += (2 * ++collison - 1);
						if (pos >= tableSize)
							pos -= tableSize;
					}

					table[pos] = HashPairEntry<Key, Value>(key, value);
				}
			}
			delete[] oldTable;
		}
		int size() const{
			return load;
		}
		const HashPairEntry<Key, Value>* begin() const{
			return table;
		}
		const HashPairEntry<Key, Value>* end() const{
			return table + tableSize;
		}

	private:
		HashPairEntry<Key, Value> *table;
		size_t tableSize;
		int primeIndex;
		size_t load;
		float factor;

		Hasher hasher;
	};
}

#endif