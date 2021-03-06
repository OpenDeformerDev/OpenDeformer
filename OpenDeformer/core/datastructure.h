#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_DATASTRUCTURE_H
#define ODER_CORE_DATASTRUCTURE_H

#include "oder.h"
#include "memory.h"
#include <type_traits>
#include <set>

namespace ODER{

	template<class T, class Compare = std::less<T>> class DisjointSets{
	private:
		template<class T> struct DisjointSet{
			T value;
			int rank;
			DisjointSet *parent;
		};
	public:
		DisjointSets(Compare comp = Compare()) :sets(SetCompare(std::move(comp))){
			sentinel = arena.Alloc();
		}
		void makeSet(T x){ 
			DisjointSet<T> *entry = arena.Alloc();
			entry->value = std::move(x); entry->rank = 0; entry->parent = entry;
			sets.insert(entry);
		}
		void Union(T x, T y){ Link(findSet(x), findSet(y));}
		bool isInSameSet(T x, T y){ return findSet(x) == findSet(y);}
	private:
		DisjointSet<T>* findSet(T x) const{
			sentinel->value = x;
			auto result = sets.find(sentinel);
			if (result == sets.end())
				return NULL;
			return localFind(*result);
		}
		void Link(DisjointSet<T>* x, DisjointSet<T>* y){
			if (x->rank > y->rank)
				y->parent = x;
			else{
				x->parent = y;
				if (x->rank == y->rank)
					y->rank++;
			}
		}
		DisjointSet<T>* localFind(DisjointSet<T>* x){
			if (x != x->parent)
				x->parent = localFind(x->parent);
			return x->parent;	
		}

		struct SetCompare{
			SetCompare(Compare comp) :comparer(std::move(comp)){}
			bool operator()(const DisjointSet<T>* left, const DisjointSet<T>* right){
				return comparer(left->value, right->value);
			}
		private:
			Compare comparer;
		};
		DisjointSet<T>* sentinel;
		std::set<DisjointSet<T>*, SetCompare> sets;
		MemoryArena<DisjointSet<T>> arena;
	};

	//not thread safe
	template<class T> class RecycledList{
	private:
		template<class T> struct ListNode{
			template<class... Args> ListNode(Args&&... vals): data(std::forward<Args>(vals)...){}
			ListNode(const ListNode&) = delete;
			ListNode& operator=(const ListNode&) = delete;
			ListNode(ListNode&&) = default;
			ListNode& operator=(ListNode&&) = default;

			~ListNode(){ data.~T(); }
			ListNode* next;
			ListNode* prev;
			T data;
		};
	public:
		RecycledList(){
			++listCount;
			node = pool->Alloc();
			node->next = node;
			node->prev = node;
		}

		RecycledList(const RecycledList &) = default;
		RecycledList& operator=(const RecycledList &) = default;
		RecycledList(RecycledList&& list) noexcept{
			node = list.node;
			list.node = NULL;
		}
		RecycledList& operator=(RecycledList&& list) noexcept{
			std::swap(node, list.node);
		}

		template<class T, class Ref, class Ptr> struct ListIterator{
			using iterator_category = std::bidirectional_iterator_tag;
			using value_type = T;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = Ref;
			using pointer = Ptr;

			ListIterator():node(nullptr){}
			ListIterator(ListNode<T>* n) :node(n){}
			ListIterator(const ListIterator&) = default;
			ListIterator(ListIterator&&) = default;
			ListIterator& operator=(const ListIterator&) = default;
			ListIterator& operator=(ListIterator&&) = default;

			bool operator==(const ListIterator& x) const{ return node == x.node; }
			bool operator!=(const ListIterator& x) const{ return node != x.node; }

			reference operator*() const{ return node->data; }
			pointer operator->() const{ return &(node->data); }

			ListIterator& operator++(){
				node = node->next;
				return *this;
			}
			ListIterator operator++(int){
				ListIterator temp = *this;
				++*this;
				return temp;
			}
			ListIterator& operator--(){
				node = node->prev;
				return *this;
			}
			ListIterator operator--(int){
				ListIterator temp = *this;
				--*this;
				return temp;
			}
			~ListIterator() = default;
			ListNode<T>* node;
		};

		using const_type = typename std::add_const<T>::type;
		using iterator = ListIterator<T, T&, T*>;
		using const_iterator = ListIterator<T, const_type&, const_type*>;

		iterator begin(){ return node->next; }
		iterator end(){ return node; }
		const_iterator cbegin() const{ return node->next; }
		const_iterator cend() const{ return node; }

		void push_back(T&& val){
			emplace_back(std::forward<T>(val));
		}
		const_iterator insert(const const_iterator& pos, T&& val){
			return emplace(pos, std::forward<T>(val));
		}
		iterator insert(const iterator& pos, T&& val){
			return emplace(pos, std::forward<T>(val));
		}
		template<class Ref, class Ptr, class... Args> ListIterator<T, Ref, Ptr> emplace(const ListIterator<T, Ref, Ptr>& pos, Args&&... val){
			ListNode<T>* postNode = pos.node;
			ListNode<T>* tmp = pool->Alloc(std::forward<Args>(val)...);
			tmp->next = postNode;
			tmp->prev = postNode->prev;
			postNode->prev->next = tmp;
			postNode->prev = tmp;

			return tmp;
		}
		template<class... Args> void emplace_back(Args&&... val){
			ListNode<T>* end = node;
			ListNode<T>* tmp = pool->Alloc(std::forward<Args>(val)...);
			tmp->next = end;
			tmp->prev = end->prev;
			end->prev->next = tmp;
			end->prev = tmp;
		}
		void clear(){
			ListNode<T>* cur = node->next;
			while (cur != node){
				ListNode<T>* tmp = cur;
				cur = cur->next;
				pool->Dealloc(tmp);
			}
			node->next = node;
			node->prev = node;
		}
		template<class Ref, class Ptr> ListIterator<T, Ref, Ptr> erase(const ListIterator<T, Ref, Ptr>& iter){
			ListNode<T>* toBeDelete = iter.node;
			ListNode<T>* prev = toBeDelete->prev;
			ListNode<T>* next = toBeDelete->next;

			prev->next = next;
			next->prev = prev;
			pool->Dealloc(toBeDelete);

			return next;
		}

		~RecycledList(){
			if (--listCount == 0)
				pool->freeAll();
		}
	private:
		ListNode<T> *node;
		static int listCount;
		static MemoryPool<ListNode<T>> *pool;
	};

	template<class T> int RecycledList<T>::listCount = 0;
	template<class T> MemoryPool<typename RecycledList<T>::ListNode<T>> *RecycledList<T>::pool = new MemoryPool<typename RecycledList<T>::ListNode<T>>();
}

#endif