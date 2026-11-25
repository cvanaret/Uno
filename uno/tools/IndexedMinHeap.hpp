// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INDEXEDMINHEAP_H
#define UNO_INDEXEDMINHEAP_H

#include <cstddef>
#include <limits>
#include <vector>

namespace uno {
   // binary min-heap of indices in {0, ..., capacity - 1}, ordered by the keys of an external array (not owned). Each index
   // is stored at most once and its slot in the heap is tracked, so that an index can be moved up after its key decreased
   // (decrease-key) or removed, in O(log size). The keys of the indices in the heap must only be modified through
   // decrease-key (push_or_decrease)
   class IndexedMinHeap {
   public:
      IndexedMinHeap(size_t capacity, const std::vector<double>& keys): keys(keys), slot_of(capacity, absent) {
         this->heap.reserve(capacity);
      }

      [[nodiscard]] bool empty() const { return this->heap.empty(); }
      // index with the smallest key (the heap must not be empty)
      [[nodiscard]] size_t top() const { return this->heap.front(); }
      [[nodiscard]] bool contains(size_t index) const { return this->slot_of[index] != absent; }

      // insert the index, or restore the heap order after its key decreased
      void push_or_decrease(size_t index) {
         if (this->slot_of[index] == absent) {
            this->slot_of[index] = this->heap.size();
            this->heap.push_back(index);
         }
         this->sift_up(this->slot_of[index]);
      }

      // remove and return the index with the smallest key (the heap must not be empty)
      size_t pop() {
         const size_t index = this->heap.front();
         this->remove(index);
         return index;
      }

      // remove an index of the heap
      void remove(size_t index) {
         const size_t slot = this->slot_of[index];
         const size_t last = this->heap.back();
         this->heap.pop_back();
         this->slot_of[index] = absent;
         if (last != index) {
            // the last index fills the hole, then moves up or down
            this->place(slot, last);
            this->sift_up(slot);
            this->sift_down(this->slot_of[last]);
         }
      }

      // empty the heap in O(size)
      void clear() {
         for (const size_t index: this->heap) {
            this->slot_of[index] = absent;
         }
         this->heap.clear();
      }

   private:
      static constexpr size_t absent{std::numeric_limits<size_t>::max()};
      const std::vector<double>& keys;
      std::vector<size_t> heap{}; // heap[slot]: index
      std::vector<size_t> slot_of; // slot of each index in the heap, or absent

      void place(size_t slot, size_t index) {
         this->heap[slot] = index;
         this->slot_of[index] = slot;
      }

      void sift_up(size_t slot) {
         const size_t index = this->heap[slot];
         const double key = this->keys[index];
         while (0 < slot) {
            const size_t parent = (slot - 1) / 2;
            if (this->keys[this->heap[parent]] <= key) {
               break;
            }
            this->place(slot, this->heap[parent]);
            slot = parent;
         }
         this->place(slot, index);
      }

      void sift_down(size_t slot) {
         const size_t index = this->heap[slot];
         const double key = this->keys[index];
         const size_t size = this->heap.size();
         while (true) {
            size_t child = 2 * slot + 1;
            if (size <= child) {
               break;
            }
            if (child + 1 < size && this->keys[this->heap[child + 1]] < this->keys[this->heap[child]]) {
               ++child;
            }
            if (key <= this->keys[this->heap[child]]) {
               break;
            }
            this->place(slot, this->heap[child]);
            slot = child;
         }
         this->place(slot, index);
      }
   };
} // namespace

#endif // UNO_INDEXEDMINHEAP_H
