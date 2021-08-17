#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H

#include <cassert>
#include <unordered_map>

using SparseVector = std::unordered_map<unsigned int, double>;

/*
class SparseVector {
public:
   SparseVector(int capacity) {
      this->indices.reserve(capacity);
      this->elements.reserve(capacity);
      //std::cout << " *** Initialized with capacity " << capacity << "\n";
   }

   SparseVector() {
   }

   void reserve(size_t capacity) {
      this->indices.reserve(capacity);
      this->elements.reserve(capacity);
   }

   size_t size() const {
      return this->number_nonzeros;
   }

   double& operator[](size_t index) {
      size_t i = 0;
      while (i < this->number_nonzeros) {
         if (this->indices[i] == index) {
            return this->elements[i];
         }
         i++;
      }
      assert(this->number_nonzeros <= this->indices.capacity());
      this->indices[this->number_nonzeros] = index;
      double& result = this->elements[this->number_nonzeros];
      this->number_nonzeros++;
      return result;
   }

   const double& operator[](size_t index) const {
      size_t i = 0;
      while (i < this->number_nonzeros) {
         if (this->indices[i] == index) {
            return this->elements[i];
         }
         i++;
      }
      assert(false && "Entry not found");
   }

   double at(size_t i) const {
      return (*this)[i];
   }

   void erase(size_t index) {
      size_t i = 0;
      while (i < this->number_nonzeros) {
         if (this->indices[i] == index) {
            const int nnz = this->number_nonzeros;
            // swap with the last elements
            this->indices[i] = this->indices[nnz - 1];
            this->elements[i] = this->elements[nnz - 1];
            this->number_nonzeros--;
            return;
         }
         i++;
      }
   }

   void clear() {
      this->number_nonzeros = 0;
   }

   class iterator {
   public:
      iterator(const SparseVector& v, int i) : vector(&v) {
         current_index = i;
      }
      bool operator!=(const iterator x) const {
         return current_index != x.current_index;
      };
      iterator& operator++() {
         ++current_index;
         return *this;
      };
      std::pair<size_t, double> operator*() const {
         return std::pair(this->vector->indices[this->current_index], this->vector->elements[this->current_index]);
      };

   private:
      const SparseVector* vector;
      int current_index;
   };

   iterator begin() {
      return iterator(*this, 0);
   };
   iterator end() {
      return iterator(*this, this->number_nonzeros);
   };
   iterator begin() const {
      return iterator(*this, 0);
   };
   iterator end() const {
      return iterator(*this, this->number_nonzeros);
   };

//protected:
   std::vector<size_t> indices;
   std::vector<double> elements;
   size_t number_nonzeros{0};
};
*/

#endif // SPARSEVECTOR_H
