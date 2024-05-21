// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COLLECTION_H
#define UNO_COLLECTION_H

#include <functional>

// iterable collection
template <typename ElementType>
class Collection {
public:
   // https://internalpointers.com/post/writing-custom-iterators-modern-cpp
   class iterator {
   public:
      iterator(const Collection& collection, size_t index): collection(collection), index(index) { }

      [[nodiscard]] std::pair<size_t, ElementType> operator*() const {
         return this->collection.dereference_iterator(this->index, 0);
      }

      iterator& operator++() {
         this->collection.increment_iterator(this->index);
         return *this;
      }

      friend bool operator!=(const iterator& a, const iterator& b) {
         return &a.collection != &b.collection || a.index != b.index;
      }

   protected:
      const Collection& collection;
      size_t index;
   };

   using value_type = ElementType;

   explicit Collection() = default;
   virtual ~Collection() = default;

   [[nodiscard]] virtual size_t size() const = 0;
   [[nodiscard]] bool is_empty() const;

   [[nodiscard]] iterator begin() const { return iterator(*this, 0); }
   [[nodiscard]] iterator end() const { return iterator(*this, this->size()); }

   [[nodiscard]] virtual std::pair<size_t, ElementType> dereference_iterator(size_t index, size_t offset) const = 0;
   virtual void increment_iterator(size_t& index) const = 0;
};

template <typename ElementType>
inline bool Collection<ElementType>::is_empty() const {
   return (this->size() == 0);
}

#endif // UNO_COLLECTION_H