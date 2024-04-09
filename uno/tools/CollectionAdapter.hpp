// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COLLECTIONADAPTER_H
#define UNO_COLLECTIONADAPTER_H

#include "Collection.hpp"

// Adapter to adapt standard iterators (e.g. std::vector) to the Collection interface
template <typename Iterator>
class CollectionAdapter: public Collection<typename std::remove_reference<Iterator>::type::value_type> {
public:
   explicit CollectionAdapter(Iterator&& iterator);
   void for_each(const std::function<void(size_t /*index*/, typename CollectionAdapter::value_type /*element*/)>& f) const override;
   [[nodiscard]] size_t size() const override;

protected:
   Iterator iterator;
};

template <typename Iterator>
CollectionAdapter<Iterator>::CollectionAdapter(Iterator&& iterator): iterator(std::forward<Iterator>(iterator)) {
}

template <typename Iterator>
void CollectionAdapter<Iterator>::for_each(const std::function<void(size_t /*index*/, typename CollectionAdapter::value_type /*element*/)>& f) const {
   size_t index = 0;
   for (auto element: this->iterator) {
      f(index, element);
      index++;
   }
}

template <typename Iterator>
size_t CollectionAdapter<Iterator>::size() const {
   return this->iterator.size();
}

// free function
template <typename Iterator>
CollectionAdapter<Iterator> adapt(Iterator&& iterator) {
   return CollectionAdapter(std::forward<Iterator>(iterator));
}

#endif // UNO_COLLECTIONADAPTER_H
