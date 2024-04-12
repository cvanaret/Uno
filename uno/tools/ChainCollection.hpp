// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CHAINCOLLECTION_H
#define UNO_CHAINCOLLECTION_H

#include "Collection.hpp"

// Chain collection based on https://stackoverflow.com/questions/24175279/how-to-store-either-rvalue-or-lvalue-references-in-template
// Aggregates two collections, not necessarily of the same type.
template <typename Collection1, typename Collection2>
class ChainCollection: public Collection<typename std::remove_reference<Collection1>::type::value_type> {
public:
   ChainCollection(Collection1&& collection1, Collection2&& collection2);
   void for_each(const std::function<void(size_t /*index*/, typename ChainCollection::value_type /*element*/)>& f) const override;
   [[nodiscard]] size_t size() const override;

protected:
   Collection1 collection1;
   Collection2 collection2;
};

template <typename Collection1, typename Collection2>
ChainCollection<Collection1, Collection2>::ChainCollection(Collection1&& collection1, Collection2&& collection2):
      Collection<typename ChainCollection::value_type>(),
      collection1(std::forward<Collection1>(collection1)), collection2(std::forward<Collection2>(collection2)) {
   static_assert(std::is_same<typename std::remove_reference<Collection1>::type::value_type,
         typename std::remove_reference<Collection2>::type::value_type>::value, "The iterators should contain the same type");
}

template <typename Collection1, typename Collection2>
void ChainCollection<Collection1, Collection2>::for_each(const std::function<void(size_t /*index*/, typename ChainCollection::value_type /*element*/)>& f) const {
   this->collection1.for_each(f);
   this->collection2.for_each(f);
}

template <typename Collection1, typename Collection2>
size_t ChainCollection<Collection1, Collection2>::size() const {
   return this->collection1.size() + this->collection2.size();
}

// free function
template <typename Collection1, typename Collection2>
ChainCollection<Collection1, Collection2> concatenate(Collection1&& collection1, Collection2&& collection2) {
   return ChainCollection(std::forward<Collection1>(collection1), std::forward<Collection2>(collection2));
}

#endif // UNO_CHAINCOLLECTION_H