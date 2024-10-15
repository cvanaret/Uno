// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_CHAINCOLLECTION_H
#define UNO_CHAINCOLLECTION_H

#include "Collection.hpp"

namespace uno {
   // Concatenation based on https://stackoverflow.com/questions/24175279/how-to-store-either-rvalue-or-lvalue-references-in-template
   // Aggregates two collections, not necessarily of the same type. However they must have the same value_type (the type of the elements)
   template <typename Collection1, typename Collection2>
   class Concatenation: public Collection<typename std::remove_reference_t<Collection1>::value_type> {
   public:
      Concatenation(Collection1&& collection1, Collection2&& collection2);
      [[nodiscard]] size_t size() const override;

      typename Concatenation::value_type dereference_iterator(size_t index) const override {
         const size_t collection1_size = this->collection1.size();
         if (index < collection1_size) {
            return collection1.dereference_iterator(index);
         }
         else {
            return collection2.dereference_iterator(index - collection1_size);
         }
      }

      void increment_iterator(size_t& index) const override {
         index++;
      }

   protected:
      const Collection1 collection1;
      const Collection2 collection2;
   };

   template <typename Collection1, typename Collection2>
   Concatenation<Collection1, Collection2>::Concatenation(Collection1&& collection1, Collection2&& collection2):
         Collection<typename Concatenation::value_type>(),
         collection1(std::forward<Collection1>(collection1)), collection2(std::forward<Collection2>(collection2)) {
      static_assert(std::is_same<typename std::remove_reference_t<Collection1>::value_type,
            typename std::remove_reference_t<Collection2>::value_type>::value, "The iterators should contain the same type");
   }

   template <typename Collection1, typename Collection2>
   size_t Concatenation<Collection1, Collection2>::size() const {
      return this->collection1.size() + this->collection2.size();
   }

   // free function
   template <typename Collection1, typename Collection2>
   Concatenation<Collection1, Collection2> concatenate(Collection1&& collection1, Collection2&& collection2) {
      return {std::forward<Collection1>(collection1), std::forward<Collection2>(collection2)};
   }
} // namespace

#endif // UNO_CHAINCOLLECTION_H
