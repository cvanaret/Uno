// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COLLECTIONADAPTER_H
#define UNO_COLLECTIONADAPTER_H

#include "Collection.hpp"

namespace uno {
   // Adapter to adapt standard Arrays (e.g. std::vector) to the Collection interface
   template <typename Array>
   class CollectionAdapter: public Collection<typename std::remove_reference_t<Array>::value_type> {
   public:
      explicit CollectionAdapter(const Array& array);
      [[nodiscard]] size_t size() const override;

      [[nodiscard]] typename CollectionAdapter::value_type dereference_iterator(size_t index) const override;
      void increment_iterator(size_t& index) const override;

   protected:
      const Array& array;
   };

   template <typename Array>
   CollectionAdapter<Array>::CollectionAdapter(const Array& array):
      Collection<typename std::remove_reference_t<Array>::value_type>(), array(array) {
   }

   template <typename Array>
   size_t CollectionAdapter<Array>::size() const {
      return this->array.size();
   }

   template <typename Array>
   typename CollectionAdapter<Array>::value_type CollectionAdapter<Array>::dereference_iterator(size_t index) const {
      return this->array[index];
   }

   template <typename Array>
   void CollectionAdapter<Array>::increment_iterator(size_t& index) const {
      index++;
   }

   template <typename Array>
   CollectionAdapter<Array> adapt(const Array& array) {
      return CollectionAdapter<Array>{array};
   }
} // namespace

#endif // UNO_COLLECTIONADAPTER_H
