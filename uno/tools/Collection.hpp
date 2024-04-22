// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_COLLECTION_H
#define UNO_COLLECTION_H

#include <functional>

// iterable collection
template <typename ElementType>
class Collection {
public:
   using value_type = ElementType;

   explicit Collection() = default;
   virtual ~Collection() = default;

   virtual void for_each(const std::function<void (size_t /*index*/, ElementType /*element*/)>& f) const = 0;
   [[nodiscard]] virtual size_t size() const = 0;
   [[nodiscard]] bool is_empty() const;
};

template <typename ElementType>
inline bool Collection<ElementType>::is_empty() const {
   return (this->size() == 0);
}

#endif // UNO_COLLECTION_H