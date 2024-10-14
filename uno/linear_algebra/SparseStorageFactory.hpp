// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SPARSESTORAGEFACTORY_H
#define UNO_SPARSESTORAGEFACTORY_H

#include "SparseStorage.hpp"
#include "COOSparseStorage.hpp"
#include "CSCSparseStorage.hpp"

namespace uno {
   template <typename IndexType, typename ElementType>
   class SparseStorageFactory {
   public:
      static std::unique_ptr<SparseStorage<IndexType, ElementType>> create(const std::string& sparse_storage_type, size_t number_rows,
            size_t number_columns, size_t capacity, bool use_regularization);
   };

   template <typename IndexType, typename ElementType>
   std::unique_ptr<SparseStorage<IndexType, ElementType>> SparseStorageFactory<IndexType, ElementType>::create(const std::string& sparse_storage_type,
         size_t number_rows, size_t number_columns, size_t capacity, bool use_regularization) {
      if (sparse_storage_type == "COO") {
         return std::make_unique<COOSparseStorage<IndexType, ElementType>>(number_rows, number_columns, capacity, use_regularization);
      }
      else if (sparse_storage_type == "CSC") {
         return std::make_unique<CSCSparseStorage<IndexType, ElementType>>(number_rows, number_columns, capacity, use_regularization);
      }
      throw std::invalid_argument("Sparse storage " + sparse_storage_type + " unknown");
   }
} // namespace

#endif // UNO_SPARSESTORAGEFACTORY_H