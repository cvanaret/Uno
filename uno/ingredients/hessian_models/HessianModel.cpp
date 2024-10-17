// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"

namespace uno {
   HessianModel::HessianModel(size_t dimension, size_t maximum_number_nonzeros, const std::string& sparse_format, bool use_regularization) :
         hessian(dimension, maximum_number_nonzeros, use_regularization, sparse_format) {
   }

   HessianModel::~HessianModel() { }
} // namespace