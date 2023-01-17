// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_FILTERFACTORY_H
#define UNO_FILTERFACTORY_H

#include "Filter.hpp"

class FilterFactory {
public:
   static std::unique_ptr<Filter> create(const Options& options);
};

#endif // UNO_FILTERFACTORY_H