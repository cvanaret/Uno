// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <stdexcept>
#include "FilterFactory.hpp"
#include "Filter.hpp"
#include "options/Options.hpp"

namespace uno {
   // FilterFactory class
   std::unique_ptr<Filter> FilterFactory::create(const Options& options) {
      const std::string& filter_type = options.get_string("filter_type");
      if (filter_type == "standard") {
         return std::make_unique<Filter>(options);
      }
      else {
         throw std::invalid_argument("Filter type " + filter_type + " does not exist");
      }
   }
} // namespace
