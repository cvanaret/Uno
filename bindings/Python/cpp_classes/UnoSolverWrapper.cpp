// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "UnoSolverWrapper.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "tools/Logger.hpp"

namespace uno {
   UnoSolverWrapper::UnoSolverWrapper() {
      DefaultOptions::load(this->options);
   }

   void UnoSolverWrapper::optimize(const PythonUserModel& user_model) {
      const PythonModel model{user_model};
      Logger::set_logger(this->options.get_string("logger"));
      Result result = this->uno_solver.solve(model, this->options);
   }
} // namespace