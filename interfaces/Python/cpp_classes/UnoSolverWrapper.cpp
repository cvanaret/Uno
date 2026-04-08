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

   void UnoSolverWrapper::set_logger_stream(py::object python_stream) {
      this->stream = python_stream;  // keep Python object alive
      this->stream_buffer = std::make_unique<PythonStreamBuffer>(python_stream);
      this->ostream = std::make_unique<std::ostream>(this->stream_buffer.get());
      uno::Logger::set_stream(*this->ostream);
   }

   Result UnoSolverWrapper::optimize(const PythonUserModel& user_model) {
      const PythonModel model{user_model};
      Logger::set_logger(this->options.get_string("logger"));
      return this->uno_solver.solve(model, this->options);
   }
} // namespace