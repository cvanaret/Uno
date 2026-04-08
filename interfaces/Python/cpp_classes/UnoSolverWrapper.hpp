// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOSOLVERWRAPPER_H
#define UNO_UNOSOLVERWRAPPER_H

#include "Uno.hpp"
#include "options/Options.hpp"
#include "PythonModel.hpp"
#include "PythonStreamBuffer.hpp"

namespace uno {
   class UnoSolverWrapper {
   public:
      Uno uno_solver{};
      Options options{};

      UnoSolverWrapper();

      void set_logger_stream(py::object py_stream);
      [[nodiscard]] Result optimize(const PythonUserModel& user_model);

   private:
      py::object stream;
      std::unique_ptr<PythonStreamBuffer> stream_buffer;
      std::unique_ptr<std::ostream> ostream;
   };
} // namespace

#endif // UNO_UNOSOLVERWRAPPER_H