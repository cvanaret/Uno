// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PythonStreamBuffer.hpp"

namespace uno {
   PythonStreamBuffer::PythonStreamBuffer(pybind11::object output_function) : output_function(output_function) {
   }

   std::streamsize PythonStreamBuffer::xsputn(const char* s, std::streamsize n) {
      output_function(std::string(s, n));
      return n;
   }

   int PythonStreamBuffer::overflow(int c) {
      if (c != EOF) {
         output_function(std::string(1, static_cast<char>(c)));
      }
      return c;
   }
} // namespace