// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PythonStreamBuffer.hpp"

namespace uno {
   PythonStreamBuffer::PythonStreamBuffer(pybind11::object stream):
         write_function(stream.attr("write")),
         flush_function(stream.attr("flush")) {
   }

   std::streamsize PythonStreamBuffer::xsputn(const char* s, std::streamsize n) {
      this->write_function(std::string(s, n));
      return n;
   }

   int PythonStreamBuffer::overflow(int c) {
      if (c != EOF) {
         this->write_function(std::string(1, static_cast<char>(c)));
      }
      return c;
   }

   int PythonStreamBuffer::sync() {
      this->flush_function();
      return 0;
   }
} // namespace