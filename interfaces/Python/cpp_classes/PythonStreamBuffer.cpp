// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cstring>
#include "PythonStreamBuffer.hpp"

namespace uno {
   PythonStreamBuffer::PythonStreamBuffer(const pybind11::object& stream):
         write_function(stream.attr("write")),
         flush_function(stream.attr("flush")) {
   }

   PythonStreamBuffer::~PythonStreamBuffer() {
      try {
         this->write_pending_text();
      }
      catch (...) {
         // the Python stream may already be closed: drop the remaining output
      }
   }

   std::streamsize PythonStreamBuffer::xsputn(const char* s, std::streamsize n) {
      this->pending_text.append(s, static_cast<size_t>(n));
      if (std::memchr(s, '\n', static_cast<size_t>(n)) != nullptr) {
         this->write_pending_text();
      }
      return n;
   }

   int PythonStreamBuffer::overflow(int c) {
      if (c != traits_type::eof()) {
         this->pending_text.push_back(traits_type::to_char_type(c));
         if (c == '\n') {
            this->write_pending_text();
         }
      }
      return traits_type::not_eof(c);
   }

   int PythonStreamBuffer::sync() {
      this->write_pending_text();
      pybind11::gil_scoped_acquire gil;
      this->flush_function();
      return 0;
   }

   void PythonStreamBuffer::write_pending_text() {
      if (!this->pending_text.empty()) {
         pybind11::gil_scoped_acquire gil;
         this->write_function(this->pending_text);
         this->pending_text.clear();
      }
   }
} // namespace
