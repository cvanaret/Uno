// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PYTHONSTREAMBUFFER_H
#define UNO_PYTHONSTREAMBUFFER_H

#include <pybind11/pybind11.h>
#include <streambuf>
#include <string>

namespace uno {
   // line-buffered stream that forwards its content to a Python stream: write() is called once per line (or on flush)
   // instead of once per inserted piece of text
   class PythonStreamBuffer : public std::streambuf {
   public:
      explicit PythonStreamBuffer(const pybind11::object& stream);
      ~PythonStreamBuffer() override;

   protected:
      std::streamsize xsputn(const char* s, std::streamsize n) override;
      int overflow(int c) override;
      int sync() override;

   private:
      pybind11::object write_function;
      pybind11::object flush_function;
      std::string pending_text{};

      void write_pending_text();
   };
} // namespace

#endif // UNO_PYTHONSTREAMBUFFER_H
