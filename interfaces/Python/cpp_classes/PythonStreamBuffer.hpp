// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PYTHONSTREAMBUFFER_H
#define UNO_PYTHONSTREAMBUFFER_H

#include <pybind11/pybind11.h>
#include <streambuf>

namespace uno {
   class PythonStreamBuffer : public std::streambuf {
   public:
      explicit PythonStreamBuffer(pybind11::object stream);

   protected:
      std::streamsize xsputn(const char* s, std::streamsize n) override;
      int overflow(int c) override;
      int sync() override;

   private:
      pybind11::object write_function;
      pybind11::object flush_function;
   };
} // namespace

#endif // UNO_PYTHONSTREAMBUFFER_H