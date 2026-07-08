// Copyright (c) 2026 Joris Gillis, Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifdef _WIN32
#include <windows.h>
#endif

#include <string>

namespace uno {
#ifdef _WIN32
   using LibraryHandle = HMODULE;
#else
   using LibraryHandle = void*;
#endif

   LibraryHandle open_library(const char* name);
   void* raw_symbol(LibraryHandle handle, const char* symbol);
   void* resolve_symbol(LibraryHandle handle, const std::string& base);

   template <typename FunctionPointer>
   void resolve(LibraryHandle handle, FunctionPointer& function_pointer, const std::string& base) {
      function_pointer = reinterpret_cast<FunctionPointer>(resolve_symbol(handle, base));
   }
} // namespace