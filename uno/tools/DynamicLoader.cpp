// Copyright (c) 2026 Joris Gillis, Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DynamicLoader.hpp"

#ifndef _WIN32
#include <dlfcn.h>
#endif

namespace uno {
#ifdef _WIN32
   LibraryHandle open_library(const char* name) {
      return LoadLibraryA(name);
   }

   void* raw_symbol(LibraryHandle handle, const char* symbol) {
      return reinterpret_cast<void*>(GetProcAddress(handle, symbol));
   }
#else
   // match upstream IPOPT: resolve now, do not export the HSL symbols globally
   LibraryHandle open_library(const char* name) {
      return dlopen(name, RTLD_NOW);
   }

   void* raw_symbol(LibraryHandle handle, const char* symbol) {
      return dlsym(handle, symbol);
   }
#endif

   // Resolve a Fortran symbol trying the manglings IPOPT tries, so the runtime
   // libhsl can have been built by any compiler regardless of how Uno was:
   // base, base_, lower_, lower, UPPER_, UPPER.
   void* resolve_symbol(LibraryHandle handle, const std::string& base) {
      std::string lower = base, upper = base;
      for (char& c: lower) {
         c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      }
      for (char& c: upper) {
         c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      }
      const std::string candidates[] = {base, base + "_", lower + "_", lower, upper + "_", upper};
      for (const std::string& candidate: candidates) {
         void* symbol = raw_symbol(handle, candidate.c_str());
         if (symbol != nullptr) {
            return symbol;
         }
      }
      return nullptr;
   }
} // namespace