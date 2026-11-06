// Copyright (c) 2026 Joris Gillis, Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "DynamicLoader.hpp"
#include <cstring>
#include <cctype>
#include <cstdlib>
#include <string>
#ifndef _WIN32
#include <dlfcn.h>
#endif

namespace uno {
#ifdef _WIN32
   // with a path, LOAD_WITH_ALTERED_SEARCH_PATH searches the dependencies of the library in its own directory
   // (e.g. libkrylov.dll and its bundled Julia runtime DLLs). A bare name keeps the standard search order
   LibraryHandle open_library(const char* name) {
      const bool has_path = (std::strchr(name, '/') != nullptr || std::strchr(name, '\\') != nullptr);
      if (!has_path) {
         return LoadLibraryA(name);
      }
      // the altered search path is derived from the path: use backslashes only (e.g. MSYS2 mixed-separator paths)
      std::string path(name);
      for (char& character: path) {
         if (character == '/') {
            character = '\\';
         }
      }
      return LoadLibraryExA(path.c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
   }

   void close_library(LibraryHandle handle) {
      FreeLibrary(handle);
   }

   void* raw_symbol(LibraryHandle handle, const char* symbol) {
      return reinterpret_cast<void*>(GetProcAddress(handle, symbol));
   }
#else
   // resolve now, do not export the HSL symbols globally. RTLD_NODELETE: the library stays mapped after
   // close_library. Unmapping it would orphan the memory that it (or its runtime, e.g. gfortran or OpenMP) allocated
   // when it was loaded, and unloading such runtimes from a running process is unsafe. close_library only releases
   // the handle, and a later open_library of the same library reuses the mapped image
   LibraryHandle open_library(const char* name) {
      return dlopen(name, RTLD_NOW | RTLD_NODELETE);
   }

   void close_library(LibraryHandle handle) {
      dlclose(handle);
   }

   void* raw_symbol(LibraryHandle handle, const char* symbol) {
      return dlsym(handle, symbol);
   }
#endif

   // resolve a Fortran symbol by trying different manglings (the runtime libhsl may have been built by any compiler):
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

   // resolve the effective library name (explicit request > env > platform default)
   std::string resolve_library_name(const std::string& user_library_name, const char* env_library_name, const char* default_library) {
      if (!user_library_name.empty()) {
         return user_library_name;
      }
      // user_library_name empty
      const char* env_path = std::getenv(env_library_name);
      if (env_path != nullptr && 0 < strlen(env_path)) {
         return env_path;
      }
      return default_library;
   }
} // namespace