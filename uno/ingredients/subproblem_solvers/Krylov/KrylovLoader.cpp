// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cmath>
#include <exception>
#include <string>
#include <utility>
#include "KrylovLoader.hpp"
#include "KrylovSolver.hpp"
#include "tools/DynamicLoader.hpp"
#include "tools/Logger.hpp"
#ifndef _WIN32
#include <fcntl.h>
#endif

#if defined(_WIN32)
#define UNO_KRYLOV_DEFAULT_LIBRARY "libkrylov.dll"
#elif defined(__APPLE__)
#define UNO_KRYLOV_DEFAULT_LIBRARY "libkrylov.dylib"
#else
#define UNO_KRYLOV_DEFAULT_LIBRARY "libkrylov.so"
#endif

namespace uno {
   namespace {
      KrylovLibrary library{};
      LibraryHandle handle = nullptr;
      std::string loaded_name{};
      std::string default_name{UNO_KRYLOV_DEFAULT_LIBRARY};
      int functional_check{-1}; // -1: not checked, 0: failed, 1: passed
      std::string error_message{}; // why libkrylov is not loaded/functional

      void report(const std::string& message, bool warning) {
         error_message = message;
         if (warning) {
            WARNING << "Uno: " << message << '\n';
         }
         else {
            DEBUG << "Uno: " << message << '\n';
         }
      }

      template <typename FunctionPointer>
      void resolve_c(FunctionPointer& function_pointer, const char* symbol) {
         // C symbols: no Fortran mangling
         function_pointer = reinterpret_cast<FunctionPointer>(raw_symbol(handle, symbol));
      }

#ifdef _WIN32
      // the Julia runtime loads some of its DLLs by name (e.g. libpcre2-8.dll) with the standard search order, which
      // does not include the directory of libkrylov.dll: prepend this directory to the PATH of the process
      void add_library_directory_to_path(const std::string& library_name) {
         std::string directory = library_name;
         for (char& character: directory) {
            if (character == '/') {
               character = '\\';
            }
         }
         const size_t separator = directory.find_last_of('\\');
         if (separator == std::string::npos) {
            return;
         }
         directory.resize(separator);
         const DWORD length = GetEnvironmentVariableA("PATH", nullptr, 0);
         std::string path(length, '\0');
         if (0 < length) {
            path.resize(GetEnvironmentVariableA("PATH", path.data(), length));
         }
         if ((";" + path + ";").find(";" + directory + ";") == std::string::npos) {
            const std::string new_path = path.empty() ? directory : directory + ";" + path;
            SetEnvironmentVariableA("PATH", new_path.c_str());
         }
      }
#endif

      // the Julia runtime is initialized by the first call into libkrylov, and its I/O layer (libuv) then sets
      // O_NONBLOCK on the standard streams. This flag belongs to the open file description, which is shared with the
      // parent process: a shell reading commands from the same pipe gets EAGAIN and exits, and writes to a full pipe
      // fail. Initialize the runtime with a harmless call and restore the original flags
      void initialize_julia_runtime() {
#ifndef _WIN32
         int saved_flags[3];
         for (int fd = 0; fd < 3; ++fd) {
            saved_flags[fd] = fcntl(fd, F_GETFL);
         }
#endif
         (void) library.default_workspace_options();
#ifndef _WIN32
         for (int fd = 0; fd < 3; ++fd) {
            if (saved_flags[fd] != -1 && fcntl(fd, F_GETFL) != saved_flags[fd]) {
               fcntl(fd, F_SETFL, saved_flags[fd]);
            }
         }
#endif
      }

      bool common_symbols_resolved() {
         return library.workspace_create && library.default_workspace_options && library.get_x && library.is_solved &&
            library.niter && library.elapsed_time && library.warm_start && library.workspace_free;
      }
   } // namespace

   const KrylovLibrary& krylov_library() {
      return library;
   }

   const std::string& krylov_library_error() {
      return error_message;
   }

   const std::string& loaded_krylov_library() {
      return loaded_name;
   }

   void set_default_krylov_library(const std::string& library_name) {
      default_name = library_name.empty() ? UNO_KRYLOV_DEFAULT_LIBRARY : library_name;
   }

   bool load_krylov_library(const std::string& user_library_name) {
      const std::string name = resolve_library_name(user_library_name, "UNO_KRYLOV_LIBRARY", default_name.c_str());
      if (handle != nullptr) {
         if (!user_library_name.empty() && name != loaded_name) {
            WARNING << "Uno: libkrylov '" << loaded_name << "' is already loaded; ignoring the request for '" << name <<
               "' (the first load wins)\n";
         }
         return library.abi != 0;
      }

#ifdef _WIN32
      add_library_directory_to_path(name);
#endif
      handle = open_library(name.c_str());
      if (handle == nullptr) {
         // an explicit request deserves a warning (e.g. a missing file, or an incompatible libstdc++ already loaded)
         report("could not load libkrylov '" + name + "': " + last_library_error(), name != UNO_KRYLOV_DEFAULT_LIBRARY);
         return false;
      }
      loaded_name = name;

      resolve_c(library.workspace_create, "krylov_workspace_create");
      resolve_c(library.default_workspace_options, "krylov_default_workspace_options");
      resolve_c(library.get_x, "krylov_get_x");
      resolve_c(library.is_solved, "krylov_is_solved");
      resolve_c(library.niter, "krylov_niter");
      resolve_c(library.elapsed_time, "krylov_elapsed_time");
      resolve_c(library.warm_start, "krylov_warm_start");
      resolve_c(library.workspace_free, "krylov_workspace_free");
      // the ABI is determined by the presence of krylov_get_version
      resolve_c(library.get_version, "krylov_get_version");
      if (library.get_version != nullptr) {
         resolve_c(library.default_options_abi2, "krylov_default_options");
         resolve_c(library.solve_abi2, "krylov_solve");
         library.abi = (common_symbols_resolved() && library.default_options_abi2 && library.solve_abi2) ? 2 : 0;
      }
      else {
         resolve_c(library.default_options_abi1, "krylov_default_options");
         resolve_c(library.solve_abi1, "krylov_solve");
         library.abi = (common_symbols_resolved() && library.default_options_abi1 && library.solve_abi1) ? 1 : 0;
      }
      if (library.abi == 0) {
         std::string missing;
         const std::pair<const void*, const char*> symbols[] = {
            {reinterpret_cast<const void*>(library.workspace_create), "krylov_workspace_create"},
            {reinterpret_cast<const void*>(library.default_workspace_options), "krylov_default_workspace_options"},
            {reinterpret_cast<const void*>(library.get_x), "krylov_get_x"},
            {reinterpret_cast<const void*>(library.is_solved), "krylov_is_solved"},
            {reinterpret_cast<const void*>(library.niter), "krylov_niter"},
            {reinterpret_cast<const void*>(library.elapsed_time), "krylov_elapsed_time"},
            {reinterpret_cast<const void*>(library.warm_start), "krylov_warm_start"},
            {reinterpret_cast<const void*>(library.workspace_free), "krylov_workspace_free"},
            {(library.get_version != nullptr) ? reinterpret_cast<const void*>(library.default_options_abi2) :
               reinterpret_cast<const void*>(library.default_options_abi1), "krylov_default_options"},
            {(library.get_version != nullptr) ? reinterpret_cast<const void*>(library.solve_abi2) :
               reinterpret_cast<const void*>(library.solve_abi1), "krylov_solve"}
         };
         for (const auto& [pointer, symbol]: symbols) {
            if (pointer == nullptr) {
               missing += std::string(missing.empty() ? "" : ", ") + symbol;
            }
         }
         report("'" + name + "' was loaded but does not export the libkrylov symbols (missing: " + missing + ")", true);
         return false;
      }
      initialize_julia_runtime();
      if (library.get_version != nullptr) {
         int major = 0, minor = 0, patch = 0;
         library.get_version(&major, &minor, &patch);
         DEBUG << "Uno: loaded libkrylov '" << name << "' (Krylov.jl v" << major << '.' << minor << '.' << patch << ")\n";
      }
      else {
         DEBUG << "Uno: loaded libkrylov '" << name << "' (ABI 1)\n";
      }
      return true;
   }

   bool is_libkrylov_functional(const std::string& user_library_name) {
      if (!load_krylov_library(user_library_name)) {
         return false;
      }
      if (functional_check == -1) {
         // solve A x = b with CG, A symmetric positive definite (tridiagonal)
         try {
            constexpr size_t n = 3;
            const auto A = [](const double* x, double* y) {
               y[0] = 4. * x[0] + x[1];
               y[1] = x[0] + 3. * x[1] + x[2];
               y[2] = x[1] + 2. * x[2];
            };
            const double b[n] = {1., 2., 3.};
            double x[n] = {0., 0., 0.};
            KrylovSolver solver(KrylovMethod::CG, n, n);
            KrylovParameters parameters;
            parameters.atol = 1e-12;
            parameters.rtol = 1e-12;
            const KrylovResult result = solver.solve(A, b, x, parameters);
            double y[n];
            A(x, y);
            const double residual = std::fabs(y[0] - b[0]) + std::fabs(y[1] - b[1]) + std::fabs(y[2] - b[2]);
            functional_check = (result.solved && residual <= 1e-8) ? 1 : 0;
            if (functional_check == 0) {
               report("libkrylov is not functional: CG solve " + std::string(result.solved ? "solved" : "not solved") +
                  " in " + std::to_string(result.iterations) + " iterations, residual " + std::to_string(residual), true);
            }
         }
         catch (const std::exception& exception) {
            report(std::string("libkrylov is not functional: ") + exception.what(), true);
            functional_check = 0;
         }
      }
      return functional_check == 1;
   }
} // namespace
