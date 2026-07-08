// Copyright (c) 2026 Joris Gillis, Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HSLLoader.hpp"
#include <cctype>
#include <cstdlib>
#include <string>
#include "tools/DynamicLoader.hpp"
#include "tools/Logger.hpp"

// default library name: libhsl.<platform shared-lib extension> (matches IPOPT's hsllib)
#if defined(_WIN32)
#define UNO_HSL_DEFAULT_LIBRARY "libhsl.dll"
#elif defined(__APPLE__)
#define UNO_HSL_DEFAULT_LIBRARY "libhsl.dylib"
#else
#define UNO_HSL_DEFAULT_LIBRARY "libhsl.so"
#endif

namespace uno {
   ma57id_fp hsl_ma57id = nullptr;
   ma57ad_fp hsl_ma57ad = nullptr;
   ma57bd_fp hsl_ma57bd = nullptr;
   ma57cd_fp hsl_ma57cd = nullptr;
   ma57dd_fp hsl_ma57dd = nullptr;
   ma57ed_fp hsl_ma57ed = nullptr;
   ma27id_fp hsl_ma27id = nullptr;
   ma27ad_fp hsl_ma27ad = nullptr;
   ma27bd_fp hsl_ma27bd = nullptr;
   ma27cd_fp hsl_ma27cd = nullptr;

   namespace {
      LibraryHandle hsl_handle = nullptr;
      std::string hsl_loaded_name{}; // the library name that was actually dlopen'd (for the mismatch warning)
   } // anonymous namespace

   bool load_hsl_library(const std::string& library_name) {
      // resolve the effective library name (explicit request > UNO_HSL_LIBRARY env > platform default),
      // done unconditionally so we can compare it against an already-loaded library below.
      std::string name = library_name;
      if (name.empty()) {
         const char* env = std::getenv("UNO_HSL_LIBRARY");
         if (env != nullptr) {
            name = env;
         }
      }
      if (name.empty()) {
         name = UNO_HSL_DEFAULT_LIBRARY;
      }

      // cache success only: a failed probe (e.g. the early available_solvers() check
      // with no name) must not block a later load with an explicit hsllib path.
      if (hsl_handle != nullptr) {
         // the first successful load wins; warn if an explicit, different library was requested too late
         if (!library_name.empty() && name != hsl_loaded_name) {
            WARNING << "Uno: the HSL library '" << hsl_loaded_name << "' is already loaded; ignoring the request for '"
               << name << "' (the first load wins)\n";
         }
         return true;
      }

      hsl_handle = open_library(name.c_str());
      if (hsl_handle == nullptr) {
         DEBUG << "Uno: could not load the HSL library '" << name << "' at runtime\n";
         return false;
      }
      hsl_loaded_name = name;
      DEBUG << "Uno: loaded the HSL library '" << name << "' at runtime\n";

      resolve(hsl_handle, hsl_ma57id, "ma57id");
      resolve(hsl_handle, hsl_ma57ad, "ma57ad");
      resolve(hsl_handle, hsl_ma57bd, "ma57bd");
      resolve(hsl_handle, hsl_ma57cd, "ma57cd");
      resolve(hsl_handle, hsl_ma57dd, "ma57dd");
      resolve(hsl_handle, hsl_ma57ed, "ma57ed");
      resolve(hsl_handle, hsl_ma27id, "ma27id");
      resolve(hsl_handle, hsl_ma27ad, "ma27ad");
      resolve(hsl_handle, hsl_ma27bd, "ma27bd");
      resolve(hsl_handle, hsl_ma27cd, "ma27cd");
      return true;
   }

   bool ma57_symbols_available() {
      load_hsl_library();
      return hsl_ma57id && hsl_ma57ad && hsl_ma57bd && hsl_ma57cd && hsl_ma57dd && hsl_ma57ed;
   }

   bool ma27_symbols_available() {
      load_hsl_library();
      return hsl_ma27id && hsl_ma27ad && hsl_ma27bd && hsl_ma27cd;
   }
} // namespace

// Mirrors the symbol the linked coinhsl meta-library exports; the
// SymmetricIndefiniteLinearSolverFactory uses it to gate MA27/MA57.
extern "C" bool LIBHSL_isfunctional() {
   return uno::ma57_symbols_available() || uno::ma27_symbols_available();
}