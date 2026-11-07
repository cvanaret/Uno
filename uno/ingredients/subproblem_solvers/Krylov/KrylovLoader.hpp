// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_KRYLOVLOADER_H
#define UNO_KRYLOVLOADER_H

// Runtime loading of libkrylov (the C interface of Krylov.jl, compiled with juliac): libkrylov is never linked, it is
// dlopen'd on first use and its symbols are resolved by name. krylov.h is not included: the C ABI is declared below.
// Two ABIs are supported, distinguished by the presence of krylov_get_version:
// - ABI 1 (Krylov.jl artifacts up to v0.0.33): krylov_solve/krylov_block_solve without right preconditioner, KrylovOptions
//   with 7 fields;
// - ABI 2 (krylov_get_version exported): right preconditioner matvec_N, KrylovOptions with 12 fields.

#include <string>

namespace uno {
   extern "C" {
      // callbacks
      using KrylovMatvec = void (*)(const void* x, void* y, void* userdata);
      using KrylovBlockMatvec = void (*)(const void* X, void* Y, int p, void* userdata);

      struct KrylovWorkspaceOptions {
         int memory;
         int window;
      };

      struct KrylovOptionsABI1 {
         double atol;
         double rtol;
         int itmax;
         int verbose;
         double lambda;
         double tau;
         double nu;
      };

      struct KrylovOptionsABI2 {
         double atol;
         double rtol;
         int itmax;
         int verbose;
         double lambda;
         double tau;
         double nu;
         double timemax;
         double radius;
         int restart;
         int reorthogonalization;
         int linesearch;
      };

      using krylov_workspace_create_fp = int (*)(int solver, int m, int n, int dtype, int device,
         const KrylovWorkspaceOptions* wopts, void** ws_out);
      using krylov_default_workspace_options_fp = KrylovWorkspaceOptions (*)();
      using krylov_default_options_abi1_fp = KrylovOptionsABI1 (*)();
      using krylov_default_options_abi2_fp = KrylovOptionsABI2 (*)();
      using krylov_get_version_fp = void (*)(int* major, int* minor, int* patch);
      using krylov_solve_abi1_fp = int (*)(void* ws, KrylovMatvec matvec_A, KrylovMatvec matvec_At, KrylovMatvec matvec_M,
         const void* b, const void* c, void* userdata, const KrylovOptionsABI1* opts);
      using krylov_solve_abi2_fp = int (*)(void* ws, KrylovMatvec matvec_A, KrylovMatvec matvec_At, KrylovMatvec matvec_M,
         KrylovMatvec matvec_N, const void* b, const void* c, void* userdata, const KrylovOptionsABI2* opts);
      using krylov_get_vector_fp = int (*)(void* ws, void* x, int n);
      using krylov_workspace_query_fp = int (*)(void* ws);
      using krylov_elapsed_time_fp = double (*)(void* ws);
      using krylov_warm_start_fp = int (*)(void* ws, const void* x0, int n);
   }

   // libkrylov enumerators (KrylovSolverType, KrylovDataType, KrylovDeviceType)
   enum class KrylovMethod: int {
      CG = 0, CR = 1, SYMMLQ = 2, MINRES = 3, MINRES_QLP = 4, DIOM = 5, DQGMRES = 6, FOM = 7, GMRES = 8, FGMRES = 9,
      BICGSTAB = 10, CGS = 11, BILQ = 12, QMR = 13, USYMLQ = 14, USYMQR = 15, TRICG = 16, TRIMR = 17, TRILQR = 18,
      BILQR = 19, LSLQ = 20, LSQR = 21, LSMR = 22, USYMLQR = 23, CGLS = 24, CRLS = 25, CGNE = 26, CRMR = 27, CRAIG = 28,
      CRAIGMR = 29, LNLQ = 30, GPMR = 31, CAR = 32, MINARES = 33
   };
   constexpr int KRYLOV_FLOAT64 = 1;
   constexpr int KRYLOV_CPU = 0;

   struct KrylovLibrary {
      int abi{0}; // 0: not loaded, 1 or 2: see above
      krylov_workspace_create_fp workspace_create{nullptr};
      krylov_default_workspace_options_fp default_workspace_options{nullptr};
      krylov_default_options_abi1_fp default_options_abi1{nullptr};
      krylov_default_options_abi2_fp default_options_abi2{nullptr};
      krylov_get_version_fp get_version{nullptr}; // ABI 2 only
      krylov_solve_abi1_fp solve_abi1{nullptr};
      krylov_solve_abi2_fp solve_abi2{nullptr};
      krylov_get_vector_fp get_x{nullptr};
      krylov_workspace_query_fp is_solved{nullptr};
      krylov_workspace_query_fp niter{nullptr};
      krylov_elapsed_time_fp elapsed_time{nullptr};
      krylov_warm_start_fp warm_start{nullptr};
      krylov_workspace_query_fp workspace_free{nullptr};
   };

   // the symbols resolved by load_krylov_library(); abi == 0 until loaded
   [[nodiscard]] const KrylovLibrary& krylov_library();

   // dlopen libkrylov and resolve its symbols. The library name is resolved in this order: user_library_name (option
   // libkrylov_path) > UNO_KRYLOV_LIBRARY environment variable > default set by set_default_krylov_library (e.g. the
   // library bundled with unopy) > platform default libkrylov.{so,dylib,dll}. Idempotent: the first successful load wins.
   // Returns true if all the required symbols were resolved
   bool load_krylov_library(const std::string& user_library_name = "");
   // why libkrylov could not be loaded or is not functional (empty if no failure)
   [[nodiscard]] const std::string& krylov_library_error();
   // name of the loaded library (empty if not loaded)
   [[nodiscard]] const std::string& loaded_krylov_library();
   // fallback library name, used when neither libkrylov_path nor UNO_KRYLOV_LIBRARY is set
   void set_default_krylov_library(const std::string& library_name);
   // loads libkrylov and checks that it is functional by solving a small linear system with CG. Cached per library
   [[nodiscard]] bool is_libkrylov_functional(const std::string& user_library_name = "");
} // namespace

#endif // UNO_KRYLOVLOADER_H
