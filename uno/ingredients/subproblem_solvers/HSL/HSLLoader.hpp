// Copyright (c) 2026 Joris Gillis, Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HSLLOADER_H
#define UNO_HSLLOADER_H

// IPOPT-style runtime loading of the HSL routines: instead of linking libhsl at
// build time, dlopen it on first use and resolve MA27/MA57 (Fortran) and MA86/MC68
// (C interface) by symbol name. This keeps Uno build- and ship-able without HSL;
// the user supplies libhsl at runtime.

#include <string>

namespace uno {
   extern "C" {
      struct ma86_control;
      struct ma86_info;
      struct mc68_control;
      struct mc68_info;

      // MA57
      using ma57id_fp = void (*)(double cntl[], int icntl[]);
      using ma57ad_fp = void (*)(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep,
         int keep[], int iwork[], int icntl[], int info[], double rinfo[]);
      using ma57bd_fp = void (*)(const int* n, int* ne, const double a[], double fact[], const int* lfact,
         int ifact[], const int* lifact, const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[],
         int info[], double rinfo[]);
      using ma57cd_fp = void (*)(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact,
         const int* nrhs, double rhs[], const int* lrhs, double work[], int* lwork, int iwork[], int icntl[], int info[]);
      using ma57dd_fp = void (*)(const int* job, const int* n, int* ne, const double a[], const int irn[],
         const int jcn[], double fact[], int* lfact, int ifact[], int* lifact, const double rhs[], double x[], double resid[],
         double work[], int iwork[], int icntl[], double cntl[], int info[], double rinfo[]);
      using ma57ed_fp = void (*)(const int* n, const int* ic, int keep[], const double fact[], const int* lfact,
         double newfac[], const int* lnew, const int ifact[], const int* lifact, int newifc[], const int* linew, int info[]);

      // MA27
      using ma27id_fp = void (*)(int ICNTL[], double CNTL[]);
      using ma27ad_fp = void (*)(int* N, int* NZ, int IRN[], int ICN[], int IW[], int* LIW, int IKEEP[], int IW1[],
         int* NSTEPS, int* IFLAG, int ICNTL[], double CNTL[], int INFO[], double* OPS);
      using ma27bd_fp = void (*)(int* N, int* NZ, int IRN[], int ICN[], double A[], int* LA, int IW[], int* LIW,
         int IKEEP[], int* NSTEPS, int* MAXFRT, int IW1[], int ICNTL[], double CNTL[], int INFO[]);
      using ma27cd_fp = void (*)(int* N, double A[], int* LA, int IW[], int* LIW, double W[], int* MAXFRT, double RHS[],
         int IW1[], int* NSTEPS, int ICNTL[], int INFO[]);

      // MA86
      using ma86_default_control_fp = void (*)(struct ma86_control* control);
      using ma86_analyse_fp = void (*)(const int n, const int ptr[], const int row[], int order[], void** keep,
         const struct ma86_control* control, struct ma86_info* info);
      using ma86_factor_fp = void (*)(const int n, const int ptr[], const int row[], const double val[], const int order[],
         void** keep, const struct ma86_control* control, struct ma86_info* info, const double scale[]);
      using ma86_solve_fp = void (*)(const int job, const int nrhs, const int ldx, double* x, const int order[], void** keep,
         const struct ma86_control* control, struct ma86_info* info, const double scale[]);
      using ma86_finalise_fp = void (*)(void** keep, const struct ma86_control* control);

      // MC68
      using mc68_default_control_fp = void (*)(struct mc68_control* control);
      using mc68_order_fp = void (*)(const int ord, const int n, const int ptr[], const int row[], int perm[],
         const struct mc68_control* control, struct mc68_info* info);
   }

   // Resolved by load_hsl_library(); nullptr until loaded or if the symbol is absent.
   extern ma57id_fp hsl_ma57id;
   extern ma57ad_fp hsl_ma57ad;
   extern ma57bd_fp hsl_ma57bd;
   extern ma57cd_fp hsl_ma57cd;
   extern ma57dd_fp hsl_ma57dd;
   extern ma57ed_fp hsl_ma57ed;
   extern ma27id_fp hsl_ma27id;
   extern ma27ad_fp hsl_ma27ad;
   extern ma27bd_fp hsl_ma27bd;
   extern ma27cd_fp hsl_ma27cd;
   extern ma86_default_control_fp hsl_ma86_default_control;
   extern ma86_analyse_fp hsl_ma86_analyse;
   extern ma86_factor_fp hsl_ma86_factor;
   extern ma86_solve_fp hsl_ma86_solve;
   extern ma86_finalise_fp hsl_ma86_finalise;
   extern mc68_default_control_fp hsl_mc68_default_control;
   extern mc68_order_fp hsl_mc68_order;

   // dlopen the HSL library (default name per platform, overridable via the
   // UNO_HSL_LIBRARY environment variable) and resolve the symbols. Idempotent.
   bool load_hsl_library(const std::string& user_library_name = "");
   bool ma57_symbols_available();
   bool ma27_symbols_available();
   bool ma86_symbols_available();
} // namespace

#endif // UNO_HSLLOADER_H
