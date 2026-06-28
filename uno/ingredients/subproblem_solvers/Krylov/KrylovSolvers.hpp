// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_KRYLOVSOLVERS_H
#define UNO_KRYLOVSOLVERS_H

#include <krylov.h>
#include <cassert>
#include <cmath>
#include <iostream>

static void test_krylov_solvers() {
   /*
   * basic_cg.c — minimal example: solve a 5x5 SPD system with CG.
   * A = tridiag(-1, 2, -1),  b = [1, 0, 0, 0, 1]^T
   * Expected output:
   *   Solved: yes   niter: 3   time: ...
   *   x = [ 1.00  1.00  1.00  1.00  1.00 ]
   */

   #define N 5

   /* Tridiagonal matrix stored as three diagonals for simplicity. */
   typedef struct {
      int    n;
      double diag[N];   /* main diagonal */
      double off[N-1];  /* sub/super diagonal */
   } TriDiag;

   /* Matvec callback:  y = A * x  */
   const auto matvec_A = [](const void *xv, void *yv, void *userdata) {
      const double *x  = (const double *) xv;
      double       *y  = (double *) yv;
      const TriDiag *A = (const TriDiag *) userdata;
      const int n = A->n;

      for (int i = 0; i < n; i++) {
         y[i] = A->diag[i] * x[i];
         if (i > 0)   y[i] += A->off[i-1] * x[i-1];
         if (i < n-1) y[i] += A->off[i]   * x[i+1];
      }
   };

   /* Build A = tridiag(-1, 2, -1) */
   TriDiag A;
   A.n = N;
   for (int i = 0; i < N;   i++) A.diag[i] = 2.0;
   for (int i = 0; i < N-1; i++) A.off[i]  = -1.0;

   /* Right-hand side */
   double b[N] = {1.0, 0.0, 0.0, 0.0, 1.0};

   /* Solution buffer */
   double x[N];

   /* -----------------------------------------------------------------------
   * Create workspace for CG, double precision, CPU
   * --------------------------------------------------------------------- */
   void *ws = NULL;
   int ret = krylov_workspace_create(KRYLOV_CG, N, N, KRYLOV_FLOAT64, KRYLOV_CPU,
                                     NULL,   /* workspace options (NULL = defaults) */
                                     &ws);
   assert(ret == 0);

   /* -----------------------------------------------------------------------
   * Solve
   * --------------------------------------------------------------------- */
   KrylovOptions opts = krylov_default_options();
   opts.atol = 1e-10;   /* absolute tolerance */
   opts.rtol = 1e-10;   /* relative tolerance */
   ret = krylov_solve(ws,
                     matvec_A,   /* y = A*x */
                     NULL,       /* y = A'*x  (CG doesn't need it) */
                     NULL,       /* no preconditioner */
                     b,          /* right-hand side b (size m) */
                     NULL,       /* c = NULL  (CG only needs one RHS) */
                     &A,         /* userdata forwarded to matvec_A */
                     &opts);     /* solver options (NULL = all defaults) */
   assert(ret == 0);

   /* -----------------------------------------------------------------------
   * Retrieve results
   * --------------------------------------------------------------------- */
   ret = krylov_get_x(ws, x, N);
   assert(ret == 0);

   std::cout << "KrylovSolvers: CG solution (should be a vector of 1):";
   for (int i = 0; i < N; i++) {
      std::cout << ' ' << x[i];
   }
   std::cout << '\n';
   krylov_workspace_free(ws);
}

#endif // UNO_KRYLOVSOLVERS_H