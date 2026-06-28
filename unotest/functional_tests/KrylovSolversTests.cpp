// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/Krylov/KrylovSolvers.hpp"

/*
 * basic_cg.c — minimal example: solve a 5x5 SPD system with CG.
 *
 * A = tridiag(-1, 2, -1),  b = [1, 0, 0, 0, 1]^T
 *
 * Compile (after building libkrylov with CMake):
 *
 *   gcc -o basic_cg basic_cg.c -I../include -L../build -lkrylov -Wl,-rpath,../build
 *
 * Expected output:
 *   Solved: yes   niter: 3   time: ...
 *   x = [ 1.00  1.00  1.00  1.00  1.00 ]
 */

/* -------------------------------------------------------------------------
 * Problem data
 * ------------------------------------------------------------------------- */

#define N 5

/* Tridiagonal matrix stored as three diagonals for simplicity. */
typedef struct {
  int    n;
  double diag[N];   /* main diagonal */
  double off[N-1];  /* sub/super diagonal */
} TriDiag;

/* Matvec callback:  y = A * x  */
static void matvec_A(const void *xv, void *yv, void *userdata) {
  const double *x  = (const double *)xv;
  double       *y  = (double *)yv;
  const TriDiag *A = (const TriDiag *)userdata;
  int n = A->n;

  for (int i = 0; i < n; i++) {
    y[i] = A->diag[i] * x[i];
    if (i > 0)   y[i] += A->off[i-1] * x[i-1];
    if (i < n-1) y[i] += A->off[i]   * x[i+1];
  }
}

TEST(Krylov, Example) {
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
  ASSERT_EQ(ret, 0);

  /* -----------------------------------------------------------------------
   * Solve
   * --------------------------------------------------------------------- */
  KrylovOptions opts = krylov_default_options();
  opts.atol = 1e-10;   /* absolute tolerance */
  opts.rtol = 1e-10;   /* relative tolerance */
  ret = krylov_solve(ws,
                     matvec_A,   /* y = A*x */
                     NULL,       /* y = Aᴴ*x  (CG doesn't need it) */
                     NULL,       /* matvec_M: left/centered preconditioner (none) */
                     NULL,       /* matvec_N: right preconditioner (none) */
                     b,          /* right-hand side b (size m) */
                     NULL,       /* c = NULL  (CG only needs one RHS) */
                     &A,         /* userdata forwarded to matvec_A */
                     &opts);     /* solver options (NULL = all defaults) */
  ASSERT_EQ(ret, 0);

  /* -----------------------------------------------------------------------
   * Retrieve results
   * --------------------------------------------------------------------- */
  ret = krylov_get_x(ws, x, N);
  ASSERT_EQ(ret, 0);

  const double tolerance = 1e-8;
  for (int i = 0; i < N; i++) {
    EXPECT_NEAR(x[i], 1., tolerance);
  }

  /* -----------------------------------------------------------------------
   * Free workspace
   * --------------------------------------------------------------------- */
  krylov_workspace_free(ws);
}