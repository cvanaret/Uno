// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_KRYLOVSOLVERS_H
#define UNO_KRYLOVSOLVERS_H

#include <krylov.h>

namespace uno {
#ifdef HAS_KRYLOV
   void test_krylov_solvers();
#endif
} // namespace

#endif // UNO_KRYLOVSOLVERS_H