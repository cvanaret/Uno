// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// First of two C translation units that include the C API header. Linking both into the same executable checks that
// the header does not define symbols with external linkage (in C, a file-scope "const" object has external linkage).

#include "../../interfaces/C/Uno_C_API.h"

uno_int uno_header_linkage_a(void) {
   return UNO_MINIMIZE + UNO_ONE_BASED_INDEXING + UNO_VERSION_MAJOR + (uno_int)(UNO_PROBLEM_NONLINEAR[0] == 'N');
}
