// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// Second of two C translation units that include the C API header (see HeaderLinkageA.c).

#include "../../interfaces/C/Uno_C_API.h"

uno_int uno_header_linkage_b(void) {
   return UNO_MAXIMIZE + UNO_ZERO_BASED_INDEXING + UNO_VERSION_MAJOR + (uno_int)(UNO_LOWER_TRIANGLE == 'L');
}
