// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <gtest/gtest.h>
#include "../../interfaces/C/Uno_C_API.h"

// defined in unotest/C_header/HeaderLinkage{A,B}.c, two C translation units that include the C API header
extern "C" {
   uno_int uno_header_linkage_a(void);
   uno_int uno_header_linkage_b(void);
}

// the executable links at all only if the header constants have internal linkage in C
TEST(CAPIHeader, ConstantsCanBeIncludedFromSeveralCTranslationUnits) {
   EXPECT_EQ(uno_header_linkage_a(), UNO_MINIMIZE + UNO_ONE_BASED_INDEXING + UNO_VERSION_MAJOR + 1);
   EXPECT_EQ(uno_header_linkage_b(), UNO_MAXIMIZE + UNO_ZERO_BASED_INDEXING + UNO_VERSION_MAJOR + 1);
}
