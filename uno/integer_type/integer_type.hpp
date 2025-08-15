// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTEGERTYPE_H
#define UNO_INTEGERTYPE_H

#include <cstdint>

#ifdef USE_64BIT_INTEGERS
using UnoInt = int64_t;
#else
using UnoInt = int32_t;
#endif

#endif // UNO_INTEGERTYPE_H