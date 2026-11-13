/* Copyright (c) 2026 Charlie Vanaret
   Licensed under the MIT license. See LICENSE file in the project directory for details. */

/* Stub HSL library for testing HSL_RUNTIME_LOADING without the real (licensed) libhsl.
   It exports only LIBHSL_version and the MA27 routines (as a libhsl built from the HSL Archive would),
   with gfortran name mangling. MA57, MA86 and MC68 are deliberately absent.
   Only ma27id is functional (it counts its calls, so that tests can check that Uno called into this library);
   the other routines abort, since a test that reaches them would be solving with a fake solver. */

#include <stdlib.h>

#if defined(_WIN32)
#define UNO_STUB_EXPORT __declspec(dllexport)
#else
#define UNO_STUB_EXPORT __attribute__((visibility("default")))
#endif

UNO_STUB_EXPORT int uno_hsl_stub_ma27id_calls = 0;

UNO_STUB_EXPORT void LIBHSL_version(int* major, int* minor, int* patch) {
   *major = 0;
   *minor = 0;
   *patch = 0;
}

// set default parameters
UNO_STUB_EXPORT void ma27id_(int icntl[], double cntl[]) {
   (void)icntl;
   (void)cntl;
   uno_hsl_stub_ma27id_calls++;
}

UNO_STUB_EXPORT void ma27ad_(void) {
   abort();
}

UNO_STUB_EXPORT void ma27bd_(void) {
   abort();
}

UNO_STUB_EXPORT void ma27cd_(void) {
   abort();
}
