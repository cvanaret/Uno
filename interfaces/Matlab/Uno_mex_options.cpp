// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "mex.h"
#include "Uno_mex_utilities.hpp"
#include "Uno.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Presets.hpp"

using namespace uno;

// gateway function
// options = uno_options([preset]);
void mexFunction( int /* nlhs */, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate arguments
    if (nrhs > 1) {
        mexErrMsgIdAndTxt("uno:error", "Too many input arguments.");
    }
    // preset (optional)
    std::string preset;
    if (nrhs == 1 && !mxIsEmpty(prhs[0])) {
        if (!mxIsChar(prhs[0]) && !mxIsClass(prhs[0],"string")) {
            mexErrMsgIdAndTxt("uno:error", "Invalid argument at position 1. Value must be of type char or string.");
        }
        preset = mxArray_to_string(prhs[0]);
    }

    // create Uno options
    Options uno_options;

    // default options
    DefaultOptions::load(uno_options);

    // add default preset
    const Options preset_options = Presets::get_preset_options(std::nullopt);
    uno_options.overwrite_with(preset_options);

    // set preset
    if (!preset.empty()) {
        try {
            Presets::set(uno_options, preset);
        }
        catch (const std::runtime_error& err) {
            mexErrMsgIdAndTxt("uno:error", "Invalid argument at position 1. %s", err.what());
        }
    }

    // options
    MxStruct options = options_to_mxStruct(uno_options);

    // output
    plhs[0] = mxStruct_to_mxArray(options);
}