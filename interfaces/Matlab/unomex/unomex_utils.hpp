// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNOMEX_UTILS_H
#define UNOMEX_UTILS_H

#include <stdint.h>
#include <vector>
#include <algorithm>
#include "linear_algebra/Vector.hpp"
#include "mex.h"

namespace uno {

    // guard on std::vector<mxArray*> for automatic memory management
    struct MxArrayVectorGuard {
        const std::vector<mxArray*>& vector;
        explicit MxArrayVectorGuard(const std::vector<mxArray*>& vector);
        ~MxArrayVectorGuard();
    };

    constexpr int mxSTRING_CLASS = 19; // cf. https://it.mathworks.com/matlabcentral/answers/2102376-how-do-i-process-a-string-class-in-a-mex-function

    template<typename T>
    mxClassID get_mxClassID();

    bool isvalid(const mxArray* arr);
    bool isempty(const mxArray* arr);
    
    template<typename T>
    bool isa(const mxArray* arr) {
        return mxGetClassID(arr) == get_mxClassID<T>();
    }
    
    bool has_size(const mxArray* arr, const size_t nrows, const size_t ncolumns);
    bool ispositive(const mxArray* arr);
    bool isinteger(const mxArray* arr);
    bool isunitary(const mxArray* arr);
    bool isscalar(const mxArray* arr);
    bool isvector(const mxArray* arr, const size_t len);

    template <typename OutType, typename InType>
    Vector<OutType> convert_vector_type(const Vector<InType>& input) {
        Vector<OutType> output(input.size());
        std::transform(input.begin(), input.end(), output.begin(),
                    [](const InType& val) { return static_cast<OutType>(val); });
        return output;
    }

}; // namespace

#endif // UNOMEX_UTILS_H