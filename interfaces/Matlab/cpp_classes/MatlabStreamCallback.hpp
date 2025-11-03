// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATLABSTREAMCALLBACK_H
#define UNO_MATLABSTREAMCALLBACK_H

#include "../../UserModel.hpp"
#include "../unomex/unomex_function.hpp"

namespace uno {

    class MatlabStreamCallback : public UserStreamCallback {
    public: 
        MatlabStreamCallback(handle_t logger_stream_callback);
        ~MatlabStreamCallback() override = default;

    int32_t operator()(const char* buf, int32_t len) const override;

    private:
        handle_t logger_stream_callback;
    };

}; // namespace

#endif // UNO_MATLABSTREAMCALLBACK_H