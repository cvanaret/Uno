// Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MATLABUSERCALLBACK_H
#define UNO_MATLABUSERCALLBACK_H

#include "tools/UserCallbacks.hpp"
#include "../unomex/unomex_function.hpp"

namespace uno {

    // Matlab user callbacks
    class MatlabUserCallbacks : public UserCallbacks {
    public:
        MatlabUserCallbacks(handle_t notify_acceptable_iterate_callback, handle_t user_termination_callback);

        void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier, double primal_feasibility, double stationarity, double complementarity) override;

        bool user_termination(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier, double primal_feasibility, double stationarity, double complementarity) override;

    private:
        handle_t notify_acceptable_iterate_callback;
        handle_t user_termination_callback;
    };

}; // namespace

#endif // UNO_MATLABUSERCALLBACK_H