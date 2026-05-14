// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>
#include "UnoSolverWrapper.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "tools/Logger.hpp"

namespace uno {
   UnopyUserCallbacks::UnopyUserCallbacks(std::optional<NotifyAcceptableIterateCallback>& notify_acceptable_iterate_callback):
         UserCallbacks(), notify_acceptable_iterate_callback(notify_acceptable_iterate_callback) { }

   void UnopyUserCallbacks::notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers,
         double objective_multiplier, double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) {
      if (this->notify_acceptable_iterate_callback.has_value()) {
         const auto primals_py = to_const_array(primals.data(), primals.size());
         const auto lower_bound_multipliers_py = to_const_array(multipliers.lower_bounds.data(), multipliers.lower_bounds.size());
         const auto upper_bound_multipliers_py = to_const_array(multipliers.upper_bounds.data(), multipliers.upper_bounds.size());
         const auto constraint_multipliers_py = to_const_array(multipliers.constraints.data(), multipliers.constraints.size());
         (*this->notify_acceptable_iterate_callback)(primals_py, lower_bound_multipliers_py, upper_bound_multipliers_py,
            constraint_multipliers_py, objective_multiplier, primal_feasibility_residual, stationarity_residual,
            complementarity_residual);
      }
   }

   bool UnopyUserCallbacks::user_termination(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/,
         double /*objective_multiplier*/, double /*primal_feasibility_residual*/, double /*stationarity_residual*/,
         double /*complementarity_residual*/) {
      return false;
   }

   UnoSolverWrapper::UnoSolverWrapper() {
      DefaultOptions::load(this->options);
   }

   void UnoSolverWrapper::set_logger_stream(py::object python_stream) {
      this->stream = python_stream;  // keep Python object alive
      this->stream_buffer = std::make_unique<PythonStreamBuffer>(python_stream);
      this->ostream = std::make_unique<std::ostream>(this->stream_buffer.get());
      uno::Logger::set_stream(*this->ostream);
   }

   void UnoSolverWrapper::set_notify_acceptable_iterate_callback(NotifyAcceptableIterateCallback notify_acceptable_iterate_callback) {
      this->notify_acceptable_iterate_callback = std::move(notify_acceptable_iterate_callback);
   }

   Result UnoSolverWrapper::optimize(const PythonUserModel& user_model) {
      const PythonModel model{user_model};
      Logger::set_logger(this->options.get_string("logger"));
      UnopyUserCallbacks callbacks{this->notify_acceptable_iterate_callback};
      return this->uno_solver.solve(model, this->options, callbacks);
   }
} // namespace
