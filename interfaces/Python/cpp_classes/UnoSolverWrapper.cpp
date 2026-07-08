// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>
#include <string>
#include "UnoSolverWrapper.hpp"
#include "model/Model.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Presets.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   UnopyUserCallbacks::UnopyUserCallbacks(std::optional<NotifyAcceptableIterateCallback>& notify_acceptable_iterate_callback,
         std::optional<TerminationCallback>& termination_callback):
      UserCallbacks(), notify_acceptable_iterate_callback(notify_acceptable_iterate_callback), termination_callback(termination_callback) { }

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

   bool UnopyUserCallbacks::termination(const Vector<double>& primals, const Multipliers& multipliers,
         double objective_multiplier, double primal_feasibility_residual, double stationarity_residual,
         double complementarity_residual) {
      if (this->termination_callback.has_value()) {
         const auto primals_py = to_const_array(primals.data(), primals.size());
         const auto lower_bound_multipliers_py = to_const_array(multipliers.lower_bounds.data(), multipliers.lower_bounds.size());
         const auto upper_bound_multipliers_py = to_const_array(multipliers.upper_bounds.data(), multipliers.upper_bounds.size());
         const auto constraint_multipliers_py = to_const_array(multipliers.constraints.data(), multipliers.constraints.size());
         return (*this->termination_callback)(primals_py, lower_bound_multipliers_py, upper_bound_multipliers_py,
            constraint_multipliers_py, objective_multiplier, primal_feasibility_residual, stationarity_residual,
            complementarity_residual);
      }
      return false;
   }

   UnoSolverWrapper::UnoSolverWrapper() {
      DefaultOptions::load(this->user_options);
   }

   UnoSolverWrapper::~UnoSolverWrapper() {
      // Logger holds a non-owning reference to *this->ostream. Repoint it at a process-lifetime sink now: the
      // destructor body runs *before* stream_buffer/ostream are destroyed. Otherwise the next Logger write
      // dereferences a freed ostream -> SIGSEGV in sentry
      Logger::set_stream(std::cout);
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

   void UnoSolverWrapper::set_termination_callback(TerminationCallback termination_callback) {
      this->termination_callback = std::move(termination_callback);
   }

   Result UnoSolverWrapper::optimize(const PythonUserModel& user_model) {
      // create an instance of PythonModel, a subclass of Model
      const PythonModel model{user_model};
      Logger::set_logger(this->user_options.get_string("logger"));

      // set the preset (default: auto)
      Options full_options;
      Presets::set(model, full_options, this->user_options.get_string("preset"));

      // copy the rest of the options
      full_options.overwrite(this->user_options);

      // solve the model
      UnopyUserCallbacks callbacks{this->notify_acceptable_iterate_callback, this->termination_callback};
      return this->uno_solver.solve(model, full_options, callbacks);
   }

   const std::string& UnoSolverWrapper::get_method_description() const {
      return this->uno_solver.get_method_description();
   }
} // namespace
