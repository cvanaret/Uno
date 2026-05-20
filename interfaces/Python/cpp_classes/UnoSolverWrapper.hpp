// Copyright (c) 2025-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOSOLVERWRAPPER_H
#define UNO_UNOSOLVERWRAPPER_H

#include <optional>
#include "Uno.hpp"
#include "options/Options.hpp"
#include "tools/UserCallbacks.hpp"
#include "PythonModel.hpp"
#include "PythonStreamBuffer.hpp"
#include "../unopy.hpp"

namespace uno {
   class UnopyUserCallbacks: public UserCallbacks {
   public:
      UnopyUserCallbacks(std::optional<NotifyAcceptableIterateCallback>& notify_acceptable_iterate_callback,
         std::optional<TerminationCallback>& termination_callback);

      void notify_acceptable_iterate(const Vector<double>& primals, const Multipliers& multipliers, double objective_multiplier,
            double primal_feasibility_residual, double stationarity_residual, double complementarity_residual) override;

      bool termination(const Vector<double>& /*primals*/, const Multipliers& /*multipliers*/, double /*objective_multiplier*/,
            double /*primal_feasibility_residual*/, double /*stationarity_residual*/, double /*complementarity_residual*/) override;

   protected:
      std::optional<NotifyAcceptableIterateCallback>& notify_acceptable_iterate_callback;
      std::optional<TerminationCallback>& termination_callback;
   };

   class UnoSolverWrapper {
   public:
      Uno uno_solver{};
      Options options{};

      UnoSolverWrapper();

      void set_logger_stream(py::object py_stream);
      void set_notify_acceptable_iterate_callback(NotifyAcceptableIterateCallback notify_acceptable_iterate_callback);
      void set_termination_callback(TerminationCallback termination_callback);
      [[nodiscard]] Result optimize(const PythonUserModel& user_model);
      [[nodiscard]] const std::string& get_method_description() const;

   private:
      py::object stream;
      std::unique_ptr<PythonStreamBuffer> stream_buffer;
      std::unique_ptr<std::ostream> ostream;
      std::optional<NotifyAcceptableIterateCallback> notify_acceptable_iterate_callback{};
      std::optional<TerminationCallback> termination_callback{};
   };
} // namespace

#endif // UNO_UNOSOLVERWRAPPER_H