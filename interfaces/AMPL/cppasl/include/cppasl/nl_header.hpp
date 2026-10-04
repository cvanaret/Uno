#pragma once

#include <cstddef>
#include <string>

namespace cppasl {

   /// Encoding of the body of an .nl file (the ten header lines are always text).
   enum class NlFormat {
      Text,               ///< 'g': one token group per line
      Binary,             ///< 'b': 4-byte integers and opcodes
      BinaryShortOpcodes, ///< 'z': 4-byte integers, 2-byte opcodes
      BinaryLongIntegers  ///< 'h': 8-byte integers, 2-byte opcodes
   };

   /// The ten header lines of an .nl file (D. M. Gay, "Writing .nl Files", 2005).
   struct NlHeader {
      NlFormat format{NlFormat::Text};
      std::string problem_name;
      int ampl_options[10]{};
      double variable_bound_tolerance{0.};

      // line 2
      int number_variables{0};
      int number_constraints{0};
      int number_objectives{0};
      int number_ranges{0};
      int number_equations{0};
      int number_logical_constraints{0};
      // line 3
      int number_nonlinear_constraints{0};
      int number_nonlinear_objectives{0};
      int number_linear_complementarities{0};
      int number_nonlinear_complementarities{0};
      int number_double_inequality_complementarities{0};
      int number_nonzero_lower_bound_complementarities{0};
      // line 4
      int number_nonlinear_network_constraints{0};
      int number_linear_network_constraints{0};
      // line 5 (AMPL orders nonlinear variables first)
      int number_nonlinear_variables_in_constraints{0};
      int number_nonlinear_variables_in_objectives{0};
      int number_nonlinear_variables_in_both{0};
      // line 6
      int number_linear_network_variables{0};
      int number_imported_functions{0};
      int arithmetic_kind{0};
      int flags{0};
      // line 7
      int number_binary_variables{0};
      int number_integer_variables{0};
      int number_nonlinear_integer_variables_in_both{0};
      int number_nonlinear_integer_variables_in_constraints{0};
      int number_nonlinear_integer_variables_in_objectives{0};
      // line 8
      std::size_t number_jacobian_nonzeros{0};
      std::size_t number_objective_gradient_nonzeros{0};
      // line 9
      int maximum_constraint_name_length{0};
      int maximum_variable_name_length{0};
      // line 10: common expressions ("defined variables")
      int number_common_expressions_in_both{0};
      int number_common_expressions_in_constraints{0};
      int number_common_expressions_in_objectives{0};
      int number_common_expressions_in_single_constraint{0};
      int number_common_expressions_in_single_objective{0};

      [[nodiscard]] int number_defined_variables() const {
         return this->number_common_expressions_in_both + this->number_common_expressions_in_constraints +
            this->number_common_expressions_in_objectives + this->number_common_expressions_in_single_constraint +
            this->number_common_expressions_in_single_objective;
      }
   };

} // namespace cppasl
