#ifndef UNO_MULTIPLIERS_H
#define UNO_MULTIPLIERS_H

struct Multipliers {
   std::vector<double> lower_bounds; /*!< Multipliers of the lower bound constraints */
   std::vector<double> upper_bounds; /*!< Multipliers of the lower bound constraints */
   std::vector<double> constraints; /*!< Multipliers of the general constraints */
   double objective{1.};

   Multipliers(size_t number_variables, size_t number_constraints);
};

inline Multipliers::Multipliers(size_t number_variables, size_t number_constraints) : lower_bounds(number_variables),
      upper_bounds(number_variables), constraints(number_constraints) {
}

#endif // UNO_MULTIPLIERS_H
