#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include <vector>
#include "Problem.hpp"
#include "Constraint.hpp"

class Preprocessing {
public:
   static void apply(const Problem& problem, std::vector<double>& x, Multipliers& multipliers);
};

#endif //PREPROCESSING_H
