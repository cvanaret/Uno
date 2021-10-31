#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include <vector>
#include "Problem.hpp"
#include "Iterate.hpp"

class Preprocessing {
public:
   static void enforce_linear_constraints(const Problem& problem, Iterate& first_iterate);
};

#endif //PREPROCESSING_H
