#ifndef UNO_PROBLEMFACTORY_H
#define UNO_PROBLEMFACTORY_H

#include <memory>
#include "Problem.hpp"
#include "tools/Options.hpp"

class ProblemFactory {
public:
   static std::unique_ptr<Problem> create(const std::string& problem_name);
};

#endif // UNO_PROBLEMFACTORY_H
