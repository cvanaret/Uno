#include "ProblemFactory.hpp"
#include "interfaces/AMPL/AMPLModel.hpp"
#include "SlackReformulation.hpp"

std::unique_ptr<Problem> ProblemFactory::create(const std::string& problem_name) {
   return std::make_unique<AMPLModel>(problem_name);
}