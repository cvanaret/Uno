#include "ModelFactory.hpp"
#include "interfaces/AMPL/AMPLModel.hpp"
#include "EqualityConstrainedModel.hpp"

std::unique_ptr<Model> ModelFactory::create(const std::string& problem_name) {
   return std::make_unique<AMPLModel>(problem_name);
}