#include "LocalApproximationFactory.hpp"
#include "QPApproximation.hpp"
#include "BQPDSolver.hpp"

std::shared_ptr<LocalApproximation> LocalApproximationFactory::create(const std::string& type, QPSolver& solver) {
	if (type == "QP") {
		return std::make_shared<QPApproximation>(solver);
	}
	else {
		throw std::invalid_argument("LocalApproximation type " + type + " does not exist");
	}
}
