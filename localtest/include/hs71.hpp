#ifndef NLPTEST_HS71_HPP
#define NLPTEST_HS71_HPP

#include <vector>
#include <limits>
#include "model/Model.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/RectangularMatrix.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "optimization/Multipliers.hpp"
#include "symbolic/CollectionAdapter.hpp"

namespace local
{

    class DataModel : public uno::Model
    {
    protected:
        std::vector<double> _variable_lower_bounds;
        std::vector<double> _variable_upper_bounds;
        std::vector<double> _constraint_lower_bounds;
        std::vector<double> _constraint_upper_bounds;

        std::vector<uno::BoundType> _variable_status;    /*!< Status of the variables (EQUALITY, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES) */
        std::vector<uno::FunctionType> _constraint_type; /*!< Types of the constraints (LINEAR, QUADRATIC, NONLINEAR) */
        std::vector<uno::BoundType> _constraint_status;  /*!< Status of the constraints (EQUAL_BOUNDS, BOUNDED_LOWER, BOUNDED_UPPER, BOUNDED_BOTH_SIDES,UNBOUNDED) */
        std::vector<size_t> _linear_constraints;
        

        uno::SparseVector<size_t> _slacks{};
    private:
        // lists of variables and constraints + corresponding collection objects
        std::vector<size_t> _equality_constraints;
        std::vector<size_t> _inequality_constraints;
        uno::CollectionAdapter<std::vector<size_t> &> _equality_constraints_collection;
        uno::CollectionAdapter<std::vector<size_t> &> _inequality_constraints_collection;
        std::vector<size_t> _lower_bounded_variables;
        uno::CollectionAdapter<std::vector<size_t> &> _lower_bounded_variables_collection;
        std::vector<size_t> _upper_bounded_variables;
        uno::CollectionAdapter<std::vector<size_t> &> _upper_bounded_variables_collection;
        std::vector<size_t> _single_lower_bounded_variables; // indices of the single lower-bounded variables
        uno::CollectionAdapter<std::vector<size_t> &> _single_lower_bounded_variables_collection;
        std::vector<size_t> _single_upper_bounded_variables; // indices of the single upper-bounded variables
        uno::CollectionAdapter<std::vector<size_t> &> _single_upper_bounded_variables_collection;

    public:
        DataModel(size_t number_variables, size_t number_constraints, const std::string &name = "DataModel")
            : uno::Model(name, number_variables, number_constraints, 1.),
              _variable_lower_bounds(number_variables),
              _variable_upper_bounds(number_variables),
              _constraint_lower_bounds(number_constraints),
              _constraint_upper_bounds(number_constraints),
              _variable_status(number_variables),
              _constraint_type(number_constraints),
              _constraint_status(number_constraints),
              _linear_constraints(0),
              _slacks(0),
              _equality_constraints(0),
              _inequality_constraints(0),
              _equality_constraints_collection(_equality_constraints),
              _inequality_constraints_collection(_inequality_constraints),
              _lower_bounded_variables(0),
              _lower_bounded_variables_collection(_lower_bounded_variables),
              _upper_bounded_variables(0),
              _upper_bounded_variables_collection(_upper_bounded_variables),
              _single_lower_bounded_variables(0),
              _single_lower_bounded_variables_collection(_single_lower_bounded_variables),
              _single_upper_bounded_variables(0),
              _single_upper_bounded_variables_collection(_single_upper_bounded_variables) {}

        virtual ~DataModel() override = default;

        const uno::Collection<size_t> &get_equality_constraints() const override
        {
            return _equality_constraints_collection;
        }

        const uno::Collection<size_t> &get_inequality_constraints() const override
        {
            return _inequality_constraints_collection;
        }

        const std::vector<size_t> &get_linear_constraints() const override
        {
            return _linear_constraints;
        }

        const uno::SparseVector<size_t> &get_slacks() const override
        {
            return _slacks;
        }

        const uno::Collection<size_t> &get_single_lower_bounded_variables() const override
        {
            return _single_lower_bounded_variables_collection;
        }

        const uno::Collection<size_t> &get_single_upper_bounded_variables() const override
        {
            return _single_upper_bounded_variables_collection;
        }

        const uno::Collection<size_t> &get_lower_bounded_variables() const override
        {
            return _lower_bounded_variables_collection;
        }

        const uno::Collection<size_t> &get_upper_bounded_variables() const override
        {
            return _upper_bounded_variables_collection;
        }
        double variable_lower_bound(size_t variable_index) const override { 
            return _variable_lower_bounds[variable_index]; 
            }

        double variable_upper_bound(size_t variable_index) const override { 
            return _variable_upper_bounds[variable_index]; 
        }

        uno::BoundType get_variable_bound_type(size_t variable_index) const override {
            return _variable_status[variable_index];
        }

        double constraint_lower_bound(size_t constraint_index) const override {
            return _constraint_lower_bounds[constraint_index];
        }

        double constraint_upper_bound(size_t constraint_index) const override {
            return _constraint_upper_bounds[constraint_index];
        }
        
        uno::BoundType get_constraint_bound_type(size_t constraint_index) const override {
            return _constraint_status[constraint_index];
        }
        uno::FunctionType get_constraint_type(size_t constraint_index) const override {
            return _constraint_type[constraint_index];
        }

        void initialise_from_data()
        {
            for (std::size_t i = 0; i < number_variables; ++i)
            {
                if (_variable_lower_bounds[i] == _variable_upper_bounds[i])
                {
                    _variable_status[i] = uno::BoundType::EQUAL_BOUNDS;
                    _lower_bounded_variables.emplace_back(i);
                    _upper_bounded_variables.emplace_back(i);
                }
                else if (std::isfinite(_variable_lower_bounds[i]) && std::isfinite(_variable_upper_bounds[i]))
                {
                    _variable_status[i] = uno::BoundType::BOUNDED_BOTH_SIDES;
                    _lower_bounded_variables.emplace_back(i);
                    _upper_bounded_variables.emplace_back(i);
                }
                else if (std::isfinite(_variable_lower_bounds[i]))
                {
                    _variable_status[i] = uno::BoundType::BOUNDED_LOWER;
                    _lower_bounded_variables.emplace_back(i);
                    _single_lower_bounded_variables.emplace_back(i);
                }
                else if (std::isfinite(_variable_upper_bounds[i]))
                {
                    _variable_status[i] = uno::BoundType::BOUNDED_UPPER;
                    _upper_bounded_variables.emplace_back(i);
                    _single_upper_bounded_variables.emplace_back(i);
                }
                else
                {
                    _variable_status[i] = uno::BoundType::UNBOUNDED;
                }
            }

            for (std::size_t i = 0; i<number_constraints; ++i) {
                if (_constraint_lower_bounds[i] == _constraint_upper_bounds[i])
                {
                    _constraint_status[i] = uno::BoundType::EQUAL_BOUNDS;
                    _equality_constraints.emplace_back(i);
                }
                else if (std::isfinite(_constraint_lower_bounds[i]) && std::isfinite(_constraint_upper_bounds[i]))
                {
                    _constraint_status[i] = uno::BoundType::BOUNDED_BOTH_SIDES;
                    _inequality_constraints.emplace_back(i);
                }
                else if (std::isfinite(_constraint_lower_bounds[i]))
                {
                    _constraint_status[i] = uno::BoundType::BOUNDED_LOWER;
                    _inequality_constraints.emplace_back(i);
                }
                else if (std::isfinite(_constraint_upper_bounds[i]))
                {
                    _constraint_status[i] = uno::BoundType::BOUNDED_UPPER;
                    _inequality_constraints.emplace_back(i);
                }
                else
                {
                    _constraint_status[i] = uno::BoundType::UNBOUNDED;
                }
            }
        }
    };

    class HS71
        : public local::DataModel
    {
        public:
        HS71() : local::DataModel(4, 2, "HS71") {
            _variable_lower_bounds = {1.0, 1.0, 1.0, 1.0};
            _variable_upper_bounds = {5.0, 5.0, 5.0, 5.0};
            _constraint_lower_bounds = {25.0, 40.0};
            _constraint_upper_bounds = {std::numeric_limits<double>::infinity(), 40.0};
            _constraint_type = {uno::FunctionType::NONLINEAR, uno::FunctionType::NONLINEAR};
            initialise_from_data();
        }
        
        double evaluate_objective(const uno::Vector<double> &x) const override { return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]; }

        void evaluate_objective_gradient(const uno::Vector<double> &x, uno::SparseVector<double> &gradient) const override
        {
            gradient.insert(0, x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]));
            gradient.insert(1, x[0] * x[3]);
            gradient.insert(2, x[0] * x[3] + 1);
            gradient.insert(3, x[0] * (x[0] + x[1] + x[2]));
        }

        void evaluate_constraints(const uno::Vector<double> &x, std::vector<double> &constraints) const override
        {
            constraints[0] = x[0] * x[1] * x[2] * x[3];
            constraints[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
        }

        void evaluate_constraint_gradient(const uno::Vector<double> &x, size_t constraint_index, uno::SparseVector<double> &gradient) const override
        {
            if (constraint_index == 0)
            {
                gradient.insert(0, x[1] * x[2] * x[3]);
                gradient.insert(1, x[0] * x[2] * x[3]);
                gradient.insert(2, x[0] * x[1] * x[3]);
                gradient.insert(3, x[0] * x[1] * x[2]);
            }
            else if (constraint_index == 1)
            {
                gradient.insert(0, 2 * x[0]);
                gradient.insert(1, 2 * x[1]);
                gradient.insert(2, 2 * x[2]);
                gradient.insert(3, 2 * x[3]);
            }
        }
        void evaluate_constraint_jacobian(const uno::Vector<double> &x, uno::RectangularMatrix<double> &constraint_jacobian) const override
        {
            evaluate_constraint_gradient(x, 0, constraint_jacobian[0]);
            evaluate_constraint_gradient(x, 1, constraint_jacobian[1]);
        }
        void evaluate_lagrangian_hessian(const uno::Vector<double> &x, double objective_multiplier, const uno::Vector<double> &multipliers,
                                         uno::SymmetricMatrix<size_t, double> &hessian) const override
        {
            double H[4][4] = {{2 * x[3]              , x[3], x[3], 2 * x[0] + x[1] + x[2]},
                              {x[3]                  , 0   , 0   , x[0]},
                              {x[3]                  , 0   , 0   , x[0]},
                              {2 * x[0] + x[1] + x[2], x[0], x[0], 0}};

            double cH[4][4] = {{0          , x[2] * x[3], x[1] * x[3], x[1] * x[2]},
                               {x[2] * x[3], 0          , x[0] * x[3], x[0] * x[2]},
                               {x[1] * x[3], x[0] * x[3], 0          , x[0] * x[1]},
                               {x[1] * x[2], x[0] * x[2], x[0] * x[1], 0          }};
            const double cH2[4][4] = {{2, 0, 0, 0},
                                      {0, 2, 0, 0},
                                      {0, 0, 2, 0},
                                      {0, 0, 0, 2}};
            hessian.reset();
            hessian.reset();
            for (size_t jCol = 0; jCol < 4; jCol++)
            {
                for (size_t iRow = jCol; iRow < 4; iRow++)
                {
                    hessian.insert(objective_multiplier * H[iRow][jCol] + multipliers[0] * cH[iRow][jCol] + multipliers[1] * cH2[iRow][jCol], iRow, jCol);
                }
                hessian.finalize_column(jCol);
            }
        }

        void initial_primal_point(uno::Vector<double> &x) const override {
            std::vector<double> x0{3., 3., 3., 3.};
            std::copy(x0.cbegin(), x0.cend(), x.begin());
        }

        void initial_dual_point(uno::Vector<double> &multipliers) const override {
            std::fill(multipliers.begin(), multipliers.end(), 0.0);
        }

        void postprocess_solution(uno::Iterate &iterate, uno::TerminationStatus termination_status) const override {

        }

        size_t number_objective_gradient_nonzeros() const override { return number_variables; }
        size_t number_jacobian_nonzeros() const override { return 2 * number_variables; }
        size_t number_hessian_nonzeros() const override { return number_variables * number_variables; }
        
    };

} // namespace local

#endif // NLPTEST_HS71_HPP