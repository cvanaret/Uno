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

    using namespace uno;

    class Options;

    class HS71
        : public Model
    {
        public:
        HS71() : Model("HS71", 4, 2, 1.0) {}
        
        double evaluate_objective(const Vector<double> &x) const override { return x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2]; }

        void evaluate_objective_gradient(const Vector<double> &x, SparseVector<double> &gradient) const override
        {
            gradient.insert(0, x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]));
            gradient.insert(1, x[0] * x[3]);
            gradient.insert(2, x[0] * x[3] + 1);
            gradient.insert(3, x[0] * (x[0] + x[1] + x[2]));
        }

        void evaluate_constraints(const Vector<double> &x, std::vector<double> &constraints) const override
        {
            constraints[0] = x[0] * x[1] * x[2] * x[3];
            constraints[1] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];
        }

        void evaluate_constraint_gradient(const Vector<double> &x, size_t constraint_index, SparseVector<double> &gradient) const override
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
        void evaluate_constraint_jacobian(const Vector<double> &x, RectangularMatrix<double> &constraint_jacobian) const override
        {
            evaluate_constraint_gradient(x, 0, constraint_jacobian[0]);
            evaluate_constraint_gradient(x, 1, constraint_jacobian[1]);
        }
        void evaluate_lagrangian_hessian(const Vector<double> &x, double objective_multiplier, const Vector<double> &multipliers,
                                         SymmetricMatrix<size_t, double> &hessian) const override
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

        double variable_lower_bound(size_t variable_index) const override { return x_l[variable_index]; }
        double variable_upper_bound(size_t variable_index) const override { return x_u[variable_index]; }
        BoundType get_variable_bound_type(size_t variable_index) const override {
            if (x_l[variable_index] == x_u[variable_index])
            {
                return BoundType::EQUAL_BOUNDS;
            }
            else if (std::isfinite(x_l[variable_index]) && std::isfinite(x_u[variable_index]))
            {
                return BoundType::BOUNDED_BOTH_SIDES;
            }
            else if (std::isfinite(x_l[variable_index]))
            {
                return BoundType::BOUNDED_LOWER;
            }
            else if (std::isfinite(x_u[variable_index]))
            {
                return BoundType::BOUNDED_UPPER;
            }
            else
            {
                return BoundType::UNBOUNDED;
            }
        }
        
        const Collection<size_t> &get_lower_bounded_variables() const override {
            static std::vector<size_t> bLower;
            for (size_t i = 0; i < x_l.size(); ++i){
                if (std::isfinite(x_l[i]))
                    bLower.emplace_back(i);

            }
            static CollectionAdapter<std::vector<size_t> &> lower_bounded_variables_collection(bLower);
            return lower_bounded_variables_collection;
        }

        const Collection<size_t> &get_upper_bounded_variables() const override {
            static std::vector<size_t> bUpper;
            for (size_t i = 0; i < x_u.size(); ++i){
                if (std::isfinite(x_u[i]))
                    bUpper.emplace_back(i);

            }
            static CollectionAdapter<std::vector<size_t> &> upper_bounded_variables_collection(bUpper);
            return upper_bounded_variables_collection;
        }

        const SparseVector<size_t> &get_slacks() const override {
            static SparseVector<size_t> slacks_vector{};
            return slacks_vector;
        }

        const Collection<size_t> &get_single_lower_bounded_variables() const override {
            static std::vector<size_t> bSingleLower;
            for (size_t i = 0; i < x_l.size(); ++i){
                if (std::isfinite(x_l[i]) && !std::isfinite(x_u[i]))
                    bSingleLower.emplace_back(i);

            }
            static CollectionAdapter<std::vector<size_t> &> single_lower_bounded_variables_collection(bSingleLower);
            return single_lower_bounded_variables_collection;
        }

        const Collection<size_t> &get_single_upper_bounded_variables() const override {
            static std::vector<size_t> bSingleUpper;
            for (size_t i = 0; i < x_u.size(); ++i){
                if (!std::isfinite(x_l[i]) && std::isfinite(x_u[i]))
                    bSingleUpper.emplace_back(i);

            }
            static CollectionAdapter<std::vector<size_t> &> single_upper_bounded_variables_collection(bSingleUpper);
            return single_upper_bounded_variables_collection;
        }

        double constraint_lower_bound(size_t constraint_index) const override {
            return g_l[constraint_index];
        }

        double constraint_upper_bound(size_t constraint_index) const override {
            return g_u[constraint_index];
        }
        FunctionType get_constraint_type(size_t constraint_index) const override {
            return FunctionType::NONLINEAR;
        }
        BoundType get_constraint_bound_type(size_t constraint_index) const override {
            if (g_l[constraint_index] == g_u[constraint_index])
            {
                return BoundType::EQUAL_BOUNDS;
            }
            else if (std::isfinite(g_l[constraint_index]) && std::isfinite(g_u[constraint_index]))
            {
                return BoundType::BOUNDED_BOTH_SIDES;
            }
            else if (std::isfinite(g_l[constraint_index]))
            {
                return BoundType::BOUNDED_LOWER;
            }
            else if (std::isfinite(g_u[constraint_index]))
            {
                return BoundType::BOUNDED_UPPER;
            }
            else
            {
                return BoundType::UNBOUNDED;
            }
        }
        const Collection<size_t> &get_equality_constraints() const override {
            static std::vector<size_t> eqConstraints;
            for (size_t i = 0; i < g_l.size(); ++i){
                if (std::isfinite(g_l[i]) && std::isfinite(g_u[i]) && g_l[i] == g_u[i])
                    eqConstraints.emplace_back(i);
            }
            static CollectionAdapter<std::vector<size_t> &> equality_constraints_collection(eqConstraints);
            return equality_constraints_collection;
        }
        const Collection<size_t> &get_inequality_constraints() const override {
            static std::vector<size_t> ineqConstraints;
            for (size_t i = 0; i < g_l.size(); ++i){
                if (g_l[i] != g_u[i])
                    ineqConstraints.emplace_back(i);
            }
            static CollectionAdapter<std::vector<size_t> &> inequality_constraints_collection(ineqConstraints);
            return inequality_constraints_collection;
        }
        const std::vector<size_t> &get_linear_constraints() const override {
            static std::vector<size_t> linConstraints(0);
            return linConstraints;
        }

        void initial_primal_point(Vector<double> &x) const override {
            std::copy(x0.cbegin(), x0.cend(), x.begin());
        }

        void initial_dual_point(Vector<double> &multipliers) const override {
            std::fill(multipliers.begin(), multipliers.end(), 0.0);
        }

        void postprocess_solution(Iterate &iterate, TerminationStatus termination_status) const override {

        }

        size_t number_objective_gradient_nonzeros() const override { return x0.size(); }
        size_t number_jacobian_nonzeros() const override { return 2 * x0.size(); }
        size_t number_hessian_nonzeros() const override { return x0.size() * x0.size(); }

    private:
        std::vector<double> x_l{1.0, 1.0, 1.0, 1.0};
        std::vector<double> x_u{5.0, 5.0, 5.0, 5.0};
        std::vector<double> g_l{25.0, 40.0};
        std::vector<double> g_u{std::numeric_limits<double>::infinity(), 40.0};
        std::vector<double> x0{3., 3., 3., 3.};
    };

} // namespace local

#endif // NLPTEST_HS71_HPP