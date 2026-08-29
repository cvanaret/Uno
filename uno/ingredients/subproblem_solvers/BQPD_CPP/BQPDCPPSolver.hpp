// Copyright (c) 2026 Uno contributors
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BQPDCPPSOLVER_H
#define UNO_BQPDCPPSOLVER_H

#include <memory>
#include <vector>
#include "../QPSolver.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "bqpd/Context.hpp"
#include "bqpd/Workspace.hpp"

namespace uno {
   // forward declarations
   class BQPDQuadraticProgram;
   class Multipliers;
   class Options;
   class QuadraticProgram;

   // QP subproblem solver that drives the C++ port BQPD-CPP (SparseA + SparseFactor), so it
   // can coexist with the Fortran BQPDSolver in the same binary. BQPD-CPP exposes no Fortran
   // symbols (no bqpd_/gdotx_), so there is no collision with libbqpd.a.
   //
   // The BQPD-native QP data (weak-CSR gradients/Jacobian, concatenated bounds, COO Hessian)
   // is packed by the SAME BQPDQuadraticProgram the Fortran solver uses, so both solvers
   // receive byte-identical data; only the solve engine differs. The Hessian reaches the port
   // through a HessianProduct adapter wrapping BQPDQuadraticProgram::compute_hessian_vector_product,
   // replacing the Fortran gdotx callback and its lws-pointer smuggle.
   //
   // NOTE (deliberate limitation): BQPD-CPP clamps warm-start modes 3-6 to mode 2, so on a
   // rejected trust-region step (where determine_mode would ask for mode 4) it re-solves the
   // same QP from the same active set instead of reusing factors. The QP optimum is unchanged;
   // pivot/iteration counts may differ from the Fortran leg. See PLAN-2026-08-29-Slice7.
   class BQPDCPPSolver : public QPSolver {
   public:
      explicit BQPDCPPSolver(const Options& options);
      ~BQPDCPPSolver() override;

      void initialize_memory(const Subproblem& subproblem) override;

      [[nodiscard]] QuadraticProgram& get_quadratic_program() override;

      void solve(Statistics& statistics, const Vector<double>& initial_point, Direction& direction,
         const WarmstartInformation& warmstart_information) override;

      [[nodiscard]] SolverWorkspace& get_workspace() override;

   private:
      // BQPD-native quadratic program (built by IQPSolver before each solve), reused for data packing
      std::unique_ptr<BQPDQuadraticProgram> quadratic_program{};

      int kmax{0};
      const int mlp{1000};
      std::size_t nprof{};
      std::vector<double> alp{};
      std::vector<int> lp{}, active_set{}; // active_set = ls, persisted across solves (warm start)
      std::vector<double> gradient_solution{}, residuals{}, w{}, e{};
      int k{0};   // reduced-space dimension, persisted across solves (in/out of solve)
      int peq{0};
      const double fmin{-1e10};

      // BQPD-CPP workspace + shared former-COMMON context, persisted across solves
      bqpd::Workspace workspace{};
      bqpd::SolverContext ctx{};
      int current_mxws{0}, current_mxlws{0};

      const bool print_subproblem;

      void allocate_workspace(std::size_t number_variables, std::size_t number_constraints,
         std::size_t number_jacobian_nonzeros, bool is_qp);
      void set_multipliers(std::size_t number_variables, Multipliers& direction_multipliers) const;
   };
} // namespace

#endif // UNO_BQPDCPPSOLVER_H
