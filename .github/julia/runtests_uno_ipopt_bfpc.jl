# Copyright (c) 2026: Charlie Vanaret and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using JuMP, AmplNLWriter, Uno_jll

function Optimizer(options)
    return AmplNLWriter.Optimizer(Uno_jll.amplexe, options)
end

# this script is parameterized by the linear solver passed as command line argument (e.g., runtests_uno_ipopt_bfpc.jl MUMPS)
length(ARGS) == 1 || error("The linear solver name is missing or you supplied more than one command line arguments")
linear_solver = ARGS[1]
println("Solving with linear solver ", linear_solver)

 Optimizer_Uno() = Optimizer(["logger=SILENT", "preset=ipopt", "linear_solver=$linear_solver","unbounded_objective_threshold=-1e15"])

function test_bfpc()
    # bfpc instance described here: https://github.com/jump-dev/Ipopt.jl/pull/546
    # this is a challenging instance numerically. Linear solvers must be equipped
    # with a robust dual regularization procedure: estimated zero eigenvalues may
    # appear *during* the regularization loop.
    model = Model(() ->  Optimizer_Uno())

    @variable(model, volume_0_1 == 100)  # volume[0,1]
    @variable(model, 0 <= volume_1_1 <= 300)  # volume[1,1]
    @variable(model, 0 <= volume_2_1 <= 300)  # volume[2,1]
    @variable(model, volume_0_2 == 100)  # volume[0,2]
    @variable(model, 0 <= volume_1_2 <= 300)  # volume[1,2]
    @variable(model, 0 <= volume_2_2 <= 300)  # volume[2,2]
    @variable(model, 0 <= hydro_1_1 <= 100)  # hydro[1,1]
    @variable(model, 0 <= hydro_2_1 <= 100)  # hydro[2,1]
    @variable(model, 0 <= hydro_1_2 <= 100)  # hydro[1,2]
    @variable(model, 0 <= hydro_2_2 <= 100)  # hydro[2,2]
    @variable(model, spill_1_1 >= 0)  # spill[1,1]
    @variable(model, spill_2_1 >= 0)  # spill[2,1]
    @variable(model, spill_1_2 >= 0)  # spill[1,2]
    @variable(model, spill_2_2 >= 0)  # spill[2,2]
    @variable(model, 0 <= termos_1_1_1 <= 10)  # termos[1,1,1]
    @variable(model, 0 <= termos_2_1_1 <= 5)  # termos[2,1,1]
    @variable(model, 0 <= termos_3_1_1 <= 5)  # termos[3,1,1]
    @variable(model, 0 <= termos_1_2_1 <= 10)  # termos[1,2,1]
    @variable(model, 0 <= termos_2_2_1 <= 5)  # termos[2,2,1]
    @variable(model, 0 <= termos_3_2_1 <= 5)  # termos[3,2,1]
    @variable(model, 0 <= termos_1_1_2 <= 10)  # termos[1,1,2]
    @variable(model, 0 <= termos_2_1_2 <= 5)  # termos[2,1,2]
    @variable(model, 0 <= termos_3_1_2 <= 5)  # termos[3,1,2]
    @variable(model, 0 <= termos_1_2_2 <= 10)  # termos[1,2,2]
    @variable(model, 0 <= termos_2_2_2 <= 5)  # termos[2,2,2]
    @variable(model, 0 <= termos_3_2_2 <= 5)  # termos[3,2,2]
    @variable(model, theta_1 >= 0)  # θ[1]
    @variable(model, theta_2 >= 0)  # θ[2]
    @variable(model, volume_ldr_1_1)  # volume_ldr[1,1]
    @variable(model, volume_ldr_2_1)  # volume_ldr[2,1]
    @variable(model, volume_ldr_1_2)  # volume_ldr[1,2]
    @variable(model, volume_ldr_2_2)  # volume_ldr[2,2]
    @variable(model, hydro_ldr_1_1)  # hydro_ldr[1,1]
    @variable(model, hydro_ldr_2_1)  # hydro_ldr[2,1]
    @variable(model, hydro_ldr_1_2)  # hydro_ldr[1,2]
    @variable(model, hydro_ldr_2_2)  # hydro_ldr[2,2]
    @variable(model, spill_ldr_1_1)  # spill_ldr[1,1]
    @variable(model, spill_ldr_2_1)  # spill_ldr[2,1]
    @variable(model, spill_ldr_1_2)  # spill_ldr[1,2]
    @variable(model, spill_ldr_2_2)  # spill_ldr[2,2]
    @variable(model, termos_ldr_1_1_1)  # termos_ldr[1,1,1]
    @variable(model, termos_ldr_2_1_1)  # termos_ldr[2,1,1]
    @variable(model, termos_ldr_3_1_1)  # termos_ldr[3,1,1]
    @variable(model, termos_ldr_1_2_1)  # termos_ldr[1,2,1]
    @variable(model, termos_ldr_2_2_1)  # termos_ldr[2,2,1]
    @variable(model, termos_ldr_3_2_1)  # termos_ldr[3,2,1]
    @variable(model, termos_ldr_1_1_2)  # termos_ldr[1,1,2]
    @variable(model, termos_ldr_2_1_2)  # termos_ldr[2,1,2]
    @variable(model, termos_ldr_3_1_2)  # termos_ldr[3,1,2]
    @variable(model, termos_ldr_1_2_2)  # termos_ldr[1,2,2]
    @variable(model, termos_ldr_2_2_2)  # termos_ldr[2,2,2]
    @variable(model, termos_ldr_3_2_2)  # termos_ldr[3,2,2]
    @variable(model, theta_ldr_1)  # θ_ldr[1]
    @variable(model, theta_ldr_2)  # θ_ldr[2]
     
    @objective(model, Min, 125*termos_1_1_1 + 125*termos_1_1_2 + 125*termos_1_2_1 + 125*termos_1_2_2 + 375*termos_2_1_1 + 375*termos_2_1_2 + 375*termos_2_2_1 + 375*termos_2_2_2 + 875*termos_3_1_1 + 875*termos_3_1_2 + 875*termos_3_2_1 + 875*termos_3_2_2 + 0.5*theta_1 + 0.5*theta_2)
     
    # Constraints
    @constraint(model, c1, - volume_0_1 + volume_1_1 + hydro_1_1 + spill_1_1 == 60)
    @constraint(model, c2, hydro_1_1 + termos_1_1_1 + termos_2_1_1 + termos_3_1_1 == 70)
    @constraint(model, c3, - volume_0_2 + volume_1_2 + hydro_1_2 + spill_1_2 == 60)
    @constraint(model, c4, hydro_1_2 + termos_1_1_2 + termos_2_1_2 + termos_3_1_2 == 70)
    @constraint(model, c5, - volume_1_1 + volume_2_1 + hydro_2_1 + spill_2_1 == 0)
    @constraint(model, c6, hydro_2_1 + termos_1_2_1 + termos_2_2_1 + termos_3_2_1 == 70)
    @constraint(model, c7, - volume_1_2 + volume_2_2 + hydro_2_2 + spill_2_2 == 120)
    @constraint(model, c8, hydro_2_2 + termos_1_2_2 + termos_2_2_2 + termos_3_2_2 == 70)
    @constraint(model, c9, volume_1_1 - volume_ldr_1_1 - 60*volume_ldr_1_2 == 0)
    @constraint(model, c10, volume_1_2 - volume_ldr_1_1 - 60*volume_ldr_1_2 == 0)
    @constraint(model, c11, hydro_1_1 - hydro_ldr_1_1 - 60*hydro_ldr_1_2 == 0)
    @constraint(model, c12, hydro_1_2 - hydro_ldr_1_1 - 60*hydro_ldr_1_2 == 0)
    @constraint(model, c13, spill_1_1 - spill_ldr_1_1 - 60*spill_ldr_1_2 == 0)
    @constraint(model, c14, spill_1_2 - spill_ldr_1_1 - 60*spill_ldr_1_2 == 0)
    @constraint(model, c15, termos_1_1_1 - termos_ldr_1_1_1 - 60*termos_ldr_1_1_2 == 0)
    @constraint(model, c16, termos_1_1_2 - termos_ldr_1_1_1 - 60*termos_ldr_1_1_2 == 0)
    @constraint(model, c17, termos_2_1_1 - termos_ldr_2_1_1 - 60*termos_ldr_2_1_2 == 0)
    @constraint(model, c18, termos_2_1_2 - termos_ldr_2_1_1 - 60*termos_ldr_2_1_2 == 0)
    @constraint(model, c19, termos_3_1_1 - termos_ldr_3_1_1 - 60*termos_ldr_3_1_2 == 0)
    @constraint(model, c20, termos_3_1_2 - termos_ldr_3_1_1 - 60*termos_ldr_3_1_2 == 0)
    @constraint(model, c21, volume_2_1 - volume_ldr_2_1 == 0)
    @constraint(model, c22, volume_2_2 - volume_ldr_2_1 - 120*volume_ldr_2_2 == 0)
    @constraint(model, c23, hydro_2_1 - hydro_ldr_2_1 == 0)
    @constraint(model, c24, hydro_2_2 - hydro_ldr_2_1 - 120*hydro_ldr_2_2 == 0)
    @constraint(model, c25, spill_2_1 - spill_ldr_2_1 == 0)
    @constraint(model, c26, spill_2_2 - spill_ldr_2_1 - 120*spill_ldr_2_2 == 0)
    @constraint(model, c27, termos_1_2_1 - termos_ldr_1_2_1 == 0)
    @constraint(model, c28, termos_1_2_2 - termos_ldr_1_2_1 - 120*termos_ldr_1_2_2 == 0)
    @constraint(model, c29, termos_2_2_1 - termos_ldr_2_2_1 == 0)
    @constraint(model, c30, termos_2_2_2 - termos_ldr_2_2_1 - 120*termos_ldr_2_2_2 == 0)
    @constraint(model, c31, termos_3_2_1 - termos_ldr_3_2_1 == 0)
    @constraint(model, c32, termos_3_2_2 - termos_ldr_3_2_1 - 120*termos_ldr_3_2_2 == 0)
    @constraint(model, c33, theta_1 - theta_ldr_1 == 0)
    @constraint(model, c34, theta_2 - theta_ldr_1 - 120*theta_ldr_2 == 0)
    @constraint(model, c1_1, 300*volume_2_1 + theta_1 >= 75000)
    @constraint(model, c2_1, 600*volume_2_1 + theta_1 >= 135000)
    @constraint(model, c3_1, 1200*volume_2_1 + theta_1 >= 210000)
    @constraint(model, c4_1, 2500*volume_2_1 + theta_1 >= 310000)
    @constraint(model, c5_1, 7800*volume_2_1 + theta_1 >= 390000)
    @constraint(model, c6_1, 300*volume_2_2 + theta_2 >= 75000)
    @constraint(model, c7_1, 600*volume_2_2 + theta_2 >= 135000)
    @constraint(model, c8_1, 1200*volume_2_2 + theta_2 >= 210000)
    @constraint(model, c9_1, 2500*volume_2_2 + theta_2 >= 310000)
    @constraint(model, c10_1, 7800*volume_2_2 + theta_2 >= 390000)
    @constraint(model, c1_2, volume_ldr_1_2 == 0)
    @constraint(model, c2_2, hydro_ldr_1_2 == 0)
    @constraint(model, c3_2, spill_ldr_1_2 == 0)
    @constraint(model, c4_2, termos_ldr_1_1_2 == 0)
    @constraint(model, c5_2, termos_ldr_2_1_2 == 0)
    @constraint(model, c6_2, termos_ldr_3_1_2 == 0)

    optimize!(model)

    tolerance = 1e-6
    @assert abs(objective_value(model) - 119250) <= tolerance
end

test_bfpc()