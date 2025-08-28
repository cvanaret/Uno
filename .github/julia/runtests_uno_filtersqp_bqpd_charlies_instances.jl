# Copyright (c) 2025: Charlie Vanaret and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using JuMP, AmplNLWriter, Uno_jll

function Optimizer(options)
    return AmplNLWriter.Optimizer(Uno_jll.amplexe, options)
end

Optimizer_Uno_filtersqp() = Optimizer(["logger=SILENT", "preset=filtersqp", "QP_solver=BQPD", "max_iterations=10000", "unbounded_objective_threshold=-1e15"])

function test_hs015()
    model = Model(() -> Optimizer_Uno_filtersqp())

    @variable(model, x, start = -2.)
    @variable(model, y, start = 1.)

    @objective(model, Min, 100*(y - x^2)^2 + (1 - x)^2)

    c1 = @constraint(model, x*y >= 1)
    c2 = @constraint(model, x + y^2 >= 0)
    c3 = @constraint(model, x <= 0.5)

    optimize!(model)

    tolerance = 1e-6
    @assert abs(objective_value(model) - 306.5) <= tolerance
    @assert abs(value(x) - 0.5) <= tolerance
    @assert abs(value(y) - 2.) <= tolerance
    @assert abs(dual(c1) - 700.) <= tolerance
    @assert abs(dual(c2) - 0.) <= tolerance
    @assert abs(dual(c3) - (-1751.)) <= tolerance
end

function test_isolated()
    model = Model(() -> Optimizer_Uno_filtersqp())

    x0 = [3., 2.]
    @variable(model, x[i=1:2], start = x0[i])

    @objective(model, Min, x[1] + x[2])

    c1 = @constraint(model, -x[1]^2 + x[2] - 1 >= 0)
    c2 = @constraint(model, -x[1]^2 - x[2] - 1 >= 0)
    c3 = @constraint(model, x[1] - x[2]^2 - 1 >= 0)
    c4 = @constraint(model, -x[1] - x[2]^2 - 1 >= 0)

    optimize!(model)

    tolerance = 1e-6
    @assert abs(objective_value(model) - 0.) <= tolerance
    @assert abs(value(x[1]) - 0.) <= tolerance
    @assert abs(value(x[2]) - 0.) <= tolerance
    @assert abs(dual(c1) - 1.) <= tolerance
    @assert abs(dual(c2) - 1.) <= tolerance
    @assert abs(dual(c3) - 1.) <= tolerance
    @assert abs(dual(c4) - 1.) <= tolerance
end

function test_nactive()
    model = Model(() -> Optimizer_Uno_filtersqp())

    x0 = [-20., 10.]
    @variable(model, x[i=1:2], start = x0[i])

    @objective(model, Min, x[1])

    c1 = @constraint(model, 0.5*(-x[1] - x[2]^2 - 1) >= 0)
    c2 = @constraint(model, x[1] - x[2]^2 >= 0)
    c3 = @constraint(model, -x[1] + x[2]^2 >= 0)

    optimize!(model)

    tolerance = 1e-6
    @assert abs(objective_value(model) - 0.) <= tolerance
    @assert abs(value(x[1]) - 0.) <= tolerance
    @assert abs(value(x[2]) - 0.) <= tolerance
    @assert abs(dual(c1) - 1.) <= tolerance
    @assert abs(dual(c2) - 0.5) <= tolerance
    @assert abs(dual(c3) - 0.) <= tolerance
end

function test_unique()
    model = Model(() -> Optimizer_Uno_filtersqp())

    x0 = [3., 2.]
    @variable(model, x[i=1:2], start = x0[i])

    @objective(model, Min, x[1] + x[2])

    c1 = @constraint(model, x[2] - x[1]^2 - 1 >= 0)
    c2 = @constraint(model, 0.3*(1 - exp(x[2])) >= 0)

    optimize!(model)

    tolerance = 1e-6
    @assert abs(objective_value(model) - 1.) <= tolerance
    @assert abs(value(x[1]) - 0.) <= tolerance
    @assert abs(value(x[2]) - 1.) <= tolerance
    @assert abs(dual(c1) - 0.8154845) <= tolerance
    @assert abs(dual(c2) - 1.) <= tolerance
end

test_hs015()
test_isolated()
test_nactive()
test_unique()