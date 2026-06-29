# Copyright (c) 2026: Charlie Vanaret and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using JuMP, AmplNLWriter, Uno_jll

function Optimizer(options)
    return AmplNLWriter.Optimizer(Uno_jll.amplexe, options)
end

function test_hs015()
    model = Model(() -> Optimizer(["logger=INFO"]))

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

function test_clnlbeam()
    model = Model(() -> Optimizer(["logger=INFO"]))

    ni    = 20000
    alpha = 350.0
    h     = 1 / ni

    @variable(model, -1.0  <= t[i in 0:ni] <= 1.0,  start = 0.05 * cos(i * h))
    @variable(model, -0.05 <= x[i in 0:ni] <= 0.05, start = 0.05 * cos(i * h))
    @variable(model, u[0:ni])     # free, no bounds

    @objective(model, Min, sum(0.5 * h * (u[i+1]^2 + u[i]^2) +
            0.5 * alpha * h * (cos(t[i+1]) + cos(t[i])) for i in 0:ni-1))

    @constraint(model, cons1[i in 0:ni-1],
        x[i+1] - x[i] - 0.5 * h * (sin(t[i+1]) + sin(t[i])) == 0)

    @constraint(model, cons2[i in 0:ni-1],
        t[i+1] - t[i] - 0.5 * h * u[i+1] - 0.5 * h * u[i] == 0)

    # fix boundary values (force=true overrides existing bounds)
    fix(x[0],  0.0; force = true)
    fix(x[ni], 0.0; force = true)
    fix(t[0],  0.0; force = true)
    fix(t[ni], 0.0; force = true)

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

test_hs015()
test_clnlbeam()