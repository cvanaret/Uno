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

function test_camshape_6400()
    model = Model(() -> Optimizer(["logger=INFO"]))

    n       = 6400                       # number of discretization points
    R_v     = 1.0                        # design parameter related to the valve shape
    R_min   = 1.0                        # minimum allowed radius of the cam
    R_max   = 2.0                        # maximum allowed radius of the cam
    alpha   = 1.5                        # curvature limit parameter
    d_theta = 2 * pi / (5 * (n + 1))     # angle between discretization points
    cos_dt  = cos(d_theta)               # constant: d_theta is a parameter, not a variable

    # radius of the cam at discretization points; start at the circle of radius (R_min+R_max)/2
    @variable(model, R_min <= r[1:n] <= R_max, start = (R_min + R_max) / 2)

    @objective(model, Max, (pi * R_v / n) * sum(r[i] for i in 1:n))

    # Convexity (bilinear in r; cos_dt is a precomputed constant)
    @constraint(model, convexity[i in 2:n-1],
        -r[i-1] * r[i] - r[i] * r[i+1] + 2 * r[i-1] * r[i+1] * cos_dt <= 0)
    @constraint(model, convex_edge1,
        -R_min * r[1] - r[1] * r[2] + 2 * R_min * r[2] * cos_dt <= 0)
    @constraint(model, convex_edge2,
        -R_min^2 - R_min * r[1] + 2 * R_min * r[1] * cos_dt <= 0)
    @constraint(model, convex_edge3,
        -r[n-1] * r[n] - r[n] * R_max + 2 * r[n-1] * R_max * cos_dt <= 0)
    @constraint(model, convex_edge4,
        -2 * R_max * r[n] + 2 * r[n]^2 * cos_dt <= 0)

    # Curvature, lower bound
    @constraint(model, curvature[i in 1:n-1], -alpha * d_theta <= r[i+1] - r[i])
    @constraint(model, curvature_edge1,       -alpha * d_theta <= r[1] - R_min)
    @constraint(model, curvature_edge2,       -alpha * d_theta <= R_max - r[n])

    # Curvature, upper bound
    @constraint(model, curvature1[i in 1:n-1], r[i+1] - r[i] <= alpha * d_theta)
    @constraint(model, curvature_edge11,       r[1] - R_min  <= alpha * d_theta)
    @constraint(model, curvature_edge21,       R_max - r[n]  <= alpha * d_theta)


    optimize!(model)

    tolerance = 1e-6
    @assert abs(objective_value(model) - 4.448378) <= tolerance
end

test_hs015()
test_camshape_6400()