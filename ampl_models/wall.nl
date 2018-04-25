g3 0 1 0	# problem wall
 6 6 1 0 6	# vars, constraints, objectives, ranges, eqns
 4 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 6 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 20 1	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
3
3
3
x6
0 1
1 1
2 1
3 1
4 1
5 1
r
4 1
4 4.8
4 0.98
4 1
4 0
4 0
C0
o2
v0
v1
C1
o3
o3
v2
v0
v3
C2
o3
o3
v4
v1
v5
C3
o2
v5
v3
C4
n0
C5
n0
O0 0
n0
k5
4
8
11
14
17
J0 2
0 0
1 0
J1 3
0 0
2 0
3 0
J2 3
1 0
4 0
5 0
J3 2
3 0
5 0
J4 4
0 1
1 -1
2 1e-07
4 -1e-05
J5 6
0 2
1 -2
2 1e-07
3 -0.01
4 -1e-05
5 0.01
G0 1
0 1
