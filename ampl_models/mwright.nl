g3 0 1 0	# problem mwright
 5 3 1 0 3	# vars, constraints, objectives, ranges, eqns
 3 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 4 5 4	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 8 5	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
3
3
x5
0 -1
1 2
2 1
3 -2
4 -2
r
4 6.242640687119286
4 0.8284271247461903
4 2
C0
o0
o5
v1
n2
o5
v2
n2
C1
o16
o5
v2
n2
C2
o2
v0
v3
O0 0
o54
5
o5
v0
n2
o5
o1
v0
v1
n2
o5
o1
v1
v2
n3
o5
o1
v2
v4
n4
o5
o1
v4
v3
n4
k4
2
4
6
7
J0 3
0 1
1 0
2 0
J1 3
1 1
2 0
4 1
J2 2
0 0
3 0
G0 5
0 0
1 0
2 0
3 0
4 0
