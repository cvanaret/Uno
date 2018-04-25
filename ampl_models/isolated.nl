g3 0 1 0	# problem isolated
 2 4 1 0 0	# vars, constraints, objectives, ranges, eqns
 4 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 2 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 8 2	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
x2
0 3
1 2
r
2 1
2 1
2 1
2 1
C0
o16
o5
v0
n2
C1
o16
o5
v0
n2
C2
o16
o5
v1
n2
C3
o16
o5
v1
n2
O0 0
n0
k1
4
J0 2
0 0
1 1
J1 2
0 0
1 -1
J2 2
0 1
1 0
J3 2
0 -1
1 0
G0 2
0 1
1 1
