g3 0 1 0	# problem makela2
 3 3 1 0 0	# vars, constraints, objectives, ranges, eqns
 3 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 2 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 9 1	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
x2
0 -1
1 5
r
1 0
1 -40
1 -60
C0
o0
o5
v0
n2
o5
v1
n2
C1
o0
o5
v0
n2
o5
v1
n2
C2
o0
o5
v0
n2
o5
v1
n2
O0 0
n0
k2
3
6
J0 3
0 0
1 0
2 -1
J1 3
0 -40
1 -10
2 -1
J2 3
0 -10
1 -20
2 -1
G0 1
2 1
