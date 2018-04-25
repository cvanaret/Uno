g3 0 1 0	# problem nactive
 2 3 1 0 0	# vars, constraints, objectives, ranges, eqns
 3 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 1 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 6 1	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
x2
0 10
1 -20
r
2 0.5
2 0
2 0
C0
o2
n-0.5
o5
v0
n2
C1
o16
o5
v0
n2
C2
o5
v0
n2
O0 0
n0
k1
3
J0 2
0 0
1 -0.5
J1 2
0 0
1 1
J2 2
0 0
1 -1
G0 1
1 1
