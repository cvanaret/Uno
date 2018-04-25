g3 0 1 0	# problem minmaxrb
 3 4 1 0 0	# vars, constraints, objectives, ranges, eqns
 2 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 1 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 10 1	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
x3
0 -1.2
1 1
2 1
r
2 0
2 0
2 1
2 -1
C0
o2
n10
o5
v0
n2
C1
o2
n-10
o5
v0
n2
C2
n0
C3
n0
O0 0
n0
k2
4
6
J0 3
0 0
1 -10
2 1
J1 3
0 0
1 10
2 1
J2 2
0 1
2 1
J3 2
0 -1
2 1
G0 1
2 1
