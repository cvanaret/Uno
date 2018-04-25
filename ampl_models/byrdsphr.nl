g3 0 1 0	# problem byrdsphr
 3 2 1 0 2	# vars, constraints, objectives, ranges, eqns
 2 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 6 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
x3
0 5
1 0.0001
2 -0.0001
r
4 9
4 9
C0
o54
3
o5
v0
n2
o5
v1
n2
o5
v2
n2
C1
o54
3
o5
o0
n-1
v0
n2
o5
v1
n2
o5
v2
n2
O0 0
n0
k2
2
4
J0 3
0 0
1 0
2 0
J1 3
0 0
1 0
2 0
G0 3
0 -1
1 -1
2 -1
