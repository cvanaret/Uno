g3 0 1 0	# problem hs040
 4 3 1 0 3	# vars, constraints, objectives, ranges, eqns
 3 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 3 4 3	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 7 4	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
3
3
3
3
x4
0 0.8
1 0.8
2 0.8
3 0.8
r
4 1
4 0
4 0
C0
o0
o5
v0
n3
o5
v1
n2
C1
o2
o5
v0
n2
v2
C2
o5
v2
n2
O0 0
o16
o2
o2
o2
v0
v1
v3
v2
k3
2
4
6
J0 2
0 0
1 0
J1 3
0 0
2 0
3 -1
J2 2
1 -1
2 0
G0 4
0 0
1 0
2 0
3 0
