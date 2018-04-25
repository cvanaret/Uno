g3 0 1 0	# problem ex8_1_7
 5 5 1 0 1	# vars, constraints, objectives, ranges, eqns
 5 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 4 5 4	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 14 5	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
0 -5 5
0 -5 5
0 -5 5
0 -5 5
0 -5 5
r
1 6.24264068711929
1 -6.24264068711929
1 0.82842712474619
1 -0.82842712474619
4 2
C0
o0
o5
v1
n2
o5
v2
n3
C1
o1
o16
o5
v1
n2
o5
v2
n3
C2
o16
o5
v2
n2
C3
o5
v2
n2
C4
o0
o2
o2
n0.5
v0
v3
o2
o2
n0.5
v0
v3
O0 0
o54
5
o5
o0
n-1
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
3
7
11
12
J0 3
0 1
1 0
2 0
J1 3
0 -1
1 0
2 0
J2 3
1 1
2 0
4 1
J3 3
1 -1
2 0
4 -1
J4 2
0 0
3 0
G0 5
0 0
1 0
2 0
3 0
4 0
