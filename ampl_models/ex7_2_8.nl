g3 0 1 0	# problem ex7_2_8
 8 4 1 0 0	# vars, constraints, objectives, ranges, eqns
 4 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 6 8 2	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 13 4	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
0 0.01 10
0 0.01 10
0 0.01 10
0 0.01 10
0 0.01 10
0 0.01 10
0 0.01 10
0 0.01 10
x8
0 0.01
1 0.01
2 0.01
3 0.01
4 0.01
5 0.01
6 0.01
7 0.01
r
1 1
1 1
1 1
1 1
C0
o2
o2
n0.0588
v0
v3
C1
o2
o2
n0.0588
v1
v5
C2
o54
3
o3
o2
n4
v2
v3
o3
n2
o2
o5
v2
n0.71
v3
o3
o2
n0.0588
v0
o5
v2
n1.3
C3
o54
3
o3
o2
n4
v4
v5
o3
n2
o2
o5
v4
n0.71
v5
o2
o2
n0.0588
o5
v4
n1.3
v1
O0 0
o0
o3
o2
n0.4
o5
v6
n0.67
o5
v0
n0.67
o3
o2
n0.4
o5
v7
n0.67
o5
v1
n0.67
k7
2
4
5
7
8
10
12
J0 3
0 0
3 0
6 0.1
J1 4
1 0
5 0
6 0.1
7 0.1
J2 3
0 0
2 0
3 0
J3 3
1 0
4 0
5 0
G0 4
0 0
1 0
6 -1
7 -1
