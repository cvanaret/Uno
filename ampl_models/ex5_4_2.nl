g3 4 1 0	# problem ex5_4_2
 8 6 1 0 0	# vars, constraints, objectives, ranges, eqns
 3 0	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 8 0 0	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 17 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
b
0 100 10000
0 1000 10000
0 1000 10000
0 10 1000
0 10 1000
0 10 1000
0 10 1000
0 10 1000
x8
0 100
1 1000
2 1000
3 10
4 10
5 10
6 10
7 10
r
1 83333.3333333333
1 0
1 -1250000
1 400
1 300
1 100
C0
o16
o2
v0
v5
C1
o1
o2
v1
v3
o2
v1
v6
C2
o1
o2
v2
v4
o2
v2
v7
C3
n0
C4
n0
C5
n0
O0 0
n0
k7
1
2
3
7
11
13
15
J0 3
0 1
3 833.333333333333
5 0
J1 4
1 0
3 -1250
4 1250
6 0
J2 3
2 0
4 -2500
7 0
J3 2
3 1
5 1
J4 3
3 -1
4 1
6 1
J5 2
4 -1
7 1
G0 3
0 1
1 1
2 1
