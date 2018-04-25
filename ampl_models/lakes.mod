#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem:
#   ********
#   A problem of water resource management in Canada, which may be 
#   formulated as
#
#   Min  SUM   SUM  (T(i,j)- R(i,j))^2 + (O(i,j)-R(N+i,j))^2)
#       i=1,N j=1,5 
#
#   subject to
#
#   T(i+1,1)-T(i,1)+O(i,1)        =  G(i,1)
#   T(i+1,2)-T(i,2)-O(i,1)+O(i,2) =  G(i,2)
#   T(i+1,3)-T(i,3)-O(i,2)+O(i,3) =  G(i,3)
#   T(i+1,4)-T(i,4)-O(i,3)+O(i,4) =  G(i,4)
#   T(i+1,5)-T(i,5)-O(i,4)+O(i,5) =  G(i,5) 
#
#   i=1,N and T(N+1,j) = T(1,j)  for j=1,5
#
#   O(i,2)-a*((T(i,2)/480.8+T(i,3)/4.6)/2-543.4)^2 * 
#   (T(i,2)/480.8-T(i,3)/4.6)^.5=0
#   O(i,3)-b*((T(i,3)/4.6-543.4)^2*(T(i,3)/4.6-T(i,4)/105.15)^0.5) = 0
#   O(i,4)-c*(T(i,4)/105.15-550.11)^2.2 = 0
#
#   where T(i,j) and O(i,j) are variables, R(i,j) are given and
#   a=.0841168  b=.1280849 and c=0.2605.
#
#   Extra variables 
#   
#   v(i,2) = T(i,2) / 961.6 + T(i,3) / 9.2 - 543.4
#   w(i,2) = T(i,2) / 480.8 - T(i,3) / 4.6
#   v(i,3) = T(i,3) / 4.6 - 543.4
#   w(i,3) = T(i,3) / 4.6 - T(i,4) / 105.15
#   v(i,4) = T(i,4) / 105.15 - 550.11
#
#   are introduced so that the nonlinear constraints may be rewritten as
#
#   O(i,2)-a*v(i,2)^2 * w(i,2)^0.5 = 0 ; w(i,2) > 0
#   O(i,3)-b*v(i,3)^2 * w(i,3)^0.5 = 0 ; w(i,3) > 0
#   O(i,4)-c*v(i,4)^2.2 = 0 ; v(i,4) > 0
#
#   Source:
#   S Jafar Sadjadi
#   Dept. of Systems Design Engineering
#   University of Waterloo
#   Ontario, N2L 3G1 Canada
#   SIF input: Nick Gould and Jafar Sadjadi, November 1995
#   classification QOR2-RN-90-78
# Target value S_hat and O_hat
# Defining the variables
#  Objective function groups
# Linear Constraints
# Nonlinear constraints
# Artificial linear constraints
	param n := 6;
	param nm1 := -1 + (6);
	param np1 := 1 + (6);
	param nn := 2 * (6);
	param s1s1 := 202761.072;
	param s1s2 := 277791.816;
	param s1s3 := 2636.996;
	param s1s4 := 59987.0235;
	param s1s5 := 19490.4;
	param s2s1 := 202703.646;
	param s2s2 := 277849.512;
	param s2s3 := 2638.1;
	param s2s4 := 59998.59;
	param s2s5 := 19555.2;
	param s3s1 := 202720.536;
	param s3s2 := 277955.288;
	param s3s3 := 2639.894;
	param s3s4 := 60046.959;
	param s3s5 := 19597.6;
	param s4s1 := 202808.364;
	param s4s2 := 278104.336;
	param s4s3 := 2640.906;
	param s4s4 := 60074.298;
	param s4s5 := 19652.8;
	param s5s1 := 202916.46;
	param s5s2 := 278224.536;
	param s5s3 := 2641.458;
	param s5s4 := 60091.122;
	param s5s5 := 19708.8;
	param s6s1 := 202953.618;
	param s6s2 := 278277.424;
	param s6s3 := 2641.458;
	param s6s4 := 60082.71;
	param s6s5 := 19706.4;
	param o1o1 := 83.728;
	param o1o2 := 174.665;
	param o1o3 := 180.539;
	param o1o4 := 211.558;
	param o1o5 := 232.252;
	param o2o1 := 83.789;
	param o2o2 := 173.255;
	param o2o3 := 179.917;
	param o2o4 := 210.585;
	param o2o5 := 215.254;
	param o3o1 := 82.916;
	param o3o2 := 173.721;
	param o3o3 := 182.676;
	param o3o4 := 207.838;
	param o3o5 := 203.855;
	param o4o1 := 80.134;
	param o4o2 := 178.654;
	param o4o3 := 185.917;
	param o4o4 := 206.416;
	param o4o5 := 186.308;
	param o5o1 := 65.345;
	param o5o2 := 188.01;
	param o5o3 := 192.568;
	param o5o4 := 204.3;
	param o5o5 := 201.1;
	param o6o1 := 72.005;
	param o6o2 := 193.833;
	param o6o3 := 196.651;
	param o6o4 := 204.25;
	param o6o5 := 241.079;
	param npi := (6) + (6);
	param kp1 := 1 + (5);
	param jm1 := -1 + (5);
	param c := 1.0 / (105.15);

	var t1_1 := 1.0;
	var o1_1 := 1.0;
	var t1_2 := 1.0;
	var o1_2 := 1.0;
	var t1_3 := 1.0;
	var o1_3 := 1.0;
	var t1_4 := 1.0;
	var o1_4 := 1.0;
	var t1_5 := 1.0;
	var o1_5 := 1.0;
	var v1_2 := 1.0;
	var w1_2 >= 0.0001 ,  := 1.0;
	var v1_3 := 1.0;
	var w1_3 >= 0.0001 ,  := 1.0;
	var v1_4 >= 0.0001 ,  := 1.0;
	var t2_1 := 1.0;
	var o2_1 := 1.0;
	var t2_2 := 1.0;
	var o2_2 := 1.0;
	var t2_3 := 1.0;
	var o2_3 := 1.0;
	var t2_4 := 1.0;
	var o2_4 := 1.0;
	var t2_5 := 1.0;
	var o2_5 := 1.0;
	var v2_2 := 1.0;
	var w2_2 >= 0.0001 ,  := 1.0;
	var v2_3 := 1.0;
	var w2_3 >= 0.0001 ,  := 1.0;
	var v2_4 >= 0.0001 ,  := 1.0;
	var t3_1 := 1.0;
	var o3_1 := 1.0;
	var t3_2 := 1.0;
	var o3_2 := 1.0;
	var t3_3 := 1.0;
	var o3_3 := 1.0;
	var t3_4 := 1.0;
	var o3_4 := 1.0;
	var t3_5 := 1.0;
	var o3_5 := 1.0;
	var v3_2 := 1.0;
	var w3_2 >= 0.0001 ,  := 1.0;
	var v3_3 := 1.0;
	var w3_3 >= 0.0001 ,  := 1.0;
	var v3_4 >= 0.0001 ,  := 1.0;
	var t4_1 := 1.0;
	var o4_1 := 1.0;
	var t4_2 := 1.0;
	var o4_2 := 1.0;
	var t4_3 := 1.0;
	var o4_3 := 1.0;
	var t4_4 := 1.0;
	var o4_4 := 1.0;
	var t4_5 := 1.0;
	var o4_5 := 1.0;
	var v4_2 := 1.0;
	var w4_2 >= 0.0001 ,  := 1.0;
	var v4_3 := 1.0;
	var w4_3 >= 0.0001 ,  := 1.0;
	var v4_4 >= 0.0001 ,  := 1.0;
	var t5_1 := 1.0;
	var o5_1 := 1.0;
	var t5_2 := 1.0;
	var o5_2 := 1.0;
	var t5_3 := 1.0;
	var o5_3 := 1.0;
	var t5_4 := 1.0;
	var o5_4 := 1.0;
	var t5_5 := 1.0;
	var o5_5 := 1.0;
	var v5_2 := 1.0;
	var w5_2 >= 0.0001 ,  := 1.0;
	var v5_3 := 1.0;
	var w5_3 >= 0.0001 ,  := 1.0;
	var v5_4 >= 0.0001 ,  := 1.0;
	var t6_1 := 1.0;
	var o6_1 := 1.0;
	var t6_2 := 1.0;
	var o6_2 := 1.0;
	var t6_3 := 1.0;
	var o6_3 := 1.0;
	var t6_4 := 1.0;
	var o6_4 := 1.0;
	var t6_5 := 1.0;
	var o6_5 := 1.0;
	var v6_2 := 1.0;
	var w6_2 >= 0.0001 ,  := 1.0;
	var v6_3 := 1.0;
	var w6_3 >= 0.0001 ,  := 1.0;
	var v6_4 >= 0.0001 ,  := 1.0;

minimize obj:
	(t1_1 - 202761.072)*(t1_1 - 202761.072) + (o1_1 - 232.252)*(o1_1 - 232.252) + 
	(t1_2 - 277791.816)*(t1_2 - 277791.816) + (o1_2 - 174.665)*(o1_2 - 174.665) + 
	(t1_3 - 2636.996)*(t1_3 - 2636.996) + (o1_3 - 180.539)*(o1_3 - 180.539) + (t1_4 
	- 59987.0235)*(t1_4 - 59987.0235) + (o1_4 - 211.558)*(o1_4 - 211.558) + (t1_5 - 
	19490.4)*(t1_5 - 19490.4) + (o1_5 - 232.252)*(o1_5 - 232.252) + (t2_1 - 
	202703.646)*(t2_1 - 202703.646) + (o2_1 - 83.789)*(o2_1 - 83.789) + (t2_2 - 
	277849.512)*(t2_2 - 277849.512) + (o2_2 - 173.255)*(o2_2 - 173.255) + (t2_3 - 
	2638.1)*(t2_3 - 2638.1) + (o2_3 - 179.917)*(o2_3 - 179.917) + (t2_4 - 
	59998.59)*(t2_4 - 59998.59) + (o2_4 - 210.585)*(o2_4 - 210.585) + (t2_5 - 
	19555.2)*(t2_5 - 19555.2) + (o2_5 - 215.254)*(o2_5 - 215.254) + (t3_1 - 
	202720.536)*(t3_1 - 202720.536) + (o3_1 - 82.916)*(o3_1 - 82.916) + (t3_2 - 
	277955.288)*(t3_2 - 277955.288) + (o3_2 - 173.721)*(o3_2 - 173.721) + (t3_3 - 
	2639.894)*(t3_3 - 2639.894) + (o3_3 - 182.676)*(o3_3 - 182.676) + (t3_4 - 
	60046.959)*(t3_4 - 60046.959) + (o3_4 - 207.838)*(o3_4 - 207.838) + (t3_5 - 
	19597.6)*(t3_5 - 19597.6) + (o3_5 - 203.855)*(o3_5 - 203.855) + (t4_1 - 
	202808.364)*(t4_1 - 202808.364) + (o4_1 - 80.134)*(o4_1 - 80.134) + (t4_2 - 
	278104.336)*(t4_2 - 278104.336) + (o4_2 - 178.654)*(o4_2 - 178.654) + (t4_3 - 
	2640.906)*(t4_3 - 2640.906) + (o4_3 - 185.917)*(o4_3 - 185.917) + (t4_4 - 
	60074.298)*(t4_4 - 60074.298) + (o4_4 - 206.416)*(o4_4 - 206.416) + (t4_5 - 
	19652.8)*(t4_5 - 19652.8) + (o4_5 - 186.308)*(o4_5 - 186.308) + (t5_1 - 
	202916.46)*(t5_1 - 202916.46) + (o5_1 - 65.345)*(o5_1 - 65.345) + (t5_2 - 
	278224.536)*(t5_2 - 278224.536) + (o5_2 - 188.01)*(o5_2 - 188.01) + (t5_3 - 
	2641.458)*(t5_3 - 2641.458) + (o5_3 - 192.568)*(o5_3 - 192.568) + (t5_4 - 
	60091.122)*(t5_4 - 60091.122) + (o5_4 - 204.3)*(o5_4 - 204.3) + (t5_5 - 
	19708.8)*(t5_5 - 19708.8) + (o5_5 - 201.1)*(o5_5 - 201.1) + (t6_1 - 
	202953.618)*(t6_1 - 202953.618) + (o6_1 - 72.005)*(o6_1 - 72.005) + (t6_2 - 
	278277.424)*(t6_2 - 278277.424) + (o6_2 - 193.833)*(o6_2 - 193.833) + (t6_3 - 
	2641.458)*(t6_3 - 2641.458) + (o6_3 - 196.651)*(o6_3 - 196.651) + (t6_4 - 
	60082.71)*(t6_4 - 60082.71) + (o6_4 - 204.25)*(o6_4 - 204.25) + (t6_5 - 
	19706.4)*(t6_5 - 19706.4) + (o6_5 - 241.079)*(o6_5 - 241.079);

subject to consg1_1 := 1.0:
	t2_1 - t1_1 + o1_1 + 22.0 = 0;
subject to consg1_2 := 1.0:
	t2_2 - t1_2 + o1_2 - o1_1 + 1.0 = 0;
subject to consg1_3 := 1.0:
	t2_3 - t1_3 + o1_3 - o1_2 - 3.0 = 0;
subject to consg1_4 := 1.0:
	t2_4 - t1_4 + o1_4 - o1_3 + 27.2 = 0;
subject to consg1_5 := 1.0:
	t2_5 - t1_5 + o1_5 - o1_4 - 51.5 = 0;
subject to consg2_1 := 1.0:
	t3_1 - t2_1 + o2_1 - 44.0 = 0;
subject to consg2_2 := 1.0:
	t3_2 - t2_2 + o2_2 - o2_1 - 162.0 = 0;
subject to consg2_3 := 1.0:
	t3_3 - t2_3 + o2_3 - o2_2 - 8.0 = 0;
subject to consg2_4 := 1.0:
	t3_4 - t2_4 + o2_4 - o2_3 - 12.5 = 0;
subject to consg2_5 := 1.0:
	t3_5 - t2_5 + o2_5 - o2_4 - 53.5 = 0;
subject to consg3_1 := 1.0:
	t4_1 - t3_1 + o3_1 + 11.0 = 0;
subject to consg3_2 := 1.0:
	t4_2 - t3_2 + o3_2 - o3_1 - 60.0 = 0;
subject to consg3_3 := 1.0:
	t4_3 - t3_3 + o3_3 - o3_2 - 10.0 = 0;
subject to consg3_4 := 1.0:
	t4_4 - t3_4 + o3_4 - o3_3 - 18.0 = 0;
subject to consg3_5 := 1.0:
	t4_5 - t3_5 + o3_5 - o3_4 - 39.0 = 0;
subject to consg4_1 := 1.0:
	t5_1 - t4_1 + o4_1 - 124.0 = 0;
subject to consg4_2 := 1.0:
	t5_2 - t4_2 + o4_2 - o4_1 - 246.0 = 0;
subject to consg4_3 := 1.0:
	t5_3 - t4_3 + o4_3 - o4_2 - 6.0 = 0;
subject to consg4_4 := 1.0:
	t5_4 - t4_4 + o4_4 - o4_3 - 9.7 = 0;
subject to consg4_5 := 1.0:
	t5_5 - t4_5 + o4_5 - o4_4 - 17.2 = 0;
subject to consg5_1 := 1.0:
	t6_1 - t5_1 + o5_1 - 127.0 = 0;
subject to consg5_2 := 1.0:
	t6_2 - t5_2 + o5_2 - o5_1 - 175.0 = 0;
subject to consg5_3 := 1.0:
	t6_3 - t5_3 + o5_3 - o5_2 - 3.0 = 0;
subject to consg5_4 := 1.0:
	t6_4 - t5_4 + o5_4 - o5_3 - 10.0 = 0;
subject to consg5_5 := 1.0:
	t6_5 - t5_5 + o5_5 - o5_4 - 30.2 = 0;
subject to consg6_1 := 1.0:
	t1_1 - t6_1 + o6_1 - 78.0 = 0;
subject to consg6_2 := 1.0:
	t1_2 - t6_2 + o6_2 - o6_1 - 156.0 = 0;
subject to consg6_3 := 1.0:
	t1_3 - t6_3 + o6_3 - o6_2 - 3.0 = 0;
subject to consg6_4 := 1.0:
	t1_4 - t6_4 + o6_4 - o6_3 - 14.0 = 0;
subject to consg6_5 := 1.0:
	t1_5 - t6_5 + o6_5 - o6_4 - 23.2 = 0;
subject to consa1_1 := 1.0:
	0.0841168 * v1_2^2 * w1_2^0.5 - o1_2 = 0;
subject to consa1_2 := 1.0:
	0.1280849 * v1_3^2 * w1_3^0.5 - o1_3 = 0;
subject to consa1_3 := 1.0:
	0.2605 * v1_4^2.2 - o1_4 = 0;
subject to consv1_2 := 1.0:
	-v1_2 + 0.0010399334442595673*t1_2 + 0.10869565217391305*t1_3 - 543.4 = 0;
subject to consw1_2 := 1.0:
	-w1_2 + 0.0020798668885191347*t1_2 - 0.2173913043478261*t1_3 = 0;
subject to consv1_3 := 1.0:
	-v1_3 + 0.2173913043478261*t1_3 - 543.4 = 0;
subject to consw1_3 := 1.0:
	-w1_3 + 0.2173913043478261*t1_3 - 0.009510223490252021*t1_4 = 0;
subject to consv1_4 := 1.0:
	-v1_4 + 0.009510223490252021*t1_4 - 550.11 = 0;
subject to consa2_1 := 1.0:
	0.0841168 * v2_2^2 * w2_2^0.5 - o2_2 = 0;
subject to consa2_2 := 1.0:
	0.1280849 * v2_3^2 * w2_3^0.5 - o2_3 = 0;
subject to consa2_3 := 1.0:
	0.2605 * v2_4^2.2 - o2_4 = 0;
subject to consv2_2 := 1.0:
	-v2_2 + 0.0010399334442595673*t2_2 + 0.10869565217391305*t2_3 - 543.4 = 0;
subject to consw2_2 := 1.0:
	-w2_2 + 0.0020798668885191347*t2_2 - 0.2173913043478261*t2_3 = 0;
subject to consv2_3 := 1.0:
	-v2_3 + 0.2173913043478261*t2_3 - 543.4 = 0;
subject to consw2_3 := 1.0:
	-w2_3 + 0.2173913043478261*t2_3 - 0.009510223490252021*t2_4 = 0;
subject to consv2_4 := 1.0:
	-v2_4 + 0.009510223490252021*t2_4 - 550.11 = 0;
subject to consa3_1 := 1.0:
	0.0841168 * v3_2^2 * w3_2^0.5 - o3_2 = 0;
subject to consa3_2 := 1.0:
	0.1280849 * v3_3^2 * w3_3^0.5 - o3_3 = 0;
subject to consa3_3 := 1.0:
	0.2605 * v3_4^2.2 - o3_4 = 0;
subject to consv3_2 := 1.0:
	-v3_2 + 0.0010399334442595673*t3_2 + 0.10869565217391305*t3_3 - 543.4 = 0;
subject to consw3_2 := 1.0:
	-w3_2 + 0.0020798668885191347*t3_2 - 0.2173913043478261*t3_3 = 0;
subject to consv3_3 := 1.0:
	-v3_3 + 0.2173913043478261*t3_3 - 543.4 = 0;
subject to consw3_3 := 1.0:
	-w3_3 + 0.2173913043478261*t3_3 - 0.009510223490252021*t3_4 = 0;
subject to consv3_4 := 1.0:
	-v3_4 + 0.009510223490252021*t3_4 - 550.11 = 0;
subject to consa4_1 := 1.0:
	0.0841168 * v4_2^2 * w4_2^0.5 - o4_2 = 0;
subject to consa4_2 := 1.0:
	0.1280849 * v4_3^2 * w4_3^0.5 - o4_3 = 0;
subject to consa4_3 := 1.0:
	0.2605 * v4_4^2.2 - o4_4 = 0;
subject to consv4_2 := 1.0:
	-v4_2 + 0.0010399334442595673*t4_2 + 0.10869565217391305*t4_3 - 543.4 = 0;
subject to consw4_2 := 1.0:
	-w4_2 + 0.0020798668885191347*t4_2 - 0.2173913043478261*t4_3 = 0;
subject to consv4_3 := 1.0:
	-v4_3 + 0.2173913043478261*t4_3 - 543.4 = 0;
subject to consw4_3 := 1.0:
	-w4_3 + 0.2173913043478261*t4_3 - 0.009510223490252021*t4_4 = 0;
subject to consv4_4 := 1.0:
	-v4_4 + 0.009510223490252021*t4_4 - 550.11 = 0;
subject to consa5_1 := 1.0:
	0.0841168 * v5_2^2 * w5_2^0.5 - o5_2 = 0;
subject to consa5_2 := 1.0:
	0.1280849 * v5_3^2 * w5_3^0.5 - o5_3 = 0;
subject to consa5_3 := 1.0:
	0.2605 * v5_4^2.2 - o5_4 = 0;
subject to consv5_2 := 1.0:
	-v5_2 + 0.0010399334442595673*t5_2 + 0.10869565217391305*t5_3 - 543.4 = 0;
subject to consw5_2 := 1.0:
	-w5_2 + 0.0020798668885191347*t5_2 - 0.2173913043478261*t5_3 = 0;
subject to consv5_3 := 1.0:
	-v5_3 + 0.2173913043478261*t5_3 - 543.4 = 0;
subject to consw5_3 := 1.0:
	-w5_3 + 0.2173913043478261*t5_3 - 0.009510223490252021*t5_4 = 0;
subject to consv5_4 := 1.0:
	-v5_4 + 0.009510223490252021*t5_4 - 550.11 = 0;
subject to consa6_1 := 1.0:
	0.0841168 * v6_2^2 * w6_2^0.5 - o6_2 = 0;
subject to consa6_2 := 1.0:
	0.1280849 * v6_3^2 * w6_3^0.5 - o6_3 = 0;
subject to consa6_3 := 1.0:
	0.2605 * v6_4^2.2 - o6_4 = 0;
subject to consv6_2 := 1.0:
	-v6_2 + 0.0010399334442595673*t6_2 + 0.10869565217391305*t6_3 - 543.4 = 0;
subject to consw6_2 := 1.0:
	-w6_2 + 0.0020798668885191347*t6_2 - 0.2173913043478261*t6_3 = 0;
subject to consv6_3 := 1.0:
	-v6_3 + 0.2173913043478261*t6_3 - 543.4 = 0;
subject to consw6_3 := 1.0:
	-w6_3 + 0.2173913043478261*t6_3 - 0.009510223490252021*t6_4 = 0;
subject to consv6_4 := 1.0:
	-v6_4 + 0.009510223490252021*t6_4 - 550.11 = 0;

solve;
	display t1_1;
	display o1_1;
	display t1_2;
	display o1_2;
	display t1_3;
	display o1_3;
	display t1_4;
	display o1_4;
	display t1_5;
	display o1_5;
	display v1_2;
	display w1_2;
	display v1_3;
	display w1_3;
	display v1_4;
	display t2_1;
	display o2_1;
	display t2_2;
	display o2_2;
	display t2_3;
	display o2_3;
	display t2_4;
	display o2_4;
	display t2_5;
	display o2_5;
	display v2_2;
	display w2_2;
	display v2_3;
	display w2_3;
	display v2_4;
	display t3_1;
	display o3_1;
	display t3_2;
	display o3_2;
	display t3_3;
	display o3_3;
	display t3_4;
	display o3_4;
	display t3_5;
	display o3_5;
	display v3_2;
	display w3_2;
	display v3_3;
	display w3_3;
	display v3_4;
	display t4_1;
	display o4_1;
	display t4_2;
	display o4_2;
	display t4_3;
	display o4_3;
	display t4_4;
	display o4_4;
	display t4_5;
	display o4_5;
	display v4_2;
	display w4_2;
	display v4_3;
	display w4_3;
	display v4_4;
	display t5_1;
	display o5_1;
	display t5_2;
	display o5_2;
	display t5_3;
	display o5_3;
	display t5_4;
	display o5_4;
	display t5_5;
	display o5_5;
	display v5_2;
	display w5_2;
	display v5_3;
	display w5_3;
	display v5_4;
	display t6_1;
	display o6_1;
	display t6_2;
	display o6_2;
	display t6_3;
	display o6_3;
	display t6_4;
	display o6_4;
	display t6_5;
	display o6_5;
	display v6_2;
	display w6_2;
	display v6_3;
	display w6_3;
	display v6_4;
display obj;
