#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   The Power Generation problem for the SSGB.
#   Source:
#   N. Gould, private communication.
#   SIF input: Nick Gould, 23 October 1989
#   classification LQR2-RN-194-96
#   period is the number of time periods
#   Define a few helpful parameters
#   set the influx z (MW) and demands d (MW)
#   day 1
#   day 2
#   day 3
#   day 4
#   day 5
#   day 6
#   day 7
#   Define the problem variables
#  power output oil (MW)
#  power output coal (MW)
#  power generated  (MW)
#  spillage  (MW)
#  Flow per unit time from pump storage - generating (MW)
#  Flow per unit time from pump storage - pumping (MW)
#  Volume of hydro power (MWh)
#  Volume of pump-storage power (MWh)
#   Define objective function group, cost of generation (MW)
#   Define the hydro-power constraints
#   Define the power-storage reservoir constraints
#   Define the demand constraints
# Nonlinear constraints
#   The constant for the hydro-power constraints are nonzero.
#   The constant for the reservoir constraints are nonzero.
#   The constant for all the demand constraints are nonzero.
# Nonlinear constraints
#   Solution
param hours := 24;
param days := 1;
param z1 := 517.0;
param d1_1 := 578.0;
param d1_2 := 517.0;
param d1_3 := 461.0;
param d1_4 := 369.0;
param d1_5 := 299.0;
param d1_6 := 269.0;
param d1_7 := 370.0;
param d1_8 := 559.0;
param d1_9 := 689.0;
param d1_10 := 728.0;
param d1_11 := 683.0;
param d1_12 := 626.0;
param d1_13 := 619.0;
param d1_14 := 586.0;
param d1_15 := 582.0;
param d1_16 := 625.0;
param d1_17 := 821.0;
param d1_18 := 883.0;
param d1_19 := 768.0;
param d1_20 := 711.0;
param d1_21 := 677.0;
param d1_22 := 630.0;
param d1_23 := 545.0;
param d1_24 := 565.0;
param z2 := 400.0;
param d2_1 := 631.0;
param d2_2 := 574.0;
param d2_3 := 521.0;
param d2_4 := 446.0;
param d2_5 := 359.0;
param d2_6 := 336.0;
param d2_7 := 420.0;
param d2_8 := 588.0;
param d2_9 := 697.0;
param d2_10 := 732.0;
param d2_11 := 713.0;
param d2_12 := 682.0;
param d2_13 := 695.0;
param d2_14 := 651.0;
param d2_15 := 645.0;
param d2_16 := 664.0;
param d2_17 := 816.0;
param d2_18 := 858.0;
param d2_19 := 760.0;
param d2_20 := 700.0;
param d2_21 := 659.0;
param d2_22 := 623.0;
param d2_23 := 517.0;
param d2_24 := 542.0;
param z3 := 1017.0;
param d3_1 := 582.0;
param d3_2 := 501.0;
param d3_3 := 443.0;
param d3_4 := 367.0;
param d3_5 := 288.0;
param d3_6 := 265.0;
param d3_7 := 349.0;
param d3_8 := 503.0;
param d3_9 := 663.0;
param d3_10 := 651.0;
param d3_11 := 625.0;
param d3_12 := 596.0;
param d3_13 := 608.0;
param d3_14 := 566.0;
param d3_15 := 555.0;
param d3_16 := 584.0;
param d3_17 := 763.0;
param d3_18 := 803.0;
param d3_19 := 710.0;
param d3_20 := 648.0;
param d3_21 := 626.0;
param d3_22 := 590.0;
param d3_23 := 486.0;
param d3_24 := 540.0;
param z4 := 667.0;
param d4_1 := 602.0;
param d4_2 := 533.0;
param d4_3 := 450.0;
param d4_4 := 378.0;
param d4_5 := 298.0;
param d4_6 := 272.0;
param d4_7 := 369.0;
param d4_8 := 539.0;
param d4_9 := 647.0;
param d4_10 := 652.0;
param d4_11 := 607.0;
param d4_12 := 585.0;
param d4_13 := 587.0;
param d4_14 := 549.0;
param d4_15 := 535.0;
param d4_16 := 564.0;
param d4_17 := 748.0;
param d4_18 := 808.0;
param d4_19 := 710.0;
param d4_20 := 646.0;
param d4_21 := 620.0;
param d4_22 := 581.0;
param d4_23 := 483.0;
param d4_24 := 514.0;
param z5 := 600.0;
param d5_1 := 579.0;
param d5_2 := 518.0;
param d5_3 := 447.0;
param d5_4 := 355.0;
param d5_5 := 284.0;
param d5_6 := 261.0;
param d5_7 := 348.0;
param d5_8 := 530.0;
param d5_9 := 644.0;
param d5_10 := 648.0;
param d5_11 := 607.0;
param d5_12 := 570.0;
param d5_13 := 577.0;
param d5_14 := 536.0;
param d5_15 := 544.0;
param d5_16 := 554.0;
param d5_17 := 716.0;
param d5_18 := 765.0;
param d5_19 := 676.0;
param d5_20 := 631.0;
param d5_21 := 576.0;
param d5_22 := 528.0;
param d5_23 := 445.0;
param d5_24 := 520.0;
param z6 := 421.0;
param d6_1 := 618.0;
param d6_2 := 547.0;
param d6_3 := 430.0;
param d6_4 := 327.0;
param d6_5 := 249.0;
param d6_6 := 211.0;
param d6_7 := 227.0;
param d6_8 := 258.0;
param d6_9 := 347.0;
param d6_10 := 491.0;
param d6_11 := 524.0;
param d6_12 := 492.0;
param d6_13 := 467.0;
param d6_14 := 418.0;
param d6_15 := 358.0;
param d6_16 := 378.0;
param d6_17 := 544.0;
param d6_18 := 666.0;
param d6_19 := 589.0;
param d6_20 := 533.0;
param d6_21 := 494.0;
param d6_22 := 460.0;
param d6_23 := 404.0;
param d6_24 := 512.0;
param z7 := 425.0;
param d7_1 := 615.0;
param d7_2 := 587.0;
param d7_3 := 450.0;
param d7_4 := 320.0;
param d7_5 := 235.0;
param d7_6 := 198.0;
param d7_7 := 195.0;
param d7_8 := 173.0;
param d7_9 := 197.0;
param d7_10 := 349.0;
param d7_11 := 441.0;
param d7_12 := 459.0;
param d7_13 := 485.0;
param d7_14 := 445.0;
param d7_15 := 410.0;
param d7_16 := 421.0;
param d7_17 := 568.0;
param d7_18 := 643.0;
param d7_19 := 596.0;
param d7_20 := 566.0;
param d7_21 := 541.0;
param d7_22 := 532.0;
param d7_23 := 454.0;
param d7_24 := 511.0;
param p := -1 + (1);
param ihm1 := -1 + (24);

var v0_24 >= 240000.0 ,  <= 240000.0 ,  := 240000.0;
var r0_24 >= 3500.0 ,  <= 3500.0 ,  := 3500.0;
var p11_1 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_1 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_1 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_1 >= 0.0;
var qg1_1 >= 0.0 ,  <= 300.0;
var qp1_1 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_1 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_1 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_2 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_2 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_2 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_2 >= 0.0;
var qg1_2 >= 0.0 ,  <= 300.0;
var qp1_2 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_2 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_2 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_3 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_3 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_3 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_3 >= 0.0;
var qg1_3 >= 0.0 ,  <= 300.0;
var qp1_3 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_3 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_3 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_4 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_4 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_4 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_4 >= 0.0;
var qg1_4 >= 0.0 ,  <= 300.0;
var qp1_4 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_4 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_4 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_5 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_5 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_5 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_5 >= 0.0;
var qg1_5 >= 0.0 ,  <= 300.0;
var qp1_5 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_5 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_5 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_6 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_6 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_6 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_6 >= 0.0;
var qg1_6 >= 0.0 ,  <= 300.0;
var qp1_6 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_6 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_6 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_7 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_7 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_7 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_7 >= 0.0;
var qg1_7 >= 0.0 ,  <= 300.0;
var qp1_7 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_7 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_7 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_8 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_8 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_8 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_8 >= 0.0;
var qg1_8 >= 0.0 ,  <= 300.0;
var qp1_8 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_8 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_8 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_9 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_9 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_9 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_9 >= 0.0;
var qg1_9 >= 0.0 ,  <= 300.0;
var qp1_9 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_9 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_9 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_10 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_10 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_10 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_10 >= 0.0;
var qg1_10 >= 0.0 ,  <= 300.0;
var qp1_10 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_10 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_10 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_11 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_11 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_11 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_11 >= 0.0;
var qg1_11 >= 0.0 ,  <= 300.0;
var qp1_11 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_11 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_11 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_12 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_12 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_12 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_12 >= 0.0;
var qg1_12 >= 0.0 ,  <= 300.0;
var qp1_12 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_12 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_12 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_13 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_13 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_13 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_13 >= 0.0;
var qg1_13 >= 0.0 ,  <= 300.0;
var qp1_13 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_13 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_13 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_14 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_14 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_14 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_14 >= 0.0;
var qg1_14 >= 0.0 ,  <= 300.0;
var qp1_14 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_14 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_14 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_15 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_15 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_15 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_15 >= 0.0;
var qg1_15 >= 0.0 ,  <= 300.0;
var qp1_15 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_15 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_15 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_16 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_16 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_16 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_16 >= 0.0;
var qg1_16 >= 0.0 ,  <= 300.0;
var qp1_16 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_16 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_16 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_17 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_17 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_17 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_17 >= 0.0;
var qg1_17 >= 0.0 ,  <= 300.0;
var qp1_17 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_17 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_17 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_18 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_18 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_18 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_18 >= 0.0;
var qg1_18 >= 0.0 ,  <= 300.0;
var qp1_18 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_18 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_18 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_19 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_19 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_19 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_19 >= 0.0;
var qg1_19 >= 0.0 ,  <= 300.0;
var qp1_19 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_19 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_19 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_20 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_20 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_20 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_20 >= 0.0;
var qg1_20 >= 0.0 ,  <= 300.0;
var qp1_20 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_20 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_20 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_21 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_21 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_21 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_21 >= 0.0;
var qg1_21 >= 0.0 ,  <= 300.0;
var qp1_21 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_21 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_21 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_22 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_22 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_22 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_22 >= 0.0;
var qg1_22 >= 0.0 ,  <= 300.0;
var qp1_22 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_22 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_22 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_23 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_23 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_23 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_23 >= 0.0;
var qg1_23 >= 0.0 ,  <= 300.0;
var qp1_23 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_23 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_23 >= 0.0 ,  <= 6000.0 ,  := 3500;
var p11_24 >= 70.0 ,  <= 325.0 ,  := 70.0;
var p21_24 >= 90.0 ,  <= 290.0 ,  := 90.0;
var qh1_24 >= 25.0 ,  <= 500.0 ,  := 25.0;
var s1_24 >= 0.0;
var qg1_24 >= 0.0 ,  <= 300.0;
var qp1_24 >= 0.0 ,  <= 225.0 ,  := 225.0;
var v1_24 >= 180000.0 ,  <= 280000.0 ,  := 240000.0;
var r1_24 >= 0.0 ,  <= 6000.0 ,  := 3500;

minimize obj:
	1000.0*p11_1 + 1500.0*p21_1 + 1200.0*qh1_1 + 1200.0*s1_1 + 1200.0*qg1_1 - 1200.0*qp1_1 + 1000.0*p11_2 + 1500.0*p21_2 + 1200.0*qh1_2 + 1200.0*s1_2 + 1200.0*qg1_2 - 1200.0*qp1_2 + 1000.0*p11_3 + 1500.0*p21_3 + 1200.0*qh1_3 + 1200.0*s1_3 + 1200.0*qg1_3 - 1200.0*qp1_3 + 1000.0*p11_4 + 1500.0*p21_4 + 1200.0*qh1_4 + 1200.0*s1_4 + 1200.0*qg1_4 - 1200.0*qp1_4 + 1000.0*p11_5 + 1500.0*p21_5 + 1200.0*qh1_5 + 1200.0*s1_5 + 1200.0*qg1_5 - 1200.0*qp1_5 + 1000.0*p11_6 + 1500.0*p21_6 + 1200.0*qh1_6 + 1200.0*s1_6 + 1200.0*qg1_6 - 1200.0*qp1_6 + 1000.0*p11_7 + 1500.0*p21_7 + 1200.0*qh1_7 + 1200.0*s1_7 + 1200.0*qg1_7 - 1200.0*qp1_7 + 1000.0*p11_8 + 1500.0*p21_8 + 1200.0*qh1_8 + 1200.0*s1_8 + 1200.0*qg1_8 - 1200.0*qp1_8 + 1000.0*p11_9 + 1500.0*p21_9 + 1200.0*qh1_9 + 1200.0*s1_9 + 1200.0*qg1_9 - 1200.0*qp1_9 + 1000.0*p11_10 + 1500.0*p21_10 + 1200.0*qh1_10 + 1200.0*s1_10 + 1200.0*qg1_10 - 1200.0*qp1_10 + 1000.0*p11_11 + 1500.0*p21_11 + 1200.0*qh1_11 + 1200.0*s1_11 + 1200.0*qg1_11 - 1200.0*qp1_11 + 1000.0*p11_12 + 1500.0*p21_12 + 1200.0*qh1_12 + 1200.0*s1_12 + 1200.0*qg1_12 - 1200.0*qp1_12 + 1000.0*p11_13 + 1500.0*p21_13 + 1200.0*qh1_13 + 1200.0*s1_13 + 1200.0*qg1_13 - 1200.0*qp1_13 + 1000.0*p11_14 + 1500.0*p21_14 + 1200.0*qh1_14 + 1200.0*s1_14 + 1200.0*qg1_14 - 1200.0*qp1_14 + 1000.0*p11_15 + 1500.0*p21_15 + 1200.0*qh1_15 + 1200.0*s1_15 + 1200.0*qg1_15 - 1200.0*qp1_15 + 1000.0*p11_16 + 1500.0*p21_16 + 1200.0*qh1_16 + 1200.0*s1_16 + 1200.0*qg1_16 - 1200.0*qp1_16 + 1000.0*p11_17 + 1500.0*p21_17 + 1200.0*qh1_17 + 1200.0*s1_17 + 1200.0*qg1_17 - 1200.0*qp1_17 + 1000.0*p11_18 + 1500.0*p21_18 + 1200.0*qh1_18 + 1200.0*s1_18 + 1200.0*qg1_18 - 1200.0*qp1_18 + 1000.0*p11_19 + 1500.0*p21_19 + 1200.0*qh1_19 + 1200.0*s1_19 + 1200.0*qg1_19 - 1200.0*qp1_19 + 1000.0*p11_20 + 1500.0*p21_20 + 1200.0*qh1_20 + 1200.0*s1_20 + 1200.0*qg1_20 - 1200.0*qp1_20 + 1000.0*p11_21 + 1500.0*p21_21 + 1200.0*qh1_21 + 1200.0*s1_21 + 1200.0*qg1_21 - 1200.0*qp1_21 + 1000.0*p11_22 + 1500.0*p21_22 + 1200.0*qh1_22 + 1200.0*s1_22 + 1200.0*qg1_22 - 1200.0*qp1_22 + 1000.0*p11_23 + 1500.0*p21_23 + 1200.0*qh1_23 + 1200.0*s1_23 + 1200.0*qg1_23 - 1200.0*qp1_23 + 1000.0*p11_24 + 1500.0*p21_24 + 1200.0*qh1_24 + 1200.0*s1_24 + 1200.0*qg1_24 - 1200.0*qp1_24;

subject to cons_h1_1:
	v1_1 - v0_24 + s1_1 + qh1_1 - 517.0 = 0;
subject to cons_h1_2:
	v1_2 - v1_1 + s1_2 + qh1_2 - 517.0 = 0;
subject to cons_h1_3:
	v1_3 - v1_2 + s1_3 + qh1_3 - 517.0 = 0;
subject to cons_h1_4:
	v1_4 - v1_3 + s1_4 + qh1_4 - 517.0 = 0;
subject to cons_h1_5:
	v1_5 - v1_4 + s1_5 + qh1_5 - 517.0 = 0;
subject to cons_h1_6:
	v1_6 - v1_5 + s1_6 + qh1_6 - 517.0 = 0;
subject to cons_h1_7:
	v1_7 - v1_6 + s1_7 + qh1_7 - 517.0 = 0;
subject to cons_h1_8:
	v1_8 - v1_7 + s1_8 + qh1_8 - 517.0 = 0;
subject to cons_h1_9:
	v1_9 - v1_8 + s1_9 + qh1_9 - 517.0 = 0;
subject to cons_h1_10:
	v1_10 - v1_9 + s1_10 + qh1_10 - 517.0 = 0;
subject to cons_h1_11:
	v1_11 - v1_10 + s1_11 + qh1_11 - 517.0 = 0;
subject to cons_h1_12:
	v1_12 - v1_11 + s1_12 + qh1_12 - 517.0 = 0;
subject to cons_h1_13:
	v1_13 - v1_12 + s1_13 + qh1_13 - 517.0 = 0;
subject to cons_h1_14:
	v1_14 - v1_13 + s1_14 + qh1_14 - 517.0 = 0;
subject to cons_h1_15:
	v1_15 - v1_14 + s1_15 + qh1_15 - 517.0 = 0;
subject to cons_h1_16:
	v1_16 - v1_15 + s1_16 + qh1_16 - 517.0 = 0;
subject to cons_h1_17:
	v1_17 - v1_16 + s1_17 + qh1_17 - 517.0 = 0;
subject to cons_h1_18:
	v1_18 - v1_17 + s1_18 + qh1_18 - 517.0 = 0;
subject to cons_h1_19:
	v1_19 - v1_18 + s1_19 + qh1_19 - 517.0 = 0;
subject to cons_h1_20:
	v1_20 - v1_19 + s1_20 + qh1_20 - 517.0 = 0;
subject to cons_h1_21:
	v1_21 - v1_20 + s1_21 + qh1_21 - 517.0 = 0;
subject to cons_h1_22:
	v1_22 - v1_21 + s1_22 + qh1_22 - 517.0 = 0;
subject to cons_h1_23:
	v1_23 - v1_22 + s1_23 + qh1_23 - 517.0 = 0;
subject to cons_h1_24:
	v1_24 - v1_23 + s1_24 + qh1_24 - 517.0 = 0;
subject to cons_r1_1:
	r1_1 - r0_24 + qg1_1 - qp1_1 - 5.17 = 0;
subject to cons_r1_2:
	r1_2 - r1_1 + qg1_2 - qp1_2 - 5.17 = 0;
subject to cons_r1_3:
	r1_3 - r1_2 + qg1_3 - qp1_3 - 5.17 = 0;
subject to cons_r1_4:
	r1_4 - r1_3 + qg1_4 - qp1_4 - 5.17 = 0;
subject to cons_r1_5:
	r1_5 - r1_4 + qg1_5 - qp1_5 - 5.17 = 0;
subject to cons_r1_6:
	r1_6 - r1_5 + qg1_6 - qp1_6 - 5.17 = 0;
subject to cons_r1_7:
	r1_7 - r1_6 + qg1_7 - qp1_7 - 5.17 = 0;
subject to cons_r1_8:
	r1_8 - r1_7 + qg1_8 - qp1_8 - 5.17 = 0;
subject to cons_r1_9:
	r1_9 - r1_8 + qg1_9 - qp1_9 - 5.17 = 0;
subject to cons_r1_10:
	r1_10 - r1_9 + qg1_10 - qp1_10 - 5.17 = 0;
subject to cons_r1_11:
	r1_11 - r1_10 + qg1_11 - qp1_11 - 5.17 = 0;
subject to cons_r1_12:
	r1_12 - r1_11 + qg1_12 - qp1_12 - 5.17 = 0;
subject to cons_r1_13:
	r1_13 - r1_12 + qg1_13 - qp1_13 - 5.17 = 0;
subject to cons_r1_14:
	r1_14 - r1_13 + qg1_14 - qp1_14 - 5.17 = 0;
subject to cons_r1_15:
	r1_15 - r1_14 + qg1_15 - qp1_15 - 5.17 = 0;
subject to cons_r1_16:
	r1_16 - r1_15 + qg1_16 - qp1_16 - 5.17 = 0;
subject to cons_r1_17:
	r1_17 - r1_16 + qg1_17 - qp1_17 - 5.17 = 0;
subject to cons_r1_18:
	r1_18 - r1_17 + qg1_18 - qp1_18 - 5.17 = 0;
subject to cons_r1_19:
	r1_19 - r1_18 + qg1_19 - qp1_19 - 5.17 = 0;
subject to cons_r1_20:
	r1_20 - r1_19 + qg1_20 - qp1_20 - 5.17 = 0;
subject to cons_r1_21:
	r1_21 - r1_20 + qg1_21 - qp1_21 - 5.17 = 0;
subject to cons_r1_22:
	r1_22 - r1_21 + qg1_22 - qp1_22 - 5.17 = 0;
subject to cons_r1_23:
	r1_23 - r1_22 + qg1_23 - qp1_23 - 5.17 = 0;
subject to cons_r1_24:
	r1_24 - r1_23 + qg1_24 - qp1_24 - 5.17 = 0;
subject to cons_d1_1:
	0 <= p11_1 + p21_1 + qh1_1 + qg1_1 - 1.33*qp1_1 - 578.0;
subject to cons_d1_2:
	0 <= p11_2 + p21_2 + qh1_2 + qg1_2 - 1.33*qp1_2 - 517.0;
subject to cons_d1_3:
	0 <= p11_3 + p21_3 + qh1_3 + qg1_3 - 1.33*qp1_3 - 461.0;
subject to cons_d1_4:
	0 <= p11_4 + p21_4 + qh1_4 + qg1_4 - 1.33*qp1_4 - 369.0;
subject to cons_d1_5:
	0 <= p11_5 + p21_5 + qh1_5 + qg1_5 - 1.33*qp1_5 - 299.0;
subject to cons_d1_6:
	0 <= p11_6 + p21_6 + qh1_6 + qg1_6 - 1.33*qp1_6 - 269.0;
subject to cons_d1_7:
	0 <= p11_7 + p21_7 + qh1_7 + qg1_7 - 1.33*qp1_7 - 370.0;
subject to cons_d1_8:
	0 <= p11_8 + p21_8 + qh1_8 + qg1_8 - 1.33*qp1_8 - 559.0;
subject to cons_d1_9:
	0 <= p11_9 + p21_9 + qh1_9 + qg1_9 - 1.33*qp1_9 - 689.0;
subject to cons_d1_10:
	0 <= p11_10 + p21_10 + qh1_10 + qg1_10 - 1.33*qp1_10 - 728.0;
subject to cons_d1_11:
	0 <= p11_11 + p21_11 + qh1_11 + qg1_11 - 1.33*qp1_11 - 683.0;
subject to cons_d1_12:
	0 <= p11_12 + p21_12 + qh1_12 + qg1_12 - 1.33*qp1_12 - 626.0;
subject to cons_d1_13:
	0 <= p11_13 + p21_13 + qh1_13 + qg1_13 - 1.33*qp1_13 - 619.0;
subject to cons_d1_14:
	0 <= p11_14 + p21_14 + qh1_14 + qg1_14 - 1.33*qp1_14 - 586.0;
subject to cons_d1_15:
	0 <= p11_15 + p21_15 + qh1_15 + qg1_15 - 1.33*qp1_15 - 582.0;
subject to cons_d1_16:
	0 <= p11_16 + p21_16 + qh1_16 + qg1_16 - 1.33*qp1_16 - 625.0;
subject to cons_d1_17:
	0 <= p11_17 + p21_17 + qh1_17 + qg1_17 - 1.33*qp1_17 - 821.0;
subject to cons_d1_18:
	0 <= p11_18 + p21_18 + qh1_18 + qg1_18 - 1.33*qp1_18 - 883.0;
subject to cons_d1_19:
	0 <= p11_19 + p21_19 + qh1_19 + qg1_19 - 1.33*qp1_19 - 768.0;
subject to cons_d1_20:
	0 <= p11_20 + p21_20 + qh1_20 + qg1_20 - 1.33*qp1_20 - 711.0;
subject to cons_d1_21:
	0 <= p11_21 + p21_21 + qh1_21 + qg1_21 - 1.33*qp1_21 - 677.0;
subject to cons_d1_22:
	0 <= p11_22 + p21_22 + qh1_22 + qg1_22 - 1.33*qp1_22 - 630.0;
subject to cons_d1_23:
	0 <= p11_23 + p21_23 + qh1_23 + qg1_23 - 1.33*qp1_23 - 545.0;
subject to cons_d1_24:
	0 <= p11_24 + p21_24 + qh1_24 + qg1_24 - 1.33*qp1_24 - 565.0;
subject to cons_qgtqp1_1:
	qg1_1 * qp1_1 = 0;
subject to cons_qgtqp1_2:
	qg1_2 * qp1_2 = 0;
subject to cons_qgtqp1_3:
	qg1_3 * qp1_3 = 0;
subject to cons_qgtqp1_4:
	qg1_4 * qp1_4 = 0;
subject to cons_qgtqp1_5:
	qg1_5 * qp1_5 = 0;
subject to cons_qgtqp1_6:
	qg1_6 * qp1_6 = 0;
subject to cons_qgtqp1_7:
	qg1_7 * qp1_7 = 0;
subject to cons_qgtqp1_8:
	qg1_8 * qp1_8 = 0;
subject to cons_qgtqp1_9:
	qg1_9 * qp1_9 = 0;
subject to cons_qgtqp1_10:
	qg1_10 * qp1_10 = 0;
subject to cons_qgtqp1_11:
	qg1_11 * qp1_11 = 0;
subject to cons_qgtqp1_12:
	qg1_12 * qp1_12 = 0;
subject to cons_qgtqp1_13:
	qg1_13 * qp1_13 = 0;
subject to cons_qgtqp1_14:
	qg1_14 * qp1_14 = 0;
subject to cons_qgtqp1_15:
	qg1_15 * qp1_15 = 0;
subject to cons_qgtqp1_16:
	qg1_16 * qp1_16 = 0;
subject to cons_qgtqp1_17:
	qg1_17 * qp1_17 = 0;
subject to cons_qgtqp1_18:
	qg1_18 * qp1_18 = 0;
subject to cons_qgtqp1_19:
	qg1_19 * qp1_19 = 0;
subject to cons_qgtqp1_20:
	qg1_20 * qp1_20 = 0;
subject to cons_qgtqp1_21:
	qg1_21 * qp1_21 = 0;
subject to cons_qgtqp1_22:
	qg1_22 * qp1_22 = 0;
subject to cons_qgtqp1_23:
	qg1_23 * qp1_23 = 0;
subject to cons_qgtqp1_24:
	qg1_24 * qp1_24 = 0;

solve;
display v0_24;
display r0_24;
display p11_1;
display p21_1;
display qh1_1;
display s1_1;
display qg1_1;
display qp1_1;
display v1_1;
display r1_1;
display p11_2;
display p21_2;
display qh1_2;
display s1_2;
display qg1_2;
display qp1_2;
display v1_2;
display r1_2;
display p11_3;
display p21_3;
display qh1_3;
display s1_3;
display qg1_3;
display qp1_3;
display v1_3;
display r1_3;
display p11_4;
display p21_4;
display qh1_4;
display s1_4;
display qg1_4;
display qp1_4;
display v1_4;
display r1_4;
display p11_5;
display p21_5;
display qh1_5;
display s1_5;
display qg1_5;
display qp1_5;
display v1_5;
display r1_5;
display p11_6;
display p21_6;
display qh1_6;
display s1_6;
display qg1_6;
display qp1_6;
display v1_6;
display r1_6;
display p11_7;
display p21_7;
display qh1_7;
display s1_7;
display qg1_7;
display qp1_7;
display v1_7;
display r1_7;
display p11_8;
display p21_8;
display qh1_8;
display s1_8;
display qg1_8;
display qp1_8;
display v1_8;
display r1_8;
display p11_9;
display p21_9;
display qh1_9;
display s1_9;
display qg1_9;
display qp1_9;
display v1_9;
display r1_9;
display p11_10;
display p21_10;
display qh1_10;
display s1_10;
display qg1_10;
display qp1_10;
display v1_10;
display r1_10;
display p11_11;
display p21_11;
display qh1_11;
display s1_11;
display qg1_11;
display qp1_11;
display v1_11;
display r1_11;
display p11_12;
display p21_12;
display qh1_12;
display s1_12;
display qg1_12;
display qp1_12;
display v1_12;
display r1_12;
display p11_13;
display p21_13;
display qh1_13;
display s1_13;
display qg1_13;
display qp1_13;
display v1_13;
display r1_13;
display p11_14;
display p21_14;
display qh1_14;
display s1_14;
display qg1_14;
display qp1_14;
display v1_14;
display r1_14;
display p11_15;
display p21_15;
display qh1_15;
display s1_15;
display qg1_15;
display qp1_15;
display v1_15;
display r1_15;
display p11_16;
display p21_16;
display qh1_16;
display s1_16;
display qg1_16;
display qp1_16;
display v1_16;
display r1_16;
display p11_17;
display p21_17;
display qh1_17;
display s1_17;
display qg1_17;
display qp1_17;
display v1_17;
display r1_17;
display p11_18;
display p21_18;
display qh1_18;
display s1_18;
display qg1_18;
display qp1_18;
display v1_18;
display r1_18;
display p11_19;
display p21_19;
display qh1_19;
display s1_19;
display qg1_19;
display qp1_19;
display v1_19;
display r1_19;
display p11_20;
display p21_20;
display qh1_20;
display s1_20;
display qg1_20;
display qp1_20;
display v1_20;
display r1_20;
display p11_21;
display p21_21;
display qh1_21;
display s1_21;
display qg1_21;
display qp1_21;
display v1_21;
display r1_21;
display p11_22;
display p21_22;
display qh1_22;
display s1_22;
display qg1_22;
display qp1_22;
display v1_22;
display r1_22;
display p11_23;
display p21_23;
display qh1_23;
display s1_23;
display qg1_23;
display qp1_23;
display v1_23;
display r1_23;
display p11_24;
display p21_24;
display qh1_24;
display s1_24;
display qg1_24;
display qp1_24;
display v1_24;
display r1_24;
display obj;
