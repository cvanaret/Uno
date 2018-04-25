#
#**************************
# SET UP THE INITIAL DATA *
#**************************
#   Problem :
#   *********
#   A totally separable nonconvex multi-commodity network problem
#   Source: p. 124 of
#   P.A. Steenbrink,
#   "Optimization of Transport Networks",
#   Wiley, 1974.
#   Note that Steenbrink does not give values for TZERO, CCR and alpha.
#   The problem has also been slightly perturbed by making NONZ >0, in
#   order to avoid undefined values for some elements and infinite
#   slope for others at the solution.
#   SIF input: Ph. Toint, June 1990.
#   classification ONR2-MY-468-108
#   The network has the following structure
#                           4
#                 4 ----------------- 2
#               / | \               / | \
#              /  |  \             /  |  \
#             /   |   \           /   |   \
#         10 /    |    \ 11    6 /    |    \ 1
#           /     |12   \       /     |5    \
#          /      |      \     /      |      \
#         /   16  |   18  \   /  17   |    3  \
#       6 ------- 9 ------- 8 ------- 7 ------- 1
#         \       |       /   \       |       /
#          \      |      /     \      |      /
#           \     |     /       \     |     /
#            \    |15  /         \    |8   / 2
#          13 \   |   / 14      9 \   |   /
#              \  |  /             \  |  /
#               \ | /               \ | /
#                 5 ----------------- 3
#                           7
#   Arcs appear in both directions.  The flows using the arc from the
#   lower numbered node to the higher numbred node are direct (D),
#   the flows in the other direction reverse (R)
#   The problem has 468 variables and 108 linear equality constraints.
#   The commodities are determined from the "trip-matrix" and are
#   numbered as as follows:
#     commodity 1 : flow of 2000 from 2 to 3
#     commodity 2 : flow of 2000 from 2 to 4
#     commodity 3 : flow of 1000 from 2 to 5
#     commodity 4 : flow of 1000 from 3 to 4
#     commodity 5 : flow of 2000 from 3 to 5
#     commodity 6 : flow of 1000 from 4 to 5
#     commodity 7 : flow of  200 from 3 to 2
#     commodity 8 : flow of  200 from 4 to 2
#     commodity 9 : flow of  100 from 5 to 2
#     commodity 10: flow of  100 from 4 to 3
#     commodity 11: flow of  200 from 5 to 3
#     commodity 12: flow of  100 from 5 to 4
#   Number of bi-directional arcs
#   Number of trips (commodities)
#   Lengths of the arcs
#   Other parameters
#   Constants
#   Arc parameters
#   Minimal capacities
#   Coefficients
#   Half investment cost for arc 4
#   Objective function
#   Flow conservation constraints
#   Node 1
#   Node 2
#   Node 3
#   Node 4
#   Node 5
#   Node 6
#   Node 7
#   Node 8
#   Node 9
#   End of the flow conservation constraints
#   Solution
	param narcs := 18;
	param ntrips := 12;
	param cost1 := 35.0;
	param cost2 := 40.0;
	param cost3 := 30.0;
	param cost4 := 100.0;
	param cost5 := 15.0;
	param cost6 := 55.0;
	param cost7 := 100.0;
	param cost8 := 25.0;
	param cost9 := 60.0;
	param cost10 := 35.0;
	param cost11 := 55.0;
	param cost12 := 15.0;
	param cost13 := 40.0;
	param cost14 := 60.0;
	param cost15 := 25.0;
	param cost16 := 30.0;
	param cost17 := 50.0;
	param cost18 := 50.0;
	param alph := 0.01;
	param tzero := 0.01;
	param ccr := 0.01;
	param nonz := 0.01;
	param mnonz := -1.0 * (0.01);
	param cmd1 := (0.0) + (0.01);
	param cmr1 := (0.0) + (0.01);
	param la1 := (35.0) * (0.01);
	param lt1 := (35.0) * (0.01);
	param lc1 := (35.0) * (0.01);
	param cmd2 := (0.0) + (0.01);
	param cmr2 := (0.0) + (0.01);
	param la2 := (40.0) * (0.01);
	param lt2 := (40.0) * (0.01);
	param lc2 := (40.0) * (0.01);
	param cmd3 := (0.0) + (0.01);
	param cmr3 := (0.0) + (0.01);
	param la3 := (30.0) * (0.01);
	param lt3 := (30.0) * (0.01);
	param lc3 := (30.0) * (0.01);
	param cmd4 := (0.0) + (0.01);
	param cmr4 := (0.0) + (0.01);
	param la4 := 0.5 * ((100.0) * (0.01));
	param lt4 := (100.0) * (0.01);
	param lc4 := (100.0) * (0.01);
	param cmd5 := (0.0) + (0.01);
	param cmr5 := (0.0) + (0.01);
	param la5 := (15.0) * (0.01);
	param lt5 := (15.0) * (0.01);
	param lc5 := (15.0) * (0.01);
	param cmd6 := (0.0) + (0.01);
	param cmr6 := (0.0) + (0.01);
	param la6 := (55.0) * (0.01);
	param lt6 := (55.0) * (0.01);
	param lc6 := (55.0) * (0.01);
	param cmd7 := (0.0) + (0.01);
	param cmr7 := (0.0) + (0.01);
	param la7 := (100.0) * (0.01);
	param lt7 := (100.0) * (0.01);
	param lc7 := (100.0) * (0.01);
	param cmd8 := (0.0) + (0.01);
	param cmr8 := (0.0) + (0.01);
	param la8 := (25.0) * (0.01);
	param lt8 := (25.0) * (0.01);
	param lc8 := (25.0) * (0.01);
	param cmd9 := (0.0) + (0.01);
	param cmr9 := (0.0) + (0.01);
	param la9 := (60.0) * (0.01);
	param lt9 := (60.0) * (0.01);
	param lc9 := (60.0) * (0.01);
	param cmd10 := (0.0) + (0.01);
	param cmr10 := (0.0) + (0.01);
	param la10 := (35.0) * (0.01);
	param lt10 := (35.0) * (0.01);
	param lc10 := (35.0) * (0.01);
	param cmd11 := (0.0) + (0.01);
	param cmr11 := (0.0) + (0.01);
	param la11 := (55.0) * (0.01);
	param lt11 := (55.0) * (0.01);
	param lc11 := (55.0) * (0.01);
	param cmd12 := (0.0) + (0.01);
	param cmr12 := (0.0) + (0.01);
	param la12 := (15.0) * (0.01);
	param lt12 := (15.0) * (0.01);
	param lc12 := (15.0) * (0.01);
	param cmd13 := (0.0) + (0.01);
	param cmr13 := (0.0) + (0.01);
	param la13 := (40.0) * (0.01);
	param lt13 := (40.0) * (0.01);
	param lc13 := (40.0) * (0.01);
	param cmd14 := (0.0) + (0.01);
	param cmr14 := (0.0) + (0.01);
	param la14 := (60.0) * (0.01);
	param lt14 := (60.0) * (0.01);
	param lc14 := (60.0) * (0.01);
	param cmd15 := (0.0) + (0.01);
	param cmr15 := (0.0) + (0.01);
	param la15 := (25.0) * (0.01);
	param lt15 := (25.0) * (0.01);
	param lc15 := (25.0) * (0.01);
	param cmd16 := (0.0) + (0.01);
	param cmr16 := (0.0) + (0.01);
	param la16 := (30.0) * (0.01);
	param lt16 := (30.0) * (0.01);
	param lc16 := (30.0) * (0.01);
	param cmd17 := (0.0) + (0.01);
	param cmr17 := (0.0) + (0.01);
	param la17 := (50.0) * (0.01);
	param lt17 := (50.0) * (0.01);
	param lc17 := (50.0) * (0.01);
	param cmd18 := (0.0) + (0.01);
	param cmr18 := (0.0) + (0.01);
	param la18 := (50.0) * (0.01);
	param lt18 := (50.0) * (0.01);
	param lc18 := (50.0) * (0.01);
	param shift := ((0.0) + (0.01)) + (-1.0 * (0.01));

	var cd1 >= 0.01 ,  := 0.1;
	var cr1 >= 0.01 ,  := 0.1;
	var d1_1 >= 0.0 ,  := 0.1;
	var r1_1 >= 0.0 ,  := 0.1;
	var d2_1 >= 0.0 ,  := 0.1;
	var r2_1 >= 0.0 ,  := 0.1;
	var d3_1 >= 0.0 ,  := 0.1;
	var r3_1 >= 0.0 ,  := 0.1;
	var d4_1 >= 0.0 ,  := 0.1;
	var r4_1 >= 0.0 ,  := 0.1;
	var d5_1 >= 0.0 ,  := 0.1;
	var r5_1 >= 0.0 ,  := 0.1;
	var d6_1 >= 0.0 ,  := 0.1;
	var r6_1 >= 0.0 ,  := 0.1;
	var d7_1 >= 0.0 ,  := 0.1;
	var r7_1 >= 0.0 ,  := 0.1;
	var d8_1 >= 0.0 ,  := 0.1;
	var r8_1 >= 0.0 ,  := 0.1;
	var d9_1 >= 0.0 ,  := 0.1;
	var r9_1 >= 0.0 ,  := 0.1;
	var d10_1 >= 0.0 ,  := 0.1;
	var r10_1 >= 0.0 ,  := 0.1;
	var d11_1 >= 0.0 ,  := 0.1;
	var r11_1 >= 0.0 ,  := 0.1;
	var d12_1 >= 0.0 ,  := 0.1;
	var r12_1 >= 0.0 ,  := 0.1;
	var cd2 >= 0.01 ,  := 0.1;
	var cr2 >= 0.01 ,  := 0.1;
	var d1_2 >= 0.0 ,  := 0.1;
	var r1_2 >= 0.0 ,  := 0.1;
	var d2_2 >= 0.0 ,  := 0.1;
	var r2_2 >= 0.0 ,  := 0.1;
	var d3_2 >= 0.0 ,  := 0.1;
	var r3_2 >= 0.0 ,  := 0.1;
	var d4_2 >= 0.0 ,  := 0.1;
	var r4_2 >= 0.0 ,  := 0.1;
	var d5_2 >= 0.0 ,  := 0.1;
	var r5_2 >= 0.0 ,  := 0.1;
	var d6_2 >= 0.0 ,  := 0.1;
	var r6_2 >= 0.0 ,  := 0.1;
	var d7_2 >= 0.0 ,  := 0.1;
	var r7_2 >= 0.0 ,  := 0.1;
	var d8_2 >= 0.0 ,  := 0.1;
	var r8_2 >= 0.0 ,  := 0.1;
	var d9_2 >= 0.0 ,  := 0.1;
	var r9_2 >= 0.0 ,  := 0.1;
	var d10_2 >= 0.0 ,  := 0.1;
	var r10_2 >= 0.0 ,  := 0.1;
	var d11_2 >= 0.0 ,  := 0.1;
	var r11_2 >= 0.0 ,  := 0.1;
	var d12_2 >= 0.0 ,  := 0.1;
	var r12_2 >= 0.0 ,  := 0.1;
	var cd3 >= 0.01 ,  := 0.1;
	var cr3 >= 0.01 ,  := 0.1;
	var d1_3 >= 0.0 ,  := 0.1;
	var r1_3 >= 0.0 ,  := 0.1;
	var d2_3 >= 0.0 ,  := 0.1;
	var r2_3 >= 0.0 ,  := 0.1;
	var d3_3 >= 0.0 ,  := 0.1;
	var r3_3 >= 0.0 ,  := 0.1;
	var d4_3 >= 0.0 ,  := 0.1;
	var r4_3 >= 0.0 ,  := 0.1;
	var d5_3 >= 0.0 ,  := 0.1;
	var r5_3 >= 0.0 ,  := 0.1;
	var d6_3 >= 0.0 ,  := 0.1;
	var r6_3 >= 0.0 ,  := 0.1;
	var d7_3 >= 0.0 ,  := 0.1;
	var r7_3 >= 0.0 ,  := 0.1;
	var d8_3 >= 0.0 ,  := 0.1;
	var r8_3 >= 0.0 ,  := 0.1;
	var d9_3 >= 0.0 ,  := 0.1;
	var r9_3 >= 0.0 ,  := 0.1;
	var d10_3 >= 0.0 ,  := 0.1;
	var r10_3 >= 0.0 ,  := 0.1;
	var d11_3 >= 0.0 ,  := 0.1;
	var r11_3 >= 0.0 ,  := 0.1;
	var d12_3 >= 0.0 ,  := 0.1;
	var r12_3 >= 0.0 ,  := 0.1;
	var cd4 >= 0.01 ,  := 0.1;
	var cr4 >= 0.01 ,  := 0.1;
	var d1_4 >= 0.0 ,  := 0.1;
	var r1_4 >= 0.0 ,  := 0.1;
	var d2_4 >= 0.0 ,  := 0.1;
	var r2_4 >= 0.0 ,  := 0.1;
	var d3_4 >= 0.0 ,  := 0.1;
	var r3_4 >= 0.0 ,  := 0.1;
	var d4_4 >= 0.0 ,  := 0.1;
	var r4_4 >= 0.0 ,  := 0.1;
	var d5_4 >= 0.0 ,  := 0.1;
	var r5_4 >= 0.0 ,  := 0.1;
	var d6_4 >= 0.0 ,  := 0.1;
	var r6_4 >= 0.0 ,  := 0.1;
	var d7_4 >= 0.0 ,  := 0.1;
	var r7_4 >= 0.0 ,  := 0.1;
	var d8_4 >= 0.0 ,  := 0.1;
	var r8_4 >= 0.0 ,  := 0.1;
	var d9_4 >= 0.0 ,  := 0.1;
	var r9_4 >= 0.0 ,  := 0.1;
	var d10_4 >= 0.0 ,  := 0.1;
	var r10_4 >= 0.0 ,  := 0.1;
	var d11_4 >= 0.0 ,  := 0.1;
	var r11_4 >= 0.0 ,  := 0.1;
	var d12_4 >= 0.0 ,  := 0.1;
	var r12_4 >= 0.0 ,  := 0.1;
	var cd5 >= 0.01 ,  := 0.1;
	var cr5 >= 0.01 ,  := 0.1;
	var d1_5 >= 0.0 ,  := 0.1;
	var r1_5 >= 0.0 ,  := 0.1;
	var d2_5 >= 0.0 ,  := 0.1;
	var r2_5 >= 0.0 ,  := 0.1;
	var d3_5 >= 0.0 ,  := 0.1;
	var r3_5 >= 0.0 ,  := 0.1;
	var d4_5 >= 0.0 ,  := 0.1;
	var r4_5 >= 0.0 ,  := 0.1;
	var d5_5 >= 0.0 ,  := 0.1;
	var r5_5 >= 0.0 ,  := 0.1;
	var d6_5 >= 0.0 ,  := 0.1;
	var r6_5 >= 0.0 ,  := 0.1;
	var d7_5 >= 0.0 ,  := 0.1;
	var r7_5 >= 0.0 ,  := 0.1;
	var d8_5 >= 0.0 ,  := 0.1;
	var r8_5 >= 0.0 ,  := 0.1;
	var d9_5 >= 0.0 ,  := 0.1;
	var r9_5 >= 0.0 ,  := 0.1;
	var d10_5 >= 0.0 ,  := 0.1;
	var r10_5 >= 0.0 ,  := 0.1;
	var d11_5 >= 0.0 ,  := 0.1;
	var r11_5 >= 0.0 ,  := 0.1;
	var d12_5 >= 0.0 ,  := 0.1;
	var r12_5 >= 0.0 ,  := 0.1;
	var cd6 >= 0.01 ,  := 0.1;
	var cr6 >= 0.01 ,  := 0.1;
	var d1_6 >= 0.0 ,  := 0.1;
	var r1_6 >= 0.0 ,  := 0.1;
	var d2_6 >= 0.0 ,  := 0.1;
	var r2_6 >= 0.0 ,  := 0.1;
	var d3_6 >= 0.0 ,  := 0.1;
	var r3_6 >= 0.0 ,  := 0.1;
	var d4_6 >= 0.0 ,  := 0.1;
	var r4_6 >= 0.0 ,  := 0.1;
	var d5_6 >= 0.0 ,  := 0.1;
	var r5_6 >= 0.0 ,  := 0.1;
	var d6_6 >= 0.0 ,  := 0.1;
	var r6_6 >= 0.0 ,  := 0.1;
	var d7_6 >= 0.0 ,  := 0.1;
	var r7_6 >= 0.0 ,  := 0.1;
	var d8_6 >= 0.0 ,  := 0.1;
	var r8_6 >= 0.0 ,  := 0.1;
	var d9_6 >= 0.0 ,  := 0.1;
	var r9_6 >= 0.0 ,  := 0.1;
	var d10_6 >= 0.0 ,  := 0.1;
	var r10_6 >= 0.0 ,  := 0.1;
	var d11_6 >= 0.0 ,  := 0.1;
	var r11_6 >= 0.0 ,  := 0.1;
	var d12_6 >= 0.0 ,  := 0.1;
	var r12_6 >= 0.0 ,  := 0.1;
	var cd7 >= 0.01 ,  := 0.1;
	var cr7 >= 0.01 ,  := 0.1;
	var d1_7 >= 0.0 ,  := 0.1;
	var r1_7 >= 0.0 ,  := 0.1;
	var d2_7 >= 0.0 ,  := 0.1;
	var r2_7 >= 0.0 ,  := 0.1;
	var d3_7 >= 0.0 ,  := 0.1;
	var r3_7 >= 0.0 ,  := 0.1;
	var d4_7 >= 0.0 ,  := 0.1;
	var r4_7 >= 0.0 ,  := 0.1;
	var d5_7 >= 0.0 ,  := 0.1;
	var r5_7 >= 0.0 ,  := 0.1;
	var d6_7 >= 0.0 ,  := 0.1;
	var r6_7 >= 0.0 ,  := 0.1;
	var d7_7 >= 0.0 ,  := 0.1;
	var r7_7 >= 0.0 ,  := 0.1;
	var d8_7 >= 0.0 ,  := 0.1;
	var r8_7 >= 0.0 ,  := 0.1;
	var d9_7 >= 0.0 ,  := 0.1;
	var r9_7 >= 0.0 ,  := 0.1;
	var d10_7 >= 0.0 ,  := 0.1;
	var r10_7 >= 0.0 ,  := 0.1;
	var d11_7 >= 0.0 ,  := 0.1;
	var r11_7 >= 0.0 ,  := 0.1;
	var d12_7 >= 0.0 ,  := 0.1;
	var r12_7 >= 0.0 ,  := 0.1;
	var cd8 >= 0.01 ,  := 0.1;
	var cr8 >= 0.01 ,  := 0.1;
	var d1_8 >= 0.0 ,  := 0.1;
	var r1_8 >= 0.0 ,  := 0.1;
	var d2_8 >= 0.0 ,  := 0.1;
	var r2_8 >= 0.0 ,  := 0.1;
	var d3_8 >= 0.0 ,  := 0.1;
	var r3_8 >= 0.0 ,  := 0.1;
	var d4_8 >= 0.0 ,  := 0.1;
	var r4_8 >= 0.0 ,  := 0.1;
	var d5_8 >= 0.0 ,  := 0.1;
	var r5_8 >= 0.0 ,  := 0.1;
	var d6_8 >= 0.0 ,  := 0.1;
	var r6_8 >= 0.0 ,  := 0.1;
	var d7_8 >= 0.0 ,  := 0.1;
	var r7_8 >= 0.0 ,  := 0.1;
	var d8_8 >= 0.0 ,  := 0.1;
	var r8_8 >= 0.0 ,  := 0.1;
	var d9_8 >= 0.0 ,  := 0.1;
	var r9_8 >= 0.0 ,  := 0.1;
	var d10_8 >= 0.0 ,  := 0.1;
	var r10_8 >= 0.0 ,  := 0.1;
	var d11_8 >= 0.0 ,  := 0.1;
	var r11_8 >= 0.0 ,  := 0.1;
	var d12_8 >= 0.0 ,  := 0.1;
	var r12_8 >= 0.0 ,  := 0.1;
	var cd9 >= 0.01 ,  := 0.1;
	var cr9 >= 0.01 ,  := 0.1;
	var d1_9 >= 0.0 ,  := 0.1;
	var r1_9 >= 0.0 ,  := 0.1;
	var d2_9 >= 0.0 ,  := 0.1;
	var r2_9 >= 0.0 ,  := 0.1;
	var d3_9 >= 0.0 ,  := 0.1;
	var r3_9 >= 0.0 ,  := 0.1;
	var d4_9 >= 0.0 ,  := 0.1;
	var r4_9 >= 0.0 ,  := 0.1;
	var d5_9 >= 0.0 ,  := 0.1;
	var r5_9 >= 0.0 ,  := 0.1;
	var d6_9 >= 0.0 ,  := 0.1;
	var r6_9 >= 0.0 ,  := 0.1;
	var d7_9 >= 0.0 ,  := 0.1;
	var r7_9 >= 0.0 ,  := 0.1;
	var d8_9 >= 0.0 ,  := 0.1;
	var r8_9 >= 0.0 ,  := 0.1;
	var d9_9 >= 0.0 ,  := 0.1;
	var r9_9 >= 0.0 ,  := 0.1;
	var d10_9 >= 0.0 ,  := 0.1;
	var r10_9 >= 0.0 ,  := 0.1;
	var d11_9 >= 0.0 ,  := 0.1;
	var r11_9 >= 0.0 ,  := 0.1;
	var d12_9 >= 0.0 ,  := 0.1;
	var r12_9 >= 0.0 ,  := 0.1;
	var cd10 >= 0.01 ,  := 0.1;
	var cr10 >= 0.01 ,  := 0.1;
	var d1_10 >= 0.0 ,  := 0.1;
	var r1_10 >= 0.0 ,  := 0.1;
	var d2_10 >= 0.0 ,  := 0.1;
	var r2_10 >= 0.0 ,  := 0.1;
	var d3_10 >= 0.0 ,  := 0.1;
	var r3_10 >= 0.0 ,  := 0.1;
	var d4_10 >= 0.0 ,  := 0.1;
	var r4_10 >= 0.0 ,  := 0.1;
	var d5_10 >= 0.0 ,  := 0.1;
	var r5_10 >= 0.0 ,  := 0.1;
	var d6_10 >= 0.0 ,  := 0.1;
	var r6_10 >= 0.0 ,  := 0.1;
	var d7_10 >= 0.0 ,  := 0.1;
	var r7_10 >= 0.0 ,  := 0.1;
	var d8_10 >= 0.0 ,  := 0.1;
	var r8_10 >= 0.0 ,  := 0.1;
	var d9_10 >= 0.0 ,  := 0.1;
	var r9_10 >= 0.0 ,  := 0.1;
	var d10_10 >= 0.0 ,  := 0.1;
	var r10_10 >= 0.0 ,  := 0.1;
	var d11_10 >= 0.0 ,  := 0.1;
	var r11_10 >= 0.0 ,  := 0.1;
	var d12_10 >= 0.0 ,  := 0.1;
	var r12_10 >= 0.0 ,  := 0.1;
	var cd11 >= 0.01 ,  := 0.1;
	var cr11 >= 0.01 ,  := 0.1;
	var d1_11 >= 0.0 ,  := 0.1;
	var r1_11 >= 0.0 ,  := 0.1;
	var d2_11 >= 0.0 ,  := 0.1;
	var r2_11 >= 0.0 ,  := 0.1;
	var d3_11 >= 0.0 ,  := 0.1;
	var r3_11 >= 0.0 ,  := 0.1;
	var d4_11 >= 0.0 ,  := 0.1;
	var r4_11 >= 0.0 ,  := 0.1;
	var d5_11 >= 0.0 ,  := 0.1;
	var r5_11 >= 0.0 ,  := 0.1;
	var d6_11 >= 0.0 ,  := 0.1;
	var r6_11 >= 0.0 ,  := 0.1;
	var d7_11 >= 0.0 ,  := 0.1;
	var r7_11 >= 0.0 ,  := 0.1;
	var d8_11 >= 0.0 ,  := 0.1;
	var r8_11 >= 0.0 ,  := 0.1;
	var d9_11 >= 0.0 ,  := 0.1;
	var r9_11 >= 0.0 ,  := 0.1;
	var d10_11 >= 0.0 ,  := 0.1;
	var r10_11 >= 0.0 ,  := 0.1;
	var d11_11 >= 0.0 ,  := 0.1;
	var r11_11 >= 0.0 ,  := 0.1;
	var d12_11 >= 0.0 ,  := 0.1;
	var r12_11 >= 0.0 ,  := 0.1;
	var cd12 >= 0.01 ,  := 0.1;
	var cr12 >= 0.01 ,  := 0.1;
	var d1_12 >= 0.0 ,  := 0.1;
	var r1_12 >= 0.0 ,  := 0.1;
	var d2_12 >= 0.0 ,  := 0.1;
	var r2_12 >= 0.0 ,  := 0.1;
	var d3_12 >= 0.0 ,  := 0.1;
	var r3_12 >= 0.0 ,  := 0.1;
	var d4_12 >= 0.0 ,  := 0.1;
	var r4_12 >= 0.0 ,  := 0.1;
	var d5_12 >= 0.0 ,  := 0.1;
	var r5_12 >= 0.0 ,  := 0.1;
	var d6_12 >= 0.0 ,  := 0.1;
	var r6_12 >= 0.0 ,  := 0.1;
	var d7_12 >= 0.0 ,  := 0.1;
	var r7_12 >= 0.0 ,  := 0.1;
	var d8_12 >= 0.0 ,  := 0.1;
	var r8_12 >= 0.0 ,  := 0.1;
	var d9_12 >= 0.0 ,  := 0.1;
	var r9_12 >= 0.0 ,  := 0.1;
	var d10_12 >= 0.0 ,  := 0.1;
	var r10_12 >= 0.0 ,  := 0.1;
	var d11_12 >= 0.0 ,  := 0.1;
	var r11_12 >= 0.0 ,  := 0.1;
	var d12_12 >= 0.0 ,  := 0.1;
	var r12_12 >= 0.0 ,  := 0.1;
	var cd13 >= 0.01 ,  := 0.1;
	var cr13 >= 0.01 ,  := 0.1;
	var d1_13 >= 0.0 ,  := 0.1;
	var r1_13 >= 0.0 ,  := 0.1;
	var d2_13 >= 0.0 ,  := 0.1;
	var r2_13 >= 0.0 ,  := 0.1;
	var d3_13 >= 0.0 ,  := 0.1;
	var r3_13 >= 0.0 ,  := 0.1;
	var d4_13 >= 0.0 ,  := 0.1;
	var r4_13 >= 0.0 ,  := 0.1;
	var d5_13 >= 0.0 ,  := 0.1;
	var r5_13 >= 0.0 ,  := 0.1;
	var d6_13 >= 0.0 ,  := 0.1;
	var r6_13 >= 0.0 ,  := 0.1;
	var d7_13 >= 0.0 ,  := 0.1;
	var r7_13 >= 0.0 ,  := 0.1;
	var d8_13 >= 0.0 ,  := 0.1;
	var r8_13 >= 0.0 ,  := 0.1;
	var d9_13 >= 0.0 ,  := 0.1;
	var r9_13 >= 0.0 ,  := 0.1;
	var d10_13 >= 0.0 ,  := 0.1;
	var r10_13 >= 0.0 ,  := 0.1;
	var d11_13 >= 0.0 ,  := 0.1;
	var r11_13 >= 0.0 ,  := 0.1;
	var d12_13 >= 0.0 ,  := 0.1;
	var r12_13 >= 0.0 ,  := 0.1;
	var cd14 >= 0.01 ,  := 0.1;
	var cr14 >= 0.01 ,  := 0.1;
	var d1_14 >= 0.0 ,  := 0.1;
	var r1_14 >= 0.0 ,  := 0.1;
	var d2_14 >= 0.0 ,  := 0.1;
	var r2_14 >= 0.0 ,  := 0.1;
	var d3_14 >= 0.0 ,  := 0.1;
	var r3_14 >= 0.0 ,  := 0.1;
	var d4_14 >= 0.0 ,  := 0.1;
	var r4_14 >= 0.0 ,  := 0.1;
	var d5_14 >= 0.0 ,  := 0.1;
	var r5_14 >= 0.0 ,  := 0.1;
	var d6_14 >= 0.0 ,  := 0.1;
	var r6_14 >= 0.0 ,  := 0.1;
	var d7_14 >= 0.0 ,  := 0.1;
	var r7_14 >= 0.0 ,  := 0.1;
	var d8_14 >= 0.0 ,  := 0.1;
	var r8_14 >= 0.0 ,  := 0.1;
	var d9_14 >= 0.0 ,  := 0.1;
	var r9_14 >= 0.0 ,  := 0.1;
	var d10_14 >= 0.0 ,  := 0.1;
	var r10_14 >= 0.0 ,  := 0.1;
	var d11_14 >= 0.0 ,  := 0.1;
	var r11_14 >= 0.0 ,  := 0.1;
	var d12_14 >= 0.0 ,  := 0.1;
	var r12_14 >= 0.0 ,  := 0.1;
	var cd15 >= 0.01 ,  := 0.1;
	var cr15 >= 0.01 ,  := 0.1;
	var d1_15 >= 0.0 ,  := 0.1;
	var r1_15 >= 0.0 ,  := 0.1;
	var d2_15 >= 0.0 ,  := 0.1;
	var r2_15 >= 0.0 ,  := 0.1;
	var d3_15 >= 0.0 ,  := 0.1;
	var r3_15 >= 0.0 ,  := 0.1;
	var d4_15 >= 0.0 ,  := 0.1;
	var r4_15 >= 0.0 ,  := 0.1;
	var d5_15 >= 0.0 ,  := 0.1;
	var r5_15 >= 0.0 ,  := 0.1;
	var d6_15 >= 0.0 ,  := 0.1;
	var r6_15 >= 0.0 ,  := 0.1;
	var d7_15 >= 0.0 ,  := 0.1;
	var r7_15 >= 0.0 ,  := 0.1;
	var d8_15 >= 0.0 ,  := 0.1;
	var r8_15 >= 0.0 ,  := 0.1;
	var d9_15 >= 0.0 ,  := 0.1;
	var r9_15 >= 0.0 ,  := 0.1;
	var d10_15 >= 0.0 ,  := 0.1;
	var r10_15 >= 0.0 ,  := 0.1;
	var d11_15 >= 0.0 ,  := 0.1;
	var r11_15 >= 0.0 ,  := 0.1;
	var d12_15 >= 0.0 ,  := 0.1;
	var r12_15 >= 0.0 ,  := 0.1;
	var cd16 >= 0.01 ,  := 0.1;
	var cr16 >= 0.01 ,  := 0.1;
	var d1_16 >= 0.0 ,  := 0.1;
	var r1_16 >= 0.0 ,  := 0.1;
	var d2_16 >= 0.0 ,  := 0.1;
	var r2_16 >= 0.0 ,  := 0.1;
	var d3_16 >= 0.0 ,  := 0.1;
	var r3_16 >= 0.0 ,  := 0.1;
	var d4_16 >= 0.0 ,  := 0.1;
	var r4_16 >= 0.0 ,  := 0.1;
	var d5_16 >= 0.0 ,  := 0.1;
	var r5_16 >= 0.0 ,  := 0.1;
	var d6_16 >= 0.0 ,  := 0.1;
	var r6_16 >= 0.0 ,  := 0.1;
	var d7_16 >= 0.0 ,  := 0.1;
	var r7_16 >= 0.0 ,  := 0.1;
	var d8_16 >= 0.0 ,  := 0.1;
	var r8_16 >= 0.0 ,  := 0.1;
	var d9_16 >= 0.0 ,  := 0.1;
	var r9_16 >= 0.0 ,  := 0.1;
	var d10_16 >= 0.0 ,  := 0.1;
	var r10_16 >= 0.0 ,  := 0.1;
	var d11_16 >= 0.0 ,  := 0.1;
	var r11_16 >= 0.0 ,  := 0.1;
	var d12_16 >= 0.0 ,  := 0.1;
	var r12_16 >= 0.0 ,  := 0.1;
	var cd17 >= 0.01 ,  := 0.1;
	var cr17 >= 0.01 ,  := 0.1;
	var d1_17 >= 0.0 ,  := 0.1;
	var r1_17 >= 0.0 ,  := 0.1;
	var d2_17 >= 0.0 ,  := 0.1;
	var r2_17 >= 0.0 ,  := 0.1;
	var d3_17 >= 0.0 ,  := 0.1;
	var r3_17 >= 0.0 ,  := 0.1;
	var d4_17 >= 0.0 ,  := 0.1;
	var r4_17 >= 0.0 ,  := 0.1;
	var d5_17 >= 0.0 ,  := 0.1;
	var r5_17 >= 0.0 ,  := 0.1;
	var d6_17 >= 0.0 ,  := 0.1;
	var r6_17 >= 0.0 ,  := 0.1;
	var d7_17 >= 0.0 ,  := 0.1;
	var r7_17 >= 0.0 ,  := 0.1;
	var d8_17 >= 0.0 ,  := 0.1;
	var r8_17 >= 0.0 ,  := 0.1;
	var d9_17 >= 0.0 ,  := 0.1;
	var r9_17 >= 0.0 ,  := 0.1;
	var d10_17 >= 0.0 ,  := 0.1;
	var r10_17 >= 0.0 ,  := 0.1;
	var d11_17 >= 0.0 ,  := 0.1;
	var r11_17 >= 0.0 ,  := 0.1;
	var d12_17 >= 0.0 ,  := 0.1;
	var r12_17 >= 0.0 ,  := 0.1;
	var cd18 >= 0.01 ,  := 0.1;
	var cr18 >= 0.01 ,  := 0.1;
	var d1_18 >= 0.0 ,  := 0.1;
	var r1_18 >= 0.0 ,  := 0.1;
	var d2_18 >= 0.0 ,  := 0.1;
	var r2_18 >= 0.0 ,  := 0.1;
	var d3_18 >= 0.0 ,  := 0.1;
	var r3_18 >= 0.0 ,  := 0.1;
	var d4_18 >= 0.0 ,  := 0.1;
	var r4_18 >= 0.0 ,  := 0.1;
	var d5_18 >= 0.0 ,  := 0.1;
	var r5_18 >= 0.0 ,  := 0.1;
	var d6_18 >= 0.0 ,  := 0.1;
	var r6_18 >= 0.0 ,  := 0.1;
	var d7_18 >= 0.0 ,  := 0.1;
	var r7_18 >= 0.0 ,  := 0.1;
	var d8_18 >= 0.0 ,  := 0.1;
	var r8_18 >= 0.0 ,  := 0.1;
	var d9_18 >= 0.0 ,  := 0.1;
	var r9_18 >= 0.0 ,  := 0.1;
	var d10_18 >= 0.0 ,  := 0.1;
	var r10_18 >= 0.0 ,  := 0.1;
	var d11_18 >= 0.0 ,  := 0.1;
	var r11_18 >= 0.0 ,  := 0.1;
	var d12_18 >= 0.0 ,  := 0.1;
	var r12_18 >= 0.0 ,  := 0.1;

minimize obj:
	0.35000000000000003 * (d11_1+d12_1) + 0.35000000000000003 * 
	((d11_1+d12_1)*(d11_1+d12_1)*d11_1+d12_1) / ((cd1)*cd1) + 0.35000000000000003 * 
	(r11_1+r12_1) + 0.35000000000000003 * ((r11_1+r12_1)*(r11_1+r12_1)*r11_1+r12_1) 
	/ ((cr1)*cr1) + 0.35000000000000003*(sqrt((cd1-0.0))) + 
	0.35000000000000003*(sqrt((cr1-0.0))) + 0.4 * (d11_2+d12_2) + 0.4 * 
	((d11_2+d12_2)*(d11_2+d12_2)*d11_2+d12_2) / ((cd2)*cd2) + 0.4 * (r11_2+r12_2) + 
	0.4 * ((r11_2+r12_2)*(r11_2+r12_2)*r11_2+r12_2) / ((cr2)*cr2) + 
	0.4*(sqrt((cd2-0.0))) + 0.4*(sqrt((cr2-0.0))) + 0.3 * (d11_3+d12_3) + 0.3 * 
	((d11_3+d12_3)*(d11_3+d12_3)*d11_3+d12_3) / ((cd3)*cd3) + 0.3 * (r11_3+r12_3) + 
	0.3 * ((r11_3+r12_3)*(r11_3+r12_3)*r11_3+r12_3) / ((cr3)*cr3) + 
	0.3*(sqrt((cd3-0.0))) + 0.3*(sqrt((cr3-0.0))) + 1.0 * (d11_4+d12_4) + 1.0 * 
	((d11_4+d12_4)*(d11_4+d12_4)*d11_4+d12_4) / ((cd4)*cd4) + 1.0 * (r11_4+r12_4) + 
	1.0 * ((r11_4+r12_4)*(r11_4+r12_4)*r11_4+r12_4) / ((cr4)*cr4) + 
	0.5*(sqrt((cd4-0.0))) + 0.5*(sqrt((cr4-0.0))) + 0.15 * (d11_5+d12_5) + 0.15 * 
	((d11_5+d12_5)*(d11_5+d12_5)*d11_5+d12_5) / ((cd5)*cd5) + 0.15 * (r11_5+r12_5) 
	+ 0.15 * ((r11_5+r12_5)*(r11_5+r12_5)*r11_5+r12_5) / ((cr5)*cr5) + 
	0.15*(sqrt((cd5-0.0))) + 0.15*(sqrt((cr5-0.0))) + 0.55 * (d11_6+d12_6) + 0.55 * 
	((d11_6+d12_6)*(d11_6+d12_6)*d11_6+d12_6) / ((cd6)*cd6) + 0.55 * (r11_6+r12_6) 
	+ 0.55 * ((r11_6+r12_6)*(r11_6+r12_6)*r11_6+r12_6) / ((cr6)*cr6) + 
	0.55*(sqrt((cd6-0.0))) + 0.55*(sqrt((cr6-0.0))) + 1.0 * (d11_7+d12_7) + 1.0 * 
	((d11_7+d12_7)*(d11_7+d12_7)*d11_7+d12_7) / ((cd7)*cd7) + 1.0 * (r11_7+r12_7) + 
	1.0 * ((r11_7+r12_7)*(r11_7+r12_7)*r11_7+r12_7) / ((cr7)*cr7) + 
	(sqrt((cd7-0.0))) + (sqrt((cr7-0.0))) + 0.25 * (d11_8+d12_8) + 0.25 * 
	((d11_8+d12_8)*(d11_8+d12_8)*d11_8+d12_8) / ((cd8)*cd8) + 0.25 * (r11_8+r12_8) 
	+ 0.25 * ((r11_8+r12_8)*(r11_8+r12_8)*r11_8+r12_8) / ((cr8)*cr8) + 
	0.25*(sqrt((cd8-0.0))) + 0.25*(sqrt((cr8-0.0))) + 0.6 * (d11_9+d12_9) + 0.6 * 
	((d11_9+d12_9)*(d11_9+d12_9)*d11_9+d12_9) / ((cd9)*cd9) + 0.6 * (r11_9+r12_9) + 
	0.6 * ((r11_9+r12_9)*(r11_9+r12_9)*r11_9+r12_9) / ((cr9)*cr9) + 
	0.6*(sqrt((cd9-0.0))) + 0.6*(sqrt((cr9-0.0))) + 0.35000000000000003 * 
	(d11_10+d12_10) + 0.35000000000000003 * 
	((d11_10+d12_10)*(d11_10+d12_10)*d11_10+d12_10) / ((cd10)*cd10) + 
	0.35000000000000003 * (r11_10+r12_10) + 0.35000000000000003 * 
	((r11_10+r12_10)*(r11_10+r12_10)*r11_10+r12_10) / ((cr10)*cr10) + 
	0.35000000000000003*(sqrt((cd10-0.0))) + 0.35000000000000003*(sqrt((cr10-0.0))) 
	+ 0.55 * (d11_11+d12_11) + 0.55 * 
	((d11_11+d12_11)*(d11_11+d12_11)*d11_11+d12_11) / ((cd11)*cd11) + 0.55 * 
	(r11_11+r12_11) + 0.55 * ((r11_11+r12_11)*(r11_11+r12_11)*r11_11+r12_11) / 
	((cr11)*cr11) + 0.55*(sqrt((cd11-0.0))) + 0.55*(sqrt((cr11-0.0))) + 0.15 * 
	(d11_12+d12_12) + 0.15 * ((d11_12+d12_12)*(d11_12+d12_12)*d11_12+d12_12) / 
	((cd12)*cd12) + 0.15 * (r11_12+r12_12) + 0.15 * 
	((r11_12+r12_12)*(r11_12+r12_12)*r11_12+r12_12) / ((cr12)*cr12) + 
	0.15*(sqrt((cd12-0.0))) + 0.15*(sqrt((cr12-0.0))) + 0.4 * (d11_13+d12_13) + 0.4 
	* ((d11_13+d12_13)*(d11_13+d12_13)*d11_13+d12_13) / ((cd13)*cd13) + 0.4 * 
	(r11_13+r12_13) + 0.4 * ((r11_13+r12_13)*(r11_13+r12_13)*r11_13+r12_13) / 
	((cr13)*cr13) + 0.4*(sqrt((cd13-0.0))) + 0.4*(sqrt((cr13-0.0))) + 0.6 * 
	(d11_14+d12_14) + 0.6 * ((d11_14+d12_14)*(d11_14+d12_14)*d11_14+d12_14) / 
	((cd14)*cd14) + 0.6 * (r11_14+r12_14) + 0.6 * 
	((r11_14+r12_14)*(r11_14+r12_14)*r11_14+r12_14) / ((cr14)*cr14) + 
	0.6*(sqrt((cd14-0.0))) + 0.6*(sqrt((cr14-0.0))) + 0.25 * (d11_15+d12_15) + 0.25 
	* ((d11_15+d12_15)*(d11_15+d12_15)*d11_15+d12_15) / ((cd15)*cd15) + 0.25 * 
	(r11_15+r12_15) + 0.25 * ((r11_15+r12_15)*(r11_15+r12_15)*r11_15+r12_15) / 
	((cr15)*cr15) + 0.25*(sqrt((cd15-0.0))) + 0.25*(sqrt((cr15-0.0))) + 0.3 * 
	(d11_16+d12_16) + 0.3 * ((d11_16+d12_16)*(d11_16+d12_16)*d11_16+d12_16) / 
	((cd16)*cd16) + 0.3 * (r11_16+r12_16) + 0.3 * 
	((r11_16+r12_16)*(r11_16+r12_16)*r11_16+r12_16) / ((cr16)*cr16) + 
	0.3*(sqrt((cd16-0.0))) + 0.3*(sqrt((cr16-0.0))) + 0.5 * (d11_17+d12_17) + 0.5 * 
	((d11_17+d12_17)*(d11_17+d12_17)*d11_17+d12_17) / ((cd17)*cd17) + 0.5 * 
	(r11_17+r12_17) + 0.5 * ((r11_17+r12_17)*(r11_17+r12_17)*r11_17+r12_17) / 
	((cr17)*cr17) + 0.5*(sqrt((cd17-0.0))) + 0.5*(sqrt((cr17-0.0))) + 0.5 * 
	(d11_18+d12_18) + 0.5 * ((d11_18+d12_18)*(d11_18+d12_18)*d11_18+d12_18) / 
	((cd18)*cd18) + 0.5 * (r11_18+r12_18) + 0.5 * 
	((r11_18+r12_18)*(r11_18+r12_18)*r11_18+r12_18) / ((cr18)*cr18) + 
	0.5*(sqrt((cd18-0.0))) + 0.5*(sqrt((cr18-0.0)));

subject to n11:
	-d1_1 + r1_1 - d1_2 + r1_2 - d1_3 + r1_3 = 0;
subject to n21:
	-d1_4 + r1_4 - d1_5 + r1_5 - d1_6 + r1_6 + d1_1 - r1_1 + 2000.0 = 0;
subject to n31:
	-d1_7 + r1_7 - d1_8 + r1_8 - d1_9 + r1_9 + d1_2 - r1_2 - 2000.0 = 0;
subject to n41:
	-d1_10 + r1_10 - d1_11 + r1_11 - d1_12 + r1_12 + d1_4 - r1_4 = 0;
subject to n51:
	-d1_13 + r1_13 - d1_14 + r1_14 - d1_15 + r1_15 + d1_7 - r1_7 = 0;
subject to n61:
	-d1_16 + r1_16 + d1_10 - r1_10 + d1_13 - r1_13 = 0;
subject to n71:
	-d1_17 + r1_17 + d1_3 - r1_3 + d1_5 - r1_5 + d1_8 - r1_8 = 0;
subject to n81:
	-d1_18 + r1_18 + d1_6 - r1_6 + d1_9 - r1_9 + d1_11 - r1_11 + d1_14 - r1_14 + 
	d1_17 - r1_17 = 0;
subject to n91:
	d1_12 - r1_12 + d1_15 - r1_15 + d1_16 - r1_16 + d1_18 - r1_18 = 0;
subject to n12:
	-d2_1 + r2_1 - d2_2 + r2_2 - d2_3 + r2_3 = 0;
subject to n22:
	-d2_4 + r2_4 - d2_5 + r2_5 - d2_6 + r2_6 + d2_1 - r2_1 + 2000.0 = 0;
subject to n32:
	-d2_7 + r2_7 - d2_8 + r2_8 - d2_9 + r2_9 + d2_2 - r2_2 = 0;
subject to n42:
	-d2_10 + r2_10 - d2_11 + r2_11 - d2_12 + r2_12 + d2_4 - r2_4 - 2000.0 = 0;
subject to n52:
	-d2_13 + r2_13 - d2_14 + r2_14 - d2_15 + r2_15 + d2_7 - r2_7 = 0;
subject to n62:
	-d2_16 + r2_16 + d2_10 - r2_10 + d2_13 - r2_13 = 0;
subject to n72:
	-d2_17 + r2_17 + d2_3 - r2_3 + d2_5 - r2_5 + d2_8 - r2_8 = 0;
subject to n82:
	-d2_18 + r2_18 + d2_6 - r2_6 + d2_9 - r2_9 + d2_11 - r2_11 + d2_14 - r2_14 + 
	d2_17 - r2_17 = 0;
subject to n92:
	d2_12 - r2_12 + d2_15 - r2_15 + d2_16 - r2_16 + d2_18 - r2_18 = 0;
subject to n13:
	-d3_1 + r3_1 - d3_2 + r3_2 - d3_3 + r3_3 = 0;
subject to n23:
	-d3_4 + r3_4 - d3_5 + r3_5 - d3_6 + r3_6 + d3_1 - r3_1 + 1000.0 = 0;
subject to n33:
	-d3_7 + r3_7 - d3_8 + r3_8 - d3_9 + r3_9 + d3_2 - r3_2 = 0;
subject to n43:
	-d3_10 + r3_10 - d3_11 + r3_11 - d3_12 + r3_12 + d3_4 - r3_4 = 0;
subject to n53:
	-d3_13 + r3_13 - d3_14 + r3_14 - d3_15 + r3_15 + d3_7 - r3_7 - 1000.0 = 0;
subject to n63:
	-d3_16 + r3_16 + d3_10 - r3_10 + d3_13 - r3_13 = 0;
subject to n73:
	-d3_17 + r3_17 + d3_3 - r3_3 + d3_5 - r3_5 + d3_8 - r3_8 = 0;
subject to n83:
	-d3_18 + r3_18 + d3_6 - r3_6 + d3_9 - r3_9 + d3_11 - r3_11 + d3_14 - r3_14 + 
	d3_17 - r3_17 = 0;
subject to n93:
	d3_12 - r3_12 + d3_15 - r3_15 + d3_16 - r3_16 + d3_18 - r3_18 = 0;
subject to n14:
	-d4_1 + r4_1 - d4_2 + r4_2 - d4_3 + r4_3 = 0;
subject to n24:
	-d4_4 + r4_4 - d4_5 + r4_5 - d4_6 + r4_6 + d4_1 - r4_1 = 0;
subject to n34:
	-d4_7 + r4_7 - d4_8 + r4_8 - d4_9 + r4_9 + d4_2 - r4_2 + 1000.0 = 0;
subject to n44:
	-d4_10 + r4_10 - d4_11 + r4_11 - d4_12 + r4_12 + d4_4 - r4_4 - 1000.0 = 0;
subject to n54:
	-d4_13 + r4_13 - d4_14 + r4_14 - d4_15 + r4_15 + d4_7 - r4_7 = 0;
subject to n64:
	-d4_16 + r4_16 + d4_10 - r4_10 + d4_13 - r4_13 = 0;
subject to n74:
	-d4_17 + r4_17 + d4_3 - r4_3 + d4_5 - r4_5 + d4_8 - r4_8 = 0;
subject to n84:
	-d4_18 + r4_18 + d4_6 - r4_6 + d4_9 - r4_9 + d4_11 - r4_11 + d4_14 - r4_14 + 
	d4_17 - r4_17 = 0;
subject to n94:
	d4_12 - r4_12 + d4_15 - r4_15 + d4_16 - r4_16 + d4_18 - r4_18 = 0;
subject to n15:
	-d5_1 + r5_1 - d5_2 + r5_2 - d5_3 + r5_3 = 0;
subject to n25:
	-d5_4 + r5_4 - d5_5 + r5_5 - d5_6 + r5_6 + d5_1 - r5_1 = 0;
subject to n35:
	-d5_7 + r5_7 - d5_8 + r5_8 - d5_9 + r5_9 + d5_2 - r5_2 + 2000.0 = 0;
subject to n45:
	-d5_10 + r5_10 - d5_11 + r5_11 - d5_12 + r5_12 + d5_4 - r5_4 = 0;
subject to n55:
	-d5_13 + r5_13 - d5_14 + r5_14 - d5_15 + r5_15 + d5_7 - r5_7 - 2000.0 = 0;
subject to n65:
	-d5_16 + r5_16 + d5_10 - r5_10 + d5_13 - r5_13 = 0;
subject to n75:
	-d5_17 + r5_17 + d5_3 - r5_3 + d5_5 - r5_5 + d5_8 - r5_8 = 0;
subject to n85:
	-d5_18 + r5_18 + d5_6 - r5_6 + d5_9 - r5_9 + d5_11 - r5_11 + d5_14 - r5_14 + 
	d5_17 - r5_17 = 0;
subject to n95:
	d5_12 - r5_12 + d5_15 - r5_15 + d5_16 - r5_16 + d5_18 - r5_18 = 0;
subject to n16:
	-d6_1 + r6_1 - d6_2 + r6_2 - d6_3 + r6_3 = 0;
subject to n26:
	-d6_4 + r6_4 - d6_5 + r6_5 - d6_6 + r6_6 + d6_1 - r6_1 = 0;
subject to n36:
	-d6_7 + r6_7 - d6_8 + r6_8 - d6_9 + r6_9 + d6_2 - r6_2 = 0;
subject to n46:
	-d6_10 + r6_10 - d6_11 + r6_11 - d6_12 + r6_12 + d6_4 - r6_4 + 1000.0 = 0;
subject to n56:
	-d6_13 + r6_13 - d6_14 + r6_14 - d6_15 + r6_15 + d6_7 - r6_7 - 1000.0 = 0;
subject to n66:
	-d6_16 + r6_16 + d6_10 - r6_10 + d6_13 - r6_13 = 0;
subject to n76:
	-d6_17 + r6_17 + d6_3 - r6_3 + d6_5 - r6_5 + d6_8 - r6_8 = 0;
subject to n86:
	-d6_18 + r6_18 + d6_6 - r6_6 + d6_9 - r6_9 + d6_11 - r6_11 + d6_14 - r6_14 + 
	d6_17 - r6_17 = 0;
subject to n96:
	d6_12 - r6_12 + d6_15 - r6_15 + d6_16 - r6_16 + d6_18 - r6_18 = 0;
subject to n17:
	-d7_1 + r7_1 - d7_2 + r7_2 - d7_3 + r7_3 = 0;
subject to n27:
	-d7_4 + r7_4 - d7_5 + r7_5 - d7_6 + r7_6 + d7_1 - r7_1 - 200.0 = 0;
subject to n37:
	-d7_7 + r7_7 - d7_8 + r7_8 - d7_9 + r7_9 + d7_2 - r7_2 + 200.0 = 0;
subject to n47:
	-d7_10 + r7_10 - d7_11 + r7_11 - d7_12 + r7_12 + d7_4 - r7_4 = 0;
subject to n57:
	-d7_13 + r7_13 - d7_14 + r7_14 - d7_15 + r7_15 + d7_7 - r7_7 = 0;
subject to n67:
	-d7_16 + r7_16 + d7_10 - r7_10 + d7_13 - r7_13 = 0;
subject to n77:
	-d7_17 + r7_17 + d7_3 - r7_3 + d7_5 - r7_5 + d7_8 - r7_8 = 0;
subject to n87:
	-d7_18 + r7_18 + d7_6 - r7_6 + d7_9 - r7_9 + d7_11 - r7_11 + d7_14 - r7_14 + 
	d7_17 - r7_17 = 0;
subject to n97:
	d7_12 - r7_12 + d7_15 - r7_15 + d7_16 - r7_16 + d7_18 - r7_18 = 0;
subject to n18:
	-d8_1 + r8_1 - d8_2 + r8_2 - d8_3 + r8_3 = 0;
subject to n28:
	-d8_4 + r8_4 - d8_5 + r8_5 - d8_6 + r8_6 + d8_1 - r8_1 - 200.0 = 0;
subject to n38:
	-d8_7 + r8_7 - d8_8 + r8_8 - d8_9 + r8_9 + d8_2 - r8_2 = 0;
subject to n48:
	-d8_10 + r8_10 - d8_11 + r8_11 - d8_12 + r8_12 + d8_4 - r8_4 + 200.0 = 0;
subject to n58:
	-d8_13 + r8_13 - d8_14 + r8_14 - d8_15 + r8_15 + d8_7 - r8_7 = 0;
subject to n68:
	-d8_16 + r8_16 + d8_10 - r8_10 + d8_13 - r8_13 = 0;
subject to n78:
	-d8_17 + r8_17 + d8_3 - r8_3 + d8_5 - r8_5 + d8_8 - r8_8 = 0;
subject to n88:
	-d8_18 + r8_18 + d8_6 - r8_6 + d8_9 - r8_9 + d8_11 - r8_11 + d8_14 - r8_14 + 
	d8_17 - r8_17 = 0;
subject to n98:
	d8_12 - r8_12 + d8_15 - r8_15 + d8_16 - r8_16 + d8_18 - r8_18 = 0;
subject to n19:
	-d9_1 + r9_1 - d9_2 + r9_2 - d9_3 + r9_3 = 0;
subject to n29:
	-d9_4 + r9_4 - d9_5 + r9_5 - d9_6 + r9_6 + d9_1 - r9_1 - 100.0 = 0;
subject to n39:
	-d9_7 + r9_7 - d9_8 + r9_8 - d9_9 + r9_9 + d9_2 - r9_2 = 0;
subject to n49:
	-d9_10 + r9_10 - d9_11 + r9_11 - d9_12 + r9_12 + d9_4 - r9_4 = 0;
subject to n59:
	-d9_13 + r9_13 - d9_14 + r9_14 - d9_15 + r9_15 + d9_7 - r9_7 + 100.0 = 0;
subject to n69:
	-d9_16 + r9_16 + d9_10 - r9_10 + d9_13 - r9_13 = 0;
subject to n79:
	-d9_17 + r9_17 + d9_3 - r9_3 + d9_5 - r9_5 + d9_8 - r9_8 = 0;
subject to n89:
	-d9_18 + r9_18 + d9_6 - r9_6 + d9_9 - r9_9 + d9_11 - r9_11 + d9_14 - r9_14 + 
	d9_17 - r9_17 = 0;
subject to n99:
	d9_12 - r9_12 + d9_15 - r9_15 + d9_16 - r9_16 + d9_18 - r9_18 = 0;
subject to n110:
	-d10_1 + r10_1 - d10_2 + r10_2 - d10_3 + r10_3 = 0;
subject to n210:
	-d10_4 + r10_4 - d10_5 + r10_5 - d10_6 + r10_6 + d10_1 - r10_1 = 0;
subject to n310:
	-d10_7 + r10_7 - d10_8 + r10_8 - d10_9 + r10_9 + d10_2 - r10_2 - 100.0 = 0;
subject to n410:
	-d10_10 + r10_10 - d10_11 + r10_11 - d10_12 + r10_12 + d10_4 - r10_4 + 100.0 = 
	0;
subject to n510:
	-d10_13 + r10_13 - d10_14 + r10_14 - d10_15 + r10_15 + d10_7 - r10_7 = 0;
subject to n610:
	-d10_16 + r10_16 + d10_10 - r10_10 + d10_13 - r10_13 = 0;
subject to n710:
	-d10_17 + r10_17 + d10_3 - r10_3 + d10_5 - r10_5 + d10_8 - r10_8 = 0;
subject to n810:
	-d10_18 + r10_18 + d10_6 - r10_6 + d10_9 - r10_9 + d10_11 - r10_11 + d10_14 - 
	r10_14 + d10_17 - r10_17 = 0;
subject to n910:
	d10_12 - r10_12 + d10_15 - r10_15 + d10_16 - r10_16 + d10_18 - r10_18 = 0;
subject to n111:
	-d11_1 + r11_1 - d11_2 + r11_2 - d11_3 + r11_3 = 0;
subject to n211:
	-d11_4 + r11_4 - d11_5 + r11_5 - d11_6 + r11_6 + d11_1 - r11_1 = 0;
subject to n311:
	-d11_7 + r11_7 - d11_8 + r11_8 - d11_9 + r11_9 + d11_2 - r11_2 - 200.0 = 0;
subject to n411:
	-d11_10 + r11_10 - d11_11 + r11_11 - d11_12 + r11_12 + d11_4 - r11_4 = 0;
subject to n511:
	-d11_13 + r11_13 - d11_14 + r11_14 - d11_15 + r11_15 + d11_7 - r11_7 + 200.0 = 
	0;
subject to n611:
	-d11_16 + r11_16 + d11_10 - r11_10 + d11_13 - r11_13 = 0;
subject to n711:
	-d11_17 + r11_17 + d11_3 - r11_3 + d11_5 - r11_5 + d11_8 - r11_8 = 0;
subject to n811:
	-d11_18 + r11_18 + d11_6 - r11_6 + d11_9 - r11_9 + d11_11 - r11_11 + d11_14 - 
	r11_14 + d11_17 - r11_17 = 0;
subject to n911:
	d11_12 - r11_12 + d11_15 - r11_15 + d11_16 - r11_16 + d11_18 - r11_18 = 0;
subject to n112:
	-d12_1 + r12_1 - d12_2 + r12_2 - d12_3 + r12_3 = 0;
subject to n212:
	-d12_4 + r12_4 - d12_5 + r12_5 - d12_6 + r12_6 + d12_1 - r12_1 = 0;
subject to n312:
	-d12_7 + r12_7 - d12_8 + r12_8 - d12_9 + r12_9 + d12_2 - r12_2 = 0;
subject to n412:
	-d12_10 + r12_10 - d12_11 + r12_11 - d12_12 + r12_12 + d12_4 - r12_4 - 100.0 = 
	0;
subject to n512:
	-d12_13 + r12_13 - d12_14 + r12_14 - d12_15 + r12_15 + d12_7 - r12_7 + 100.0 = 
	0;
subject to n612:
	-d12_16 + r12_16 + d12_10 - r12_10 + d12_13 - r12_13 = 0;
subject to n712:
	-d12_17 + r12_17 + d12_3 - r12_3 + d12_5 - r12_5 + d12_8 - r12_8 = 0;
subject to n812:
	-d12_18 + r12_18 + d12_6 - r12_6 + d12_9 - r12_9 + d12_11 - r12_11 + d12_14 - 
	r12_14 + d12_17 - r12_17 = 0;
subject to n912:
	d12_12 - r12_12 + d12_15 - r12_15 + d12_16 - r12_16 + d12_18 - r12_18 = 0;

solve;
	display cd1;
	display cr1;
	display d1_1;
	display r1_1;
	display d2_1;
	display r2_1;
	display d3_1;
	display r3_1;
	display d4_1;
	display r4_1;
	display d5_1;
	display r5_1;
	display d6_1;
	display r6_1;
	display d7_1;
	display r7_1;
	display d8_1;
	display r8_1;
	display d9_1;
	display r9_1;
	display d10_1;
	display r10_1;
	display d11_1;
	display r11_1;
	display d12_1;
	display r12_1;
	display cd2;
	display cr2;
	display d1_2;
	display r1_2;
	display d2_2;
	display r2_2;
	display d3_2;
	display r3_2;
	display d4_2;
	display r4_2;
	display d5_2;
	display r5_2;
	display d6_2;
	display r6_2;
	display d7_2;
	display r7_2;
	display d8_2;
	display r8_2;
	display d9_2;
	display r9_2;
	display d10_2;
	display r10_2;
	display d11_2;
	display r11_2;
	display d12_2;
	display r12_2;
	display cd3;
	display cr3;
	display d1_3;
	display r1_3;
	display d2_3;
	display r2_3;
	display d3_3;
	display r3_3;
	display d4_3;
	display r4_3;
	display d5_3;
	display r5_3;
	display d6_3;
	display r6_3;
	display d7_3;
	display r7_3;
	display d8_3;
	display r8_3;
	display d9_3;
	display r9_3;
	display d10_3;
	display r10_3;
	display d11_3;
	display r11_3;
	display d12_3;
	display r12_3;
	display cd4;
	display cr4;
	display d1_4;
	display r1_4;
	display d2_4;
	display r2_4;
	display d3_4;
	display r3_4;
	display d4_4;
	display r4_4;
	display d5_4;
	display r5_4;
	display d6_4;
	display r6_4;
	display d7_4;
	display r7_4;
	display d8_4;
	display r8_4;
	display d9_4;
	display r9_4;
	display d10_4;
	display r10_4;
	display d11_4;
	display r11_4;
	display d12_4;
	display r12_4;
	display cd5;
	display cr5;
	display d1_5;
	display r1_5;
	display d2_5;
	display r2_5;
	display d3_5;
	display r3_5;
	display d4_5;
	display r4_5;
	display d5_5;
	display r5_5;
	display d6_5;
	display r6_5;
	display d7_5;
	display r7_5;
	display d8_5;
	display r8_5;
	display d9_5;
	display r9_5;
	display d10_5;
	display r10_5;
	display d11_5;
	display r11_5;
	display d12_5;
	display r12_5;
	display cd6;
	display cr6;
	display d1_6;
	display r1_6;
	display d2_6;
	display r2_6;
	display d3_6;
	display r3_6;
	display d4_6;
	display r4_6;
	display d5_6;
	display r5_6;
	display d6_6;
	display r6_6;
	display d7_6;
	display r7_6;
	display d8_6;
	display r8_6;
	display d9_6;
	display r9_6;
	display d10_6;
	display r10_6;
	display d11_6;
	display r11_6;
	display d12_6;
	display r12_6;
	display cd7;
	display cr7;
	display d1_7;
	display r1_7;
	display d2_7;
	display r2_7;
	display d3_7;
	display r3_7;
	display d4_7;
	display r4_7;
	display d5_7;
	display r5_7;
	display d6_7;
	display r6_7;
	display d7_7;
	display r7_7;
	display d8_7;
	display r8_7;
	display d9_7;
	display r9_7;
	display d10_7;
	display r10_7;
	display d11_7;
	display r11_7;
	display d12_7;
	display r12_7;
	display cd8;
	display cr8;
	display d1_8;
	display r1_8;
	display d2_8;
	display r2_8;
	display d3_8;
	display r3_8;
	display d4_8;
	display r4_8;
	display d5_8;
	display r5_8;
	display d6_8;
	display r6_8;
	display d7_8;
	display r7_8;
	display d8_8;
	display r8_8;
	display d9_8;
	display r9_8;
	display d10_8;
	display r10_8;
	display d11_8;
	display r11_8;
	display d12_8;
	display r12_8;
	display cd9;
	display cr9;
	display d1_9;
	display r1_9;
	display d2_9;
	display r2_9;
	display d3_9;
	display r3_9;
	display d4_9;
	display r4_9;
	display d5_9;
	display r5_9;
	display d6_9;
	display r6_9;
	display d7_9;
	display r7_9;
	display d8_9;
	display r8_9;
	display d9_9;
	display r9_9;
	display d10_9;
	display r10_9;
	display d11_9;
	display r11_9;
	display d12_9;
	display r12_9;
	display cd10;
	display cr10;
	display d1_10;
	display r1_10;
	display d2_10;
	display r2_10;
	display d3_10;
	display r3_10;
	display d4_10;
	display r4_10;
	display d5_10;
	display r5_10;
	display d6_10;
	display r6_10;
	display d7_10;
	display r7_10;
	display d8_10;
	display r8_10;
	display d9_10;
	display r9_10;
	display d10_10;
	display r10_10;
	display d11_10;
	display r11_10;
	display d12_10;
	display r12_10;
	display cd11;
	display cr11;
	display d1_11;
	display r1_11;
	display d2_11;
	display r2_11;
	display d3_11;
	display r3_11;
	display d4_11;
	display r4_11;
	display d5_11;
	display r5_11;
	display d6_11;
	display r6_11;
	display d7_11;
	display r7_11;
	display d8_11;
	display r8_11;
	display d9_11;
	display r9_11;
	display d10_11;
	display r10_11;
	display d11_11;
	display r11_11;
	display d12_11;
	display r12_11;
	display cd12;
	display cr12;
	display d1_12;
	display r1_12;
	display d2_12;
	display r2_12;
	display d3_12;
	display r3_12;
	display d4_12;
	display r4_12;
	display d5_12;
	display r5_12;
	display d6_12;
	display r6_12;
	display d7_12;
	display r7_12;
	display d8_12;
	display r8_12;
	display d9_12;
	display r9_12;
	display d10_12;
	display r10_12;
	display d11_12;
	display r11_12;
	display d12_12;
	display r12_12;
	display cd13;
	display cr13;
	display d1_13;
	display r1_13;
	display d2_13;
	display r2_13;
	display d3_13;
	display r3_13;
	display d4_13;
	display r4_13;
	display d5_13;
	display r5_13;
	display d6_13;
	display r6_13;
	display d7_13;
	display r7_13;
	display d8_13;
	display r8_13;
	display d9_13;
	display r9_13;
	display d10_13;
	display r10_13;
	display d11_13;
	display r11_13;
	display d12_13;
	display r12_13;
	display cd14;
	display cr14;
	display d1_14;
	display r1_14;
	display d2_14;
	display r2_14;
	display d3_14;
	display r3_14;
	display d4_14;
	display r4_14;
	display d5_14;
	display r5_14;
	display d6_14;
	display r6_14;
	display d7_14;
	display r7_14;
	display d8_14;
	display r8_14;
	display d9_14;
	display r9_14;
	display d10_14;
	display r10_14;
	display d11_14;
	display r11_14;
	display d12_14;
	display r12_14;
	display cd15;
	display cr15;
	display d1_15;
	display r1_15;
	display d2_15;
	display r2_15;
	display d3_15;
	display r3_15;
	display d4_15;
	display r4_15;
	display d5_15;
	display r5_15;
	display d6_15;
	display r6_15;
	display d7_15;
	display r7_15;
	display d8_15;
	display r8_15;
	display d9_15;
	display r9_15;
	display d10_15;
	display r10_15;
	display d11_15;
	display r11_15;
	display d12_15;
	display r12_15;
	display cd16;
	display cr16;
	display d1_16;
	display r1_16;
	display d2_16;
	display r2_16;
	display d3_16;
	display r3_16;
	display d4_16;
	display r4_16;
	display d5_16;
	display r5_16;
	display d6_16;
	display r6_16;
	display d7_16;
	display r7_16;
	display d8_16;
	display r8_16;
	display d9_16;
	display r9_16;
	display d10_16;
	display r10_16;
	display d11_16;
	display r11_16;
	display d12_16;
	display r12_16;
	display cd17;
	display cr17;
	display d1_17;
	display r1_17;
	display d2_17;
	display r2_17;
	display d3_17;
	display r3_17;
	display d4_17;
	display r4_17;
	display d5_17;
	display r5_17;
	display d6_17;
	display r6_17;
	display d7_17;
	display r7_17;
	display d8_17;
	display r8_17;
	display d9_17;
	display r9_17;
	display d10_17;
	display r10_17;
	display d11_17;
	display r11_17;
	display d12_17;
	display r12_17;
	display cd18;
	display cr18;
	display d1_18;
	display r1_18;
	display d2_18;
	display r2_18;
	display d3_18;
	display r3_18;
	display d4_18;
	display r4_18;
	display d5_18;
	display r5_18;
	display d6_18;
	display r6_18;
	display d7_18;
	display r7_18;
	display d8_18;
	display r8_18;
	display d9_18;
	display r9_18;
	display d10_18;
	display r10_18;
	display d11_18;
	display r11_18;
	display d12_18;
	display r12_18;
display obj;
