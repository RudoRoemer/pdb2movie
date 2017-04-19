# Rigid cluster decomposition coloring script for PyMol
#
# Created by Dan Farrell, Brandon Hespenheide.
# Department of Physics and Astronomy
# Biophysics Theory Group
# Arizona State University
############################################################

from pymol import cmd
from pymol.cgo import *

bg_color white
# script to make a movie of the FRODA generated structures.

from glob import glob
filelist = glob ("lyso_froda_*.pdb") 
filelist.sort()
cmd.load("lyso_RCD.pdb", "lyso_froda")
for file in filelist: cmd.load( file, "lyso_froda") 
show lines, lyso_froda
color black
color 0x0000b2, ( b > 0.99 and b < 1.01)
color 0x57eb0f, ( b > 1.99 and b < 2.01)
color 0xe75ab7, ( b > 2.99 and b < 3.01)
color 0x23999b, ( b > 3.99 and b < 4.01)
color 0xc19b44, ( b > 4.99 and b < 5.01)
# Draw hbonds, hydrophoic tethers and stacked rings as distance objects
set dash_gap, 0.1
distance hbonds = id    18 , id    80
distance hbonds = id    82 , id    27
distance hbonds = id    96 , id    23
distance hbonds = id   113 , id    37
distance hbonds = id   135 , id    56
distance hbonds = id   175 , id    91
distance hbonds = id   188 , id   108
distance hbonds = id   198 , id   127
distance hbonds = id   293 , id   887
distance hbonds = id   316 , id   440
distance hbonds = id   444 , id   307
distance hbonds = id   447 , id   398
distance hbonds = id   459 , id   506
distance hbonds = id   474 , id   223
distance hbonds = id   513 , id   454
distance hbonds = id   517 , id  1113
distance hbonds = id   548 , id   419
distance hbonds = id   601 , id   565
distance hbonds = id   626 , id   600
distance hbonds = id   679 , id   623
distance hbonds = id   693 , id   634
distance hbonds = id   712 , id   653
distance hbonds = id   726 , id   667
distance hbonds = id   740 , id   677
distance hbonds = id   759 , id   687
distance hbonds = id   772 , id   709
distance hbonds = id   782 , id   724
distance hbonds = id   790 , id   720
distance hbonds = id   803 , id   735
distance hbonds = id   818 , id   754
distance hbonds = id   840 , id   981
distance hbonds = id   840 , id   982
distance hbonds = id   915 , id   247
distance hbonds = id   960 , id  1007
distance hbonds = id   983 , id   932
distance hbonds = id   994 , id   929
distance hbonds = id  1008 , id   943
distance hbonds = id  1023 , id   965
distance hbonds = id  1044 , id   977
distance hbonds = id  1066 , id   992
distance hbonds = id  1083 , id  1002
distance hbonds = id  1098 , id  1017
distance hbonds = id  1114 , id  1039
distance hbonds = id  1125 , id  1058
distance hbonds = id  1142 , id  1078
distance hbonds = id  1151 , id  1092
distance hbonds = id  1161 , id  1109
distance hbonds = id  1173 , id  1121
distance hbonds = id  1193 , id  1137
distance hbonds = id  1210 , id  1149
distance hbonds = id  1221 , id  1159
distance hbonds = id  1240 , id  1169
distance hbonds = id  1262 , id  1185
distance hbonds = id  1349 , id  1278
distance hbonds = id  1404 , id  1324
distance hbonds = id  1421 , id  1343
distance hbonds = id  1504 , id  1461
distance hbonds = id  1513 , id  1462
distance hbonds = id  1528 , id  1458
distance hbonds = id  1536 , id  1395
distance hbonds = id  1570 , id  1496
distance hbonds = id  1589 , id  1520
distance hbonds = id  1608 , id  1544
distance hbonds = id  1612 , id    94
distance hbonds = id  1613 , id   174
distance hbonds = id  1622 , id  1555
distance hbonds = id  1638 , id  1565
distance hbonds = id  1658 , id  1584
distance hbonds = id  1676 , id  1603
distance hbonds = id  1692 , id  1617
distance hbonds = id  1705 , id  1634
distance hbonds = id  1751 , id  1704
distance hbonds = id  1765 , id  1711
distance hbonds = id  1788 , id  1747
distance hbonds = id  1844 , id  1780
distance hbonds = id  1866 , id  1800
distance hbonds = id  1887 , id  1814
distance hbonds = id  1904 , id  1828
distance hbonds = id  1924 , id  1839
distance hbonds = id  1968 , id  1428
distance hbonds = id  1982 , id  1882
distance hbonds = id  1993 , id  2038
distance hbonds = id  2040 , id  2039
distance hbonds = id  2051 , id  1974
distance hbonds = id  2061 , id  1998
distance hbonds = id  2073 , id  2022
distance hbonds = id  2090 , id  2034
distance hbonds = id  2104 , id  2049
distance hbonds = id  2120 , id  2059
distance hbonds = id  2134 , id  2085
distance hbonds = id  2157 , id  2099
distance hbonds = id  2179 , id   381
distance hbonds = id  2218 , id  2152
distance hbonds = id  2235 , id  2150
distance hbonds = id  2250 , id  2161
distance hbonds = id  2265 , id  2185
distance hbonds = id  2321 , id   186
distance hbonds = id  2322 , id  1607
distance hbonds = id  2323 , id   186
distance hbonds = id  2329 , id  2261
distance hbonds = id  2343 , id  2275
distance hbonds = id  2367 , id  2289
distance hbonds = id  2375 , id   173
distance hbonds = id  2376 , id   173
distance hbonds = id  2387 , id  2303
distance hbonds = id  2404 , id  2327
distance hbonds = id  2422 , id  2337
distance hbonds = id  2425 , id  2337
distance hbonds = id  2439 , id  2359
distance hbonds = id  2454 , id  2383
distance hbonds = id  2474 , id  2399
distance hbonds = id  2494 , id  2418
distance hbonds = id  2497 , id  2418
distance hbonds = id  2505 , id  2432
distance hbonds = id  2515 , id  2493
distance hbonds = id  2518 , id  2493
distance hbonds = id  2554 , id  2514
distance hbonds = id  2580 , id  2525
distance hbonds = id  2588 , id   174
color red, hbonds
hide labels, hbonds
disable hbonds
distance hydrophobics = id     5 , id  2529
distance hydrophobics = id     6 , id  2531
distance hydrophobics = id     7 , id  2572
distance hydrophobics = id     8 , id    77
distance hydrophobics = id    40 , id    59
distance hydrophobics = id    41 , id  1585
distance hydrophobics = id    41 , id  1545
distance hydrophobics = id    60 , id  1063
distance hydrophobics = id    62 , id  1124
distance hydrophobics = id    63 , id  1123
distance hydrophobics = id    93 , id  2576
distance hydrophobics = id    94 , id  1545
distance hydrophobics = id    94 , id  2534
distance hydrophobics = id    95 , id  2577
distance hydrophobics = id   109 , id  1065
distance hydrophobics = id   111 , id  1652
distance hydrophobics = id   112 , id  1653
distance hydrophobics = id   183 , id   483
distance hydrophobics = id   207 , id   482
distance hydrophobics = id   207 , id   993
distance hydrophobics = id   224 , id   311
distance hydrophobics = id   225 , id   313
distance hydrophobics = id   290 , id   458
distance hydrophobics = id   291 , id   423
distance hydrophobics = id   292 , id   689
distance hydrophobics = id   364 , id  2263
distance hydrophobics = id   403 , id   579
distance hydrophobics = id   420 , id   678
distance hydrophobics = id   456 , id   914
distance hydrophobics = id   457 , id   738
distance hydrophobics = id   458 , id   739
distance hydrophobics = id   481 , id  1062
distance hydrophobics = id   483 , id  1064
distance hydrophobics = id   484 , id  1061
distance hydrophobics = id   637 , id   690
distance hydrophobics = id   738 , id   800
distance hydrophobics = id   738 , id   913
distance hydrophobics = id   739 , id   863
distance hydrophobics = id   800 , id   913
distance hydrophobics = id   801 , id   864
distance hydrophobics = id   801 , id   826
distance hydrophobics = id   826 , id   864
distance hydrophobics = id   911 , id   979
distance hydrophobics = id  1064 , id  1124
distance hydrophobics = id  1064 , id  1657
distance hydrophobics = id  1110 , id  1656
distance hydrophobics = id  1171 , id  1238
distance hydrophobics = id  1172 , id  1588
distance hydrophobics = id  1218 , id  1637
distance hydrophobics = id  1219 , id  1396
distance hydrophobics = id  1220 , id  1586
distance hydrophobics = id  1401 , id  1521
distance hydrophobics = id  1401 , id  1588
distance hydrophobics = id  1443 , id  2008
distance hydrophobics = id  1443 , id  1497
distance hydrophobics = id  1482 , id  2434
distance hydrophobics = id  1545 , id  2535
distance hydrophobics = id  1546 , id  2533
distance hydrophobics = id  1556 , id  2452
distance hydrophobics = id  1618 , id  1750
distance hydrophobics = id  1620 , id  2451
distance hydrophobics = id  1621 , id  2103
distance hydrophobics = id  1621 , id  2194
distance hydrophobics = id  1689 , id  1750
distance hydrophobics = id  1690 , id  1786
distance hydrophobics = id  1691 , id  2187
distance hydrophobics = id  1783 , id  2103
distance hydrophobics = id  1883 , id  2050
distance hydrophobics = id  1885 , id  2035
distance hydrophobics = id  1886 , id  2086
distance hydrophobics = id  1902 , id  2452
distance hydrophobics = id  1903 , id  2007
distance hydrophobics = id  1903 , id  2050
distance hydrophobics = id  1903 , id  2450
distance hydrophobics = id  1936 , id  1976
distance hydrophobics = id  1953 , id  2006
distance hydrophobics = id  2060 , id  2468
distance hydrophobics = id  2060 , id  2402
distance hydrophobics = id  2100 , id  2403
distance hydrophobics = id  2102 , id  2449
distance hydrophobics = id  2103 , id  2194
distance hydrophobics = id  2119 , id  2212
distance hydrophobics = id  2193 , id  2384
distance hydrophobics = id  2193 , id  2328
distance hydrophobics = id  2194 , id  2403
distance hydrophobics = id  2195 , id  2401
distance hydrophobics = id  2212 , id  2403
distance hydrophobics = id  2360 , id  2577
distance hydrophobics = id  2385 , id  2451
distance hydrophobics = id  2420 , id  2562
color green, hydrophobics
hide labels, hydrophobics
disable hydrophobics
# Rigid Cluster 1 has 658 atoms.
create RC1, ( b > 0.99 and b < 1.01)
show sticks, RC1
set line_width = 3, RC1
color 0x0000b2, RC1

# Rigid Cluster 2 has 131 atoms.
create RC2, ( b > 1.99 and b < 2.01)
show sticks, RC2
set line_width = 3, RC2
color 0x57eb0f, RC2

# Rigid Cluster 3 has 76 atoms.
create RC3, ( b > 2.99 and b < 3.01)
show sticks, RC3
set line_width = 3, RC3
color 0xe75ab7, RC3

# Rigid Cluster 4 has 45 atoms.
create RC4, ( b > 3.99 and b < 4.01)
show sticks, RC4
set line_width = 3, RC4
color 0x23999b, RC4

# Rigid Cluster 5 has 25 atoms.
create RC5, ( b > 4.99 and b < 5.01)
show sticks, RC5
set line_width = 3, RC5
color 0xc19b44, RC5

# Rigid Cluster BIN2
create BIN2, ( b > 5.99 and b < 12.01)
show sticks, BIN2
set line_width = 3, BIN2
color gray, BIN2
disable BIN2

# Rigid Cluster BIN1
create BIN1, ( b > 12.99 and b < 549.01)
show sticks, BIN1
set line_width = 3, BIN1
color gray, BIN1
disable BIN1

