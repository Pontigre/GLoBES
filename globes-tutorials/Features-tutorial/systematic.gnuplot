#
# $
#
#
set terminal postscript eps enhanced color 
set output 'systematic.eps'
set title "Output of systematic"

set key top left
set style data line
set xrange [0.001:0.1]
set logscale x
set yrange [0:10]
set xlabel "sin^22{/Symbol q}_{13}"
set ylabel "{/Symbol D} {/Symbol c}^2"
set multiplot
plot "sys-data0" using (exp($1*log(10))):($2)\
title "no systematics",\
 "sys-data1" using (exp($1*log(10))):($2)\
title "norm + calibration",\
 "sys-data2" using (exp($1*log(10))):($2)\
title "+ shape",\
 "sys-data3" using (exp($1*log(10))):($2)\
title "+ bin-to-bin"
unset key
plot 2.7
