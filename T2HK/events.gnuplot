#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "rates.eps"

# Layout used by all
unset key
set cntrparam cubicspline 
set mxtics
set mytics
set xlabel "Energy [GeV]"
set ylabel "Flux"
set xrange [0.0:60]
set yrange [0.0:600.]

plot "rates_e.dat" with linespoints ls 1
set yrange [0.0:7000.]
plot "rates_mu.dat" with linespoints ls 1

# Show resulting EPS figure
system "evince rates.eps"

