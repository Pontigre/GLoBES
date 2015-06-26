#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "deltasig.eps"

# Layout
unset key
set cntrparam cubicspline
set xrange [-1.0:1.0]
set mxtics
set yrange [0.0:10.0]
set mytics

# Legend
set title "Significance to exculde {/Symbol d}_{CP} = 0"
set xlabel "{/Symbol d}_{CP} [{/Symbol p]"
set ylabel "{/Symbol s} = sqrt{{/Symbol c}^2}"

# Do the actual plotting
plot "deltasig.dat" with linespoints ls 1

# Show resulting EPS figure
system "evince deltasig.eps"

