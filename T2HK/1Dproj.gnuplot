#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "1Dproj.eps"

# Layout
unset key
set mytics
set mxtics

# Theta 13
set xrange [0.01:0.2]
set yrange [0.0:20.0]
set title "Projection onto the {/Symbol q}_{13}-plane"
set xlabel "sin^2(2{/Symbol q}_{13})"
set ylabel "{/Symbol c}^2"

#plot "th13.dat" with lines

# Theta 23
set xrange [0.22:0.24]
set yrange [0.0:20.0]
set title "Projection onto the {/Symbol q}_{23}-plane"
set xlabel "sin^2({/Symbol q}_{23})"
set ylabel "{/Symbol c}^2"

plot "th23.dat" with lines
system "evince 1Dproj.eps"