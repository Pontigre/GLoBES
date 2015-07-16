#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "test_all.eps"

# Layout 
unset key
set mxtics
set mytics

# Test 2
set title "Projection onto sin^2(2{/Symbol q}_{13}) axis"
set xrange [-2.0:0.0]
set xlabel "log(sin^2(2{/Symbol q}_{13}))"
set yrange [0.0:20.] 
set ylabel "{/Symbol c}^2"
plot "test2.dat" using 1:2

set title "Minimize over all but {/Symbol q}_{13}"
plot "test2.dat" using 1:3

# Test 3
set title "Test 3"
set xrange [-4.0:-2.0]
set xlabel "Test3-x"
set yrange [0.:200.] 
set ylabel "Test3-y"

# Test 4a
set title "Test 2"
set xrange [-7.0:-2.0]
set xlabel "Test2-x"
set yrange [4500.:7500.] 
set ylabel "Test2-y"

# Test 4b
set title "Test 2"
set xrange [-7.0:-2.0]
set xlabel "Test2-x"
set yrange [1000.:1800.] 
set ylabel "Test2-y"

# Test 4c
set title "Test 2"
set xrange [-7.0:-2.0]
set xlabel "Test2-x"
set yrange [350.:750.] 
set ylabel "Test2-y"

# Test 4d
set title "Test 2"
set xrange [-7.0:-2.0]
set xlabel "Test2-x"
set yrange [350.:750.] 
set ylabel "Test2-y"

# Show resulting EPS figure
system "evince test_all.eps"

