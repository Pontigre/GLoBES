#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "main.eps"

# Layout used by all
unset key
set cntrparam cubicspline 
set mxtics
set mytics

set xlabel "Energy [GeV]"
set ylabel "Flux"

#Initial Fluxes
set xrange [0.0:3.0]

set title "Initial Electron Flux - Inverted Hierarchy"
set yrange [0.0:11.]
plot "init_flux_e.dat" with linespoints ls 1

set title "Initial Muon Flux - Inverted Hierarchy"
set yrange [0.0:5000.]
plot "init_flux_mu.dat" with linespoints ls 1

#Event Rates
set xrange [0:60]

set title "Electron CC Event Rates - Inverted Hierarchy"
set yrange [0.0:400.]
plot "events_e_inv.dat" with linespoints ls 1

set title "Muon CC Event Rates - Inverted Hierarchy"
set yrange [0.0:7000.]
plot "events_mu.dat" with linespoints ls 1

#Delta significance
set title "Significance to exclude sin({/Symbol d}) = 0.0 - Inverted Hierarchy"
set xrange [-180.0:180.0]
set yrange [0.0:10.0]
set xlabel "{/Symbol d}_{CP} [{/Symbol p]"
set ylabel "{/Symbol s} = sqrt{{/Symbol c}^2}"
plot "deltasig.dat" with linespoints ls 1


#th13-deltacp scans
set contour base 
unset surface
set xrange [0.01:0.2]
set yrange [-180.0:180.0] 

set label "sin^2(2{/Symbol q}_{13})" at graph 0.5,-0.14 center
set label "{/Symbol d}_{CP} [degree]" at graph 1.14,0.5 center rotate by -90
set view 0,0,1.4 
set cntrparam levels discrete 4.6, 11.83   # Draw contours at C.L. 90% and 3 sigma

set title "Confidence regions in the {/Symbol q}_{13}-{/Symbol d}_{CP} - Inverted Hierarchy"
splot "th13delta.dat" with lines lt -2 lw 2
set title "Confidence regions in the {/Symbol q}_{13}-{/Symbol d}_{CP} - Half Errors - Inverted Hierarchy"
splot "th13delta_05xerror.dat" with lines lt -2 lw 2

# Show resulting EPS figure
system "evince main.eps"

