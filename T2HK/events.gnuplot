#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "rates.eps"

# Layout used by all
unset key
set cntrparam cubicspline 
set mxtics
set mytics
set ylabel "Number of Events"
set xlabel "Energy [GeV]"

set xrange [0.4:1.2]
set yrange [0.0:0.03]
set title "nuPRISM Electron events"
plot "nuPRISM_events_e.dat" with linespoints ls 1
set yrange [0.0:4]
set title "nuPRISM Muon events"
plot "nuPRISM_events_mu.dat" with linespoints ls 1


set xrange [0.0:10]
set title "nuPRISM Initial Electron Flux"
set yrange [0.0:200000]
plot "nuPRISM_flux_e.dat" with linespoints ls 1
set title "nuPRISM Initial Muon Flux"
set yrange [0.0:1e8]
plot "nuPRISM_flux_mu.dat" with linespoints ls 1

set xrange [0.4:1.2]
set yrange [0.0:3000]
set title "T2HK Electron events"
plot "T2HK_events_e.dat" with linespoints ls 1
set title "T2HK Muon events"
set yrange [0.0:50000]
plot "T2HK_events_mu.dat" with linespoints ls 1

set xrange [0.0:10]
set title "T2HK Initial Electron Flux"
set yrange [0.0:20]
plot "T2HK_flux_e.dat" with linespoints ls 1
set title "T2HK Initial Muon Flux"
set yrange [0.0:8000]
plot "T2HK_flux_mu.dat" with linespoints ls 1

# Show resulting EPS figure
system "evince rates.eps"

