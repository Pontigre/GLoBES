#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color

# Layout
set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface

set output 'oscillation.eps'
set title "Output of oscillation", -4
unset surface
unset key
set size square
set style data line

set xrange [0:0.05]            # Plot range in the horizontal direction
set yrange [0:10]           # Plot range in the vertical direction
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam cubicspline    # Smooth contours
set cntrparam levels discrete 2.3, 6.18, 11.83   # Draw contours at C.L. 1, 2 and 3 sigma

# Legend
set label "sin^22{/Symbol q}_{13}" at graph 0.5,-0.14 center
set label  "{/Symbol s}_E [MeV]" at graph 1.14,0.5 center rotate by -90
set key graph 0.35, 0.25 spacing 1.3   # Position and vertical spacing of legend
set noclabel # do not label individual contours
# Do the actual plotting
splot "osc-data2" using ($1):(1000*$2):($3) title "with correlations" ,\
 "osc-data1" using ($1):(1000*$2):($3) title "w/o correlations" 

