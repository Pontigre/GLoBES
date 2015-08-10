#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "th13delta_nuPRISM.eps"

# Layout
set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
unset key
set xrange [0.01:0.2]            # Plot range in the horizontal direction
set mxtics
set yrange [-180.0:180.0]           # Plot range in the vertical direction
set mytics
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam cubicspline    # Smooth contours
set cntrparam levels discrete 4.6, 11.83   # Draw contours at C.L. 90% and 3 sigma

# Legend
set title "Confidence regions in the {/Symbol q}_{13}-{/Symbol d}_{CP}"
set label "sin^2(2{/Symbol q}_{13})" at graph 0.5,-0.14 center
set label "{/Symbol d}_{CP} [degree]" at graph 1.14,0.5 center rotate by -90

# Do the actual plotting
#splot "th13delta_nuPRISM.dat" with lines lt -2 lw 2

set key 
set cntrparam levels discrete 4.6
splot "th13delta_nuPRISM.dat" with lines lt 1 lw 2 title "nuPRISM constraint", \
"th13delta.dat" with lines lt 0 lw 2 title "normal T2HK"

# Show resulting EPS figure
system "evince th13delta_nuPRISM.eps"

