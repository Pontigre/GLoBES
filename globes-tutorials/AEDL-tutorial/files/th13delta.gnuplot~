#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "th13delta.eps"

# Layout
set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
set xrange [8:11]            # Plot range in the horizontal direction
set yrange [70:110]           # Plot range in the vertical direction
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam cubicspline    # Smooth contours
set cntrparam levels discrete 4.6, 11.83   # Draw contours at C.L. 90% and 3 sigma

# Legend
set title "Confidence regions in the {/Symbol q}_{13}-{/Symbol d}_{CP} plane" ,-5
set label "{/Symbol q}_{13} [degrees]" at graph 0.5,-0.14 center
set label "{/Symbol d}_{CP} [degrees]" at graph 1.14,0.5 center rotate by -90
set key graph 0.75, 0.25 spacing 1.3   # Position and vertical spacing of legend

# Do the actual plotting
splot "th13delta.dat" with lines lt -2 lw 2

# Show resulting EPS figure
system "gv th13delta.eps &"

