#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps monochrome
set output "th23dm31.eps"

# Layout
set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
set xrange [34:56]           # Plot range in the horizontal direction
set yrange [2.3e-3:3.0e-3]   # Plot range in the vertical direction
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam cubicspline    # Smooth contours
set cntrparam levels discrete 2.3, 6.18, 11.83   # Draw contours at 1, 2, 3 sigma

# Legend
set title "Confidence regions in the {/Symbol q}_{23}-{/Symbol D}m_{31}^2 plane" ,-3
set label "{/Symbol q}_{23} [degrees]" at graph 0.5,-0.14 center
set label "{/Symbol D}m_{31}^2 [eV^2]" at graph 1.14,0.5 center rotate by -90
set key graph 0.12, 0.23 spacing 1.3    # Position and vertical spacing of legend

# Do the actual plotting
splot "th23dm31.dat" with lines lt -2 lw 2

# Show resulting EPS figure
system "gv th23dm31.eps"

