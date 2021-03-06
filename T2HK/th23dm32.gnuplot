#!/usr/bin/gnuplot

# Redirect output to file
set terminal postscript enhanced eps color
set output "th23dm32.eps"

# Layout
set contour base             # Draw contours on the base plane
unset surface                # Do not draw 3D surface
unset key    
set mxtics
set xrange [0.4:0.6] 
set yrange [0.0022:0.0026]           
set mytics
set view 0,0,1.4             # View 3D plot from above to obtain effective 2d contour plot
set cntrparam cubicspline    # Smooth contours
set cntrparam levels discrete 4.6, 11.83   # Draw contours at C.L. 90% and 3 sigma

# Legend
set title "Confidence regions in the {/Symbol q}_{13}-{/Symbol d}_{CP}"
set label "sin^2(2{/Symbol q}_{23})" at graph 0.5,-0.14 center
set label "{/Symbol D}m^2_{32}" at graph 1.14,0.5 center rotate by -90

# Do the actual plotting
splot "th23dm32.dat" with lines lt -2 lw 2

#set key outside
#set xrange [0.4:0.6] 
#set yrange [0.0022:0.0026]   
#splot "th23dm32.dat" with lines lt -2 lw 2 title "Normal Uncertainties", \
#      "th13delta_05xerror.dat" with lines lt 1 lw 2 title "Half Uncertainties"


# Show resulting EPS figure
system "evince th23dm32.eps"

