set terminal pngcairo size 2400,1000 enhanced font "Sans,12"
set xlabel "x"
set ylabel "y"
set zlabel "f(x,y)"
set key top left
set ticslevel 0
set view 60,35,1,1
set hidden3d

set palette defined (0 "#f7fbff", 0.5 "#6baed6", 1 "#08306b")
set colorbox

do for [caseName in "none log log10 loginv mixed"] {
    dataFile = sprintf("test_%s.dat", caseName)
    outputFile = sprintf("test_%s.png", caseName)

    set output outputFile
    set multiplot layout 1,2 title sprintf("biQuadraticInterpolation: %s", caseName)

    set title "Exact and interpolated surfaces"
    set xlabel "x"
    set ylabel "y"
    set zlabel "f(x,y)"
    set view 60,35,1,1
    set hidden3d
    splot dataFile every 3:3 using 1:2:3 with lines lc rgb "#34495e" lw 1.5 title "exact", \
          dataFile every 3:3 using 1:2:4 with lines lc rgb "#e4572e" lw 1.5 title "interpolated"

    set title "Absolute interpolation error"
    set xlabel "x"
    set ylabel "y"
    unset zlabel
    unset key
    set view map
    unset hidden3d
    set pm3d map
    set cblabel "|interpolated - exact|"
    splot dataFile using 1:2:6 with pm3d

    unset pm3d
    unset multiplot
    unset output
}
