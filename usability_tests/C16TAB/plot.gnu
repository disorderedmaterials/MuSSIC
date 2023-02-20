
set encoding iso_8859_1
set terminal pngcairo  enhanced color lw 2 size 500,500 font 'Arial,15'


# The file we'll write to
set logscale x

set logscale y
set yrange [0.7:20.0]
set output 'HD-CG-C16TAB.png'
# The graphic title
#set title 'HD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F_{CG}(Q) (barn sr^{-1} atom^{-1})"
#plot the graphic
plot   'HD-CG-weighted-total.sq' using 1:($2*0.1+1.0)  with linespoints  title "   3 carbon atoms per CG bead "  lw 2 linecolor 'blue' pointtype 7 pointsize 1.0, \
       'SANS-2d.dat' using 1:($2+1.0)  with linespoints  title " SANS-2d "  lw 2 linecolor 'red' pointtype 7 pointsize 1.0

