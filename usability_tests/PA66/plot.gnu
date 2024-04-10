
set encoding iso_8859_1
set terminal pngcairo  enhanced color lw 2 size 500,500 font 'Arial,15'


# The file we'll write to
set logscale x

set logscale y
set yrange [0.7:15]
set output 'HD-CG-PA66-4C.png'
# The graphic title
#set title 'HD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q) (barn sr^{-1} atom^{-1})"
#plot the graphic
plot   'HD-CG-weighted-total.sq' using 1:($2+1.0)  with linespoints  title "\\~ 4 carbon atoms per CG bead "  lw 2 linecolor 'blue' pointtype 7 pointsize 1.0, \
       'HD-weighted-total.sq' using 1:($2+1.0)  with linespoints  title " atomistic "  lw 2 linecolor 'red' pointtype 7 pointsize 1.0

