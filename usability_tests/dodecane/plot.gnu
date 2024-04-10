
set encoding iso_8859_1
set terminal pngcairo  enhanced color lw 1 size 500,500 font 'Arial,15'


# The file we'll write to
set logscale x
#set logscale y
set output 'HD-CG-dodecane-3C-unbound.png'
# The graphic title
#set title 'HD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F^{inter}_{CG}(Q) (barn sr^{-1} atom^{-1})"
#plot the graphic
plot   'HD-CG-weighted-total.sq' using 1:($4+1.0)  with linespoints  title "3C  per CG bead - F1_{CG}(Q)"  lw 1 linecolor 'blue' pointtype 7 pointsize 1.0, \
       'HD-weighted-total.sq' using 1:($4+1.0)  with linespoints  title " atomistic "  lw 1 linecolor 'red' pointtype 7 pointsize 1.0

set encoding iso_8859_1
set terminal pngcairo  enhanced color lw 1 size 500,500 font 'Arial,15'


# The file we'll write to
set logscale x
set logscale y
set output 'HD-CG-dodecane-3C.png'
# The graphic title
#set title 'HD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F_{CG}(Q) (barn sr^{-1} atom^{-1})"
#plot the graphic
plot   'HD-CG-weighted-total.sq' using 1:($2+1.0)  with linespoints  title "3C  per CG bead - F1_{CG}(Q)"  lw 1 linecolor 'blue' pointtype 7 pointsize 1.0, \
       'HD-weighted-total.sq' using 1:($2+1.0)  with linespoints  title " atomistic "  lw 1 linecolor 'red' pointtype 7 pointsize 1.0



set encoding iso_8859_1
set terminal pngcairo  enhanced color lw 1 size 500,500 font 'Arial,15'


# The file we'll write to
set logscale x
set logscale y
set output 'HD-CG-dodecane-3C-bound.png'
# The graphic title
#set title 'HD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F^{intra}_{CG}(Q) (barn sr^{-1} atom^{-1})"
#plot the graphic
plot   'HD-CG-weighted-total.sq' using 1:($3+$5+1.0)  with linespoints  title "3C  per CG bead - F1_{CG}(Q)"  lw 1 linecolor 'blue' pointtype 7 pointsize 1.0, \
       'HD-weighted-total.sq' using 1:($3+1.0)  with linespoints  title " atomistic "  lw 1 linecolor 'red' pointtype 7 pointsize 1.0
