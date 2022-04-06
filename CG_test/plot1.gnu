set encoding iso_8859_1
set terminal png size 500,500


# The file we'll write to
set output 'HH-totalsq-compare.png'
# The graphic title
set title 'HH-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
#plot the graphic
plot  "<(sed -n '1,200p' 'HH-weighted-total.sq')" using 1:2  with lines title " Bindu-atomistic" lw 3, \
      "<(sed -n '1,600p' 'Lupe_PA66/06_FQ_atom.dat')" using 1:2 with lines title " Lupe-atomistic" lw 3
# The file we'll write to
set output 'DD-totalsq-compare.png'
# The graphic title
set title 'DD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
#plot the graphic
plot  "<(sed -n '1,200p' 'DD-weighted-total.sq')" using 1:2  with lines title " Bindu-atomistic" lw 3, \
      "<(sed -n '1,600p' 'Lupe_PA66/06_FQ_atom.dat')" using 1:4 with lines title " Lupe-atomistic" lw 3
# The file we'll write to
set output 'DH-totalsq-compare.png'
# The graphic title
set title 'DH-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
#plot the graphic
plot  "<(sed -n '1,200p' 'DH-weighted-total.sq')" using 1:2  with lines title " Bindu-atomistic" lw 3, \
      "<(sed -n '1,600p' 'Lupe_PA66/06_FQ_atom.dat')" using 1:3 with lines title " Lupe-atomistic" lw 3
# The file we'll write to
set output 'DH-CG-compare.png'
# The graphic title
set title 'DH-CG-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
#plot the graphic
plot  "<(sed -n '1,100p' 'CG_test/HH-weighted-total.sq')" using 1:2  with lines title " Bindu-CG" lw 3, \
      "<(sed -n '1,100p' 'DH-weighted-total.sq')" using 1:2 with lines title " DD-atomistic" lw 3, \
      "<(sed -n '1,600p' 'Lupe_PA66/06_FQ_cg_HD.dat')" using 1:2 with lines title " Lupe-CG" lw 3
