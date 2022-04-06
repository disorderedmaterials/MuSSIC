set encoding iso_8859_1
set terminal png size 500,500


# The file we'll write to
set logscale x
set logscale y
set output 'HH-CG-compare.png'
# The graphic title
set title 'HH-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
#plot the graphic
plot  '../HH-weighted-total.sq' using 1:($2+1.0) with lines title " atomistic" lw 3 linecolor rgb "blue", \
      '../norm1/HH-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      'HH-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 1C per bead" lw 3 linecolor rgb "red"

set encoding iso_8859_1
set terminal png size 500,500


# The file we'll write to
set logscale x

set logscale y
set output 'DD-CG-compare.png'
# The graphic title
set title 'DD-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
#plot the graphic
plot   'DD-weighted-total.sq' using 1:($2+1.0)  with lines title "CG- 1C per bead" lw 3 linecolor rgb "red", \
      '../norm1/DD-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      '../DD-weighted-total.sq' using 1:($2+1.0) with lines title " atomistic" lw 3 linecolor rgb "blue"

set encoding iso_8859_1
set terminal png size 500,500


# The file we'll write to
set logscale x
set logscale y
set output 'DH-CG-compare.png'
# The graphic title
set title 'DH-Weighted total S(Q)' font ",15"
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
set xtics font ",15"
set ytics font ",15"
set key font ",15"

#plot the graphic

plot  'DH-weighted-total.sq' using 1:($2+1.0) with lines title "CG-1C per bead" lw 3 linecolor rgb "red", \
      '../norm1/DH-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      '../DH-weighted-total.sq' using 1:($2+1.0) with lines title "atomistic" lw 3 linecolor rgb "blue"


set encoding iso_8859_1
set terminal png size 500,500
# The file we'll write to
set logscale x
set logscale y
set output '50H-CG-compare.png'
# The graphic title
set title '50H-Weighted total S(Q)'
set xlabel "Q(1/\305)"
set ylabel " F(Q)"
set xtics font ",15"
set ytics font ",15"
set key font ",15"

#plot the graphic

plot  '../50H-weighted-total.sq' using 1:($2+1.0) with lines title "atomistic" lw 3 linecolor rgb "blue", \
      '../norm1/50H-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      '50H-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 1C per bead" lw 3 linecolor rgb "red"

set encoding iso_8859_1
set terminal png size 500,500
# The file we'll write to
set logscale x
set logscale y
set output 'HD-CG-compare.png'
# The graphic title
set title 'HD-Weighted total S(Q)' 
set xlabel "Q(\305^{-1})"
set ylabel " F(Q)"
set xtics font ",15"
set ytics font ",15"
set key font ",15"

#plot the graphic

plot  'HD-weighted-total.sq' using 1:($2+1.0)  with lines title "CG- 1C per bead " lw 3 linecolor rgb "red", \
      '../HD-weighted-total.sq' using 1:($2+1.0) with lines title " atomistic" lw 3 linecolor rgb "blue", \
      '../norm1/HD-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \




set encoding iso_8859_1
set terminal png size 500,500
# The file we'll write to
set logscale x
set logscale y
set output '75H-CG-compare.png'
# The graphic title
set title '75H-Weighted total S(Q)' 
set xlabel "Q(\305^{-1})"
set ylabel " F(Q)"
set xtics font ",15"
set ytics font ",15"
set key font ",15"

#plot the graphic

plot  '75H-weighted-total.sq' using 1:($2+1.0)  with lines title "CG- 1C per bead " lw 3 linecolor rgb "red", \
      '../norm1/75H-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      '../75H-weighted-total.sq' using 1:($2+1.0) with lines title " atomistic" lw 3 linecolor rgb "blue"


set encoding iso_8859_1
set terminal png size 500,500
# The file we'll write to
set logscale x
set logscale y
set output '25H-CG-compare.png'
# The graphic title
set title '25H-Weighted total S(Q)' 
set xlabel "Q(\305^{-1})"
set ylabel " F(Q)"
set xtics font ",15"
set ytics font ",15"
set key font ",15"

#plot the graphic

plot  '25D-weighted-total.sq' using 1:($2+1.0)  with lines title "CG- 1C per bead " lw 3 linecolor rgb "red", \
      '../norm1/25H-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      '../25H-weighted-total.sq' using 1:($2+1.0) with lines title " atomistic" lw 3 linecolor rgb "blue"


set encoding iso_8859_1
set terminal png size 500,500
# The file we'll write to
set logscale x
set logscale y
set output 'H25-CG-compare.png'
# The graphic title
set title 'H25-Weighted total S(Q)' 
set xlabel "Q(\305^{-1})"
set ylabel " F(Q)"
set xtics font ",15"
set ytics font ",15"
set key font ",15"

#plot the graphic

plot  'H25-weighted-total.sq' using 1:($2+1.0)  with lines title "CG- 1C per bead " lw 3 linecolor rgb "red", \
      '../norm1/H25-weighted-total.sq' using 1:($2+1.0) with lines title " CG- 3C per bead" lw 3 linecolor rgb "green", \
      '../H25-weighted-total.sq' using 1:($2+1.0) with lines title " atomistic" lw 3 linecolor rgb "blue"
