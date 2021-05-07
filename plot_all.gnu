set encoding iso_8859_1
set terminal png size 500,500


# The file we'll write to
set output 'partialgr_OW-OW.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{OW-OW}(r)"
#plot the graphic
plot  "<(sed -n '1,591p' unweighted.gr)" using 7:4  with lines title "OW-OW scattcg" lw 3, \
      'Mix-unweighted-OW-OW.rdf' using 1:2 with lines title "OW-OW Dissolve" lw 3

set output 'partialgr_OW-HW.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{OW-HW}(r)"
#plot the graphic
plot  "<(sed -n '3553,4143p' unweighted.gr)" using 7:4  with lines title "OW-HW scattcg" lw 3, \
      'Mix-unweighted-OW-HW.rdf' using 1:2 with lines title "OW-HW Dissolve" lw 3

set output 'partialgr_HW-HW.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HW-HW}(r)"
#plot the graphic
plot  "<(sed -n '4145,4735p' unweighted.gr)" using 7:4  with lines title "HW-HW scattcg" lw 3, \
      'Mix-unweighted-HW-HW.rdf' using 1:2 with lines title "HW-HW Dissolve" lw 3

# The file we'll write to
set output 'partialgr_OW-CT.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{OW-CT}(r)"
#plot the graphic
plot  "<(sed -n '7105,7695p' unweighted.gr)" using 7:4  with lines title "OW-CT scattcg" lw 3, \
      'Mix-unweighted-OW-CT.rdf' using 1:2 with lines title "OW-CT Dissolve" lw 3

set output 'partialgr_HW-CT.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HW-CT}(r)"
#plot the graphic
plot  "<(sed -n '7697,8287p' unweighted.gr)" using 7:4  with lines title "HW-CT scattcg" lw 3, \
      'Mix-unweighted-HW-CT.rdf' using 1:2 with lines title "HW-CT Dissolve" lw 3

set output 'partialgr_CT-CT.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{CT-CT}(r)"
#plot the graphic
plot  "<(sed -n '8289,8879p' unweighted.gr)" using 7:4  with lines title "CT-CT scattcg" lw 3, \
      'Mix-unweighted-CT-CT.rdf' using 1:2 with lines title "CT-CT Dissolve" lw 3

# The file we'll write to
set output 'partialgr_OW-HC.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{OW-HC}(r)"
#plot the graphic
plot  "<(sed -n '10657,11247p' unweighted.gr)" using 7:4  with lines title "OW-HC scattcg" lw 3, \
      'Mix-unweighted-OW-HC.rdf' using 1:2 with lines title "OW-HC Dissolve" lw 3

set output 'partialgr_HW-HC.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HW-HC}(r)"
#plot the graphic
plot  "<(sed -n '11249,11839p' unweighted.gr)" using 7:4  with lines title "HW-HC scattcg" lw 3, \
      'Mix-unweighted-HW-HC.rdf' using 1:2 with lines title "HW-HC Dissolve" lw 3

set output 'partialgr_CT-HC.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{CT-HC}(r)"
#plot the graphic
plot  "<(sed -n '11841,12431p' unweighted.gr)" using 7:4  with lines title "CT-HC scattcg" lw 3, \
      'Mix-unweighted-CT-HC.rdf' using 1:2 with lines title "CT-HC Dissolve" lw 3

# The file we'll write to
set output 'partialgr_HC-HC.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{HC-HC}(r)"
#plot the graphic
plot  "<(sed -n '12433,13023p' unweighted.gr)" using 7:4  with lines title "HC-HC scattcg" lw 3, \
      'Mix-unweighted-HC-HC.rdf' using 1:2 with lines title "HC-HC Dissolve" lw 3

set output 'partialgr_OW-OH.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{OW-OH}(r)"
#plot the graphic
plot  "<(sed -n '14209,14799p' unweighted.gr)" using 7:4  with lines title "OW-OH scattcg" lw 3, \
      'Mix-unweighted-OW-OH.rdf' using 1:2 with lines title "OW-OH Dissolve" lw 3

set output 'partialgr_HW-OH.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HW-OH}(r)"
#plot the graphic
plot  "<(sed -n '14801,15391p' unweighted.gr)" using 7:4  with lines title "HW-OH scattcg" lw 3, \
      'Mix-unweighted-HW-OH.rdf' using 1:2 with lines title "HW-OH Dissolve" lw 3

# The file we'll write to
set output 'partialgr_CT-OH.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{CT-OH}(r)"
#plot the graphic
plot  "<(sed -n '15393,15983p' unweighted.gr)" using 7:4  with lines title "CT-OH scattcg" lw 3, \
      'Mix-unweighted-CT-OH.rdf' using 1:2 with lines title "CT-OH Dissolve" lw 3

set output 'partialgr_HC-OH.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HC-OH}(r)"
#plot the graphic
plot  "<(sed -n '15985,16575p' unweighted.gr)" using 7:4  with lines title "HC-OH scattcg" lw 3, \
      'Mix-unweighted-HC-OH.rdf' using 1:2 with lines title "HC-OH Dissolve" lw 3

set output 'partialgr_OH-OH.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{OH-OH}(r)"
#plot the graphic
plot  "<(sed -n '16577,17167p' unweighted.gr)" using 7:4  with lines title "OH-OH scattcg" lw 3, \
      'Mix-unweighted-OH-OH.rdf' using 1:2 with lines title "OH-OH Dissolve" lw 3

# The file we'll write to
set output 'partialgr_OW-HO.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{OW-HO}(r)"
#plot the graphic
plot  "<(sed -n '17761,18351p' unweighted.gr)" using 7:4  with lines title "OW-HO scattcg" lw 3, \
      'Mix-unweighted-OW-HO.rdf' using 1:2 with lines title "OW-HO Dissolve" lw 3

set output 'partialgr_HW-HO.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HW-HO}(r)"
#plot the graphic
plot  "<(sed -n '18353,18943p' unweighted.gr)" using 7:4  with lines title "HW-HO scattcg" lw 3, \
      'Mix-unweighted-HW-HO.rdf' using 1:2 with lines title "HW-HO Dissolve" lw 3

set output 'partialgr_CT-HO.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{CT-HO}(r)"
#plot the graphic
plot  "<(sed -n '18945,19535p' unweighted.gr)" using 7:4  with lines title "CT-HO scattcg" lw 3, \
      'Mix-unweighted-CT-HO.rdf' using 1:2 with lines title "CT-HO Dissolve" lw 3

# The file we'll write to
set output 'partialgr_HC-HO.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel " g_{HC-HO}(r)"
#plot the graphic
plot  "<(sed -n '19537,20127p' unweighted.gr)" using 7:4  with lines title "HC-HO scattcg" lw 3, \
      'Mix-unweighted-HC-HO.rdf' using 1:2 with lines title "HC-HO Dissolve" lw 3

set output 'partialgr_OH-HO.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{OH-HO}(r)"
#plot the graphic
plot  "<(sed -n '20129,20719p' unweighted.gr)" using 7:4  with lines title "OH-HO scattcg" lw 3, \
      'Mix-unweighted-OH-HO.rdf' using 1:2 with lines title "OH-HO Dissolve" lw 3

set output 'partialgr_HO-HO.png'
# The graphic title
set title 'partial g(r)'
set xlabel "r (\305)"
set ylabel "g_{HO-HO}(r)"
#plot the graphic
plot  "<(sed -n '20721,21311p' unweighted.gr)" using 7:4  with lines title "HO-HO scattcg" lw 3, \
      'Mix-unweighted-HO-HO.rdf' using 1:2 with lines title "HO-HO Dissolve" lw 3



# The file we'll write to
set output 'partialsq_OW-OW.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-OW}(Q)"
#plot the graphic
plot  "<(sed -n '1,598p' unweighted.sq)" using 7:4  with lines title "OW-OW scattcg" lw 3, \
      'SQs-unweighted-OW-OW.sq' using 1:2 with lines title "OW-OW Dissolve" lw 3

set output 'partialsq_OW-HW.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OW-HW}(Q)"
#plot the graphic
plot  "<(sed -n '3589,4186p' unweighted.sq)" using 7:4  with lines title "OW-HW scattcg" lw 3, \
      'SQs-unweighted-OW-HW.sq' using 1:2 with lines title "OW-HW Dissolve" lw 3

set output 'partialsq_HW-HW.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HW}(Q)"
#plot the graphic
plot  "<(sed -n '4187,4784p' unweighted.sq)" using 7:4  with lines title "HW-HW scattcg" lw 3, \
      'SQs-unweighted-HW-HW.sq' using 1:2 with lines title "HW-HW Dissolve" lw 3

# The file we'll write to
set output 'partialsq_OW-CT.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q(1/\305)"
set ylabel " S_{OW-CT}(Q)"
#plot the graphic
plot  "<(sed -n '7177,7774p' unweighted.sq)" using 7:4  with lines title "OW-CT scattcg" lw 3, \
      'SQs-unweighted-OW-CT.sq' using 1:2 with lines title "OW-CT Dissolve" lw 3

set output 'partialsq_HW-CT.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-CT}(Q)"
#plot the graphic
plot  "<(sed -n '7775,8372p' unweighted.sq)" using 7:4  with lines title "HW-CT scattcg" lw 3, \
      'SQs-unweighted-HW-CT.sq' using 1:2 with lines title "HW-CT Dissolve" lw 3

set output 'partialsq_CT-CT.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{CT-CT}(Q)"
#plot the graphic
plot  "<(sed -n '8373,8970p' unweighted.sq)" using 7:4  with lines title "CT-CT scattcg" lw 3, \
      'SQs-unweighted-CT-CT.sq' using 1:2 with lines title "CT-CT Dissolve" lw 3

# The file we'll write to
set output 'partialsq_OW-HC.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-HC}(Q)"
#plot the graphic
plot  "<(sed -n '10765,11362p' unweighted.sq)" using 7:4  with lines title "OW-HC scattcg" lw 3, \
      'SQs-unweighted-OW-HC.sq' using 1:2 with lines title "OW-HC Dissolve" lw 3

set output 'partialsq_HW-HC.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HC}(Q)"
#plot the graphic
plot  "<(sed -n '11363,11960p' unweighted.sq)" using 7:4  with lines title "HW-HC scattcg" lw 3, \
      'SQs-unweighted-HW-HC.sq' using 1:2 with lines title "HW-HC Dissolve" lw 3

set output 'partialsq_CT-HC.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (\305)"
set ylabel "S_{CT-HC}(Q)"
#plot the graphic
plot  "<(sed -n '11961,12558p' unweighted.sq)" using 7:4  with lines title "CT-HC scattcg" lw 3, \
      'SQs-unweighted-CT-HC.sq' using 1:2 with lines title "CT-HC Dissolve" lw 3

# The file we'll write to
set output 'partialsq_HC-HC.sq'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{HC-HC}(Q)"
#plot the graphic
plot  "<(sed -n '12559,13156p' unweighted.sq)" using 7:4  with lines title "HC-HC scattcg" lw 3, \
      'SQs-unweighted-HC-HC.sq' using 1:2 with lines title "HC-HC Dissolve" lw 3

set output 'partialsq_OW-OH.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OW-OH}(Q)"
#plot the graphic
plot  "<(sed -n '14353,14950p' unweighted.sq)" using 7:4  with lines title "OW-OH scattcg" lw 3, \
      'SQs-unweighted-OW-OH.sq' using 1:2 with lines title "OW-OH Dissolve" lw 3

set output 'partialsq_HW-OH.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-OH}(Q)"
#plot the graphic
plot  "<(sed -n '14951,15548p' unweighted.sq)" using 7:4  with lines title "HW-OH scattcg" lw 3, \
      'SQs-unweighted-HW-OH.sq' using 1:2 with lines title "HW-OH Dissolve" lw 3

# The file we'll write to
set output 'partialsq_CT-OH.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (q/\305)"
set ylabel " S_{CT-OH}(Q)"
#plot the graphic
plot  "<(sed -n '15549,16146p' unweighted.sq)" using 7:4  with lines title "CT-OH scattcg" lw 3, \
      'SQs-unweighted-CT-OH.sq' using 1:2 with lines title "CT-OH Dissolve" lw 3

set output 'partialsq_HC-OH.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HC-OH}(Q)"
#plot the graphic
plot  "<(sed -n '16147,16744p' unweighted.sq)" using 7:4  with lines title "HC-OH scattcg" lw 3, \
      'SQs-unweighted-HC-OH.sq' using 1:2 with lines title "HC-OH Dissolve" lw 3

set output 'partialsq_OH-OH.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OH-OH}(Q)"
#plot the graphic
plot  "<(sed -n '16745,17342p' unweighted.sq)" using 7:4  with lines title "OH-OH scattcg" lw 3, \
      'SQs-unweighted-OH-OH.sq' using 1:2 with lines title "OH-OH Dissolve" lw 3

# The file we'll write to
set output 'partialsq_OW-HO.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-HO}(Q)"
#plot the graphic
plot  "<(sed -n '17941,18538p' unweighted.sq)" using 7:4  with lines title "OW-HO scattcg" lw 3, \
      'SQs-unweighted-OW-HO.sq' using 1:2 with lines title "OW-HO Dissolve" lw 3

set output 'partialsq_HW-HO.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HO}(Q)"
#plot the graphic
plot  "<(sed -n '18539,19136p' unweighted.sq)" using 7:4  with lines title "HW-HO scattcg" lw 3, \
      'SQs-unweighted-HW-HO.sq' using 1:2 with lines title "HW-HO Dissolve" lw 3

set output 'partialsq_CT-HO.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{CT-HO}(Q)"
#plot the graphic
plot  "<(sed -n '19137,19734p' unweighted.sq)" using 7:4  with lines title "CT-HO scattcg" lw 3, \
      'SQs-unweighted-CT-HO.sq' using 1:2 with lines title "CT-HO Dissolve" lw 3

# The file we'll write to
set output 'partialsq_HC-HO.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{HC-HO}(Q)"
#plot the graphic
plot  "<(sed -n '19735,20332p' unweighted.sq)" using 7:4  with lines title "HC-HO scattcg" lw 3, \
      'SQs-unweighted-HC-HO.sq' using 1:2 with lines title "HC-HO Dissolve" lw 3

set output 'partialsq_OH-HO.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OH-HO}(Q)"
#plot the graphic
plot  "<(sed -n '20333,20930p' unweighted.sq)" using 7:4  with lines title "OH-HO scattcg" lw 3, \
      'SQs-unweighted-OH-HO.sq' using 1:2 with lines title "OH-HO Dissolve" lw 3

set output 'partialsq_HO-HO.png'
# The graphic title
set title 'partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HO-HO}(Q)"
#plot the graphic
plot  "<(sed -n '20931,21528p' unweighted.sq)" using 7:4  with lines title "HO-HO scattcg" lw 3, \
      'SQs-unweighted-HO-HO.sq' using 1:2 with lines title "HO-HO Dissolve" lw 3



# The file we'll write to
set output 'weighted-partialsq_OW-OW.png'
# The graphic title
set title 'weigthed-partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-OW}(Q)"
#plot the graphic
plot  "<(sed -n '1,598p' weighted.sq)" using 7:4  with lines title "OW-OW scattcg" lw 3, \
      'HHH-weighted-OW-OW.sq' using 1:2 with lines title "OW-OW Dissolve" lw 3

set output 'weighted-partialsq_OW-HW.png'
# The graphic title
set title 'weigthed-partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OW-HW}(Q)"
#plot the graphic
plot  "<(sed -n '3589,4186p' weighted.sq)" using 7:4  with lines title "OW-HW scattcg" lw 3, \
      'HHH-weighted-OW-HW.sq' using 1:2 with lines title "OW-HW Dissolve" lw 3

set output 'weighted-partialsq_HW-HW.png'
# The graphic title
set title 'weigthed partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HW}(Q)"
#plot the graphic
plot  "<(sed -n '4187,4784p' weighted.sq)" using 7:4  with lines title "HW-HW scattcg" lw 3, \
      'HHH-weighted-HW-HW.sq' using 1:2 with lines title "HW-HW Dissolve" lw 3

# The file we'll write to
set output 'weighted-partialsq_OW-CT.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q(1/\305)"
set ylabel " S_{OW-CT}(Q)"
#plot the graphic
plot  "<(sed -n '7177,7774p' weighted.sq)" using 7:4  with lines title "OW-CT scattcg" lw 3, \
      'HHH-weighted-OW-CT.sq' using 1:2 with lines title "OW-CT Dissolve" lw 3

set output 'weighted-partialsq_HW-CT.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-CT}(Q)"
#plot the graphic
plot  "<(sed -n '7775,8372p' weighted.sq)" using 7:4  with lines title "HW-CT scattcg" lw 3, \
      'HHH-weighted-HW-CT.sq' using 1:2 with lines title "HW-CT Dissolve" lw 3

set output 'weighted-partialsq_CT-CT.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{CT-CT}(Q)"
#plot the graphic
plot  "<(sed -n '8373,8970p' weighted.sq)" using 7:4  with lines title "CT-CT scattcg" lw 3, \
      'HHH-weighted-CT-CT.sq' using 1:2 with lines title "CT-CT Dissolve" lw 3

# The file we'll write to
set output 'weighted-partialsq_OW-HC.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-HC}(Q)"
#plot the graphic
plot  "<(sed -n '10765,11362p' weighted.sq)" using 7:4  with lines title "OW-HC scattcg" lw 3, \
      'HHH-weighted-OW-HC.sq' using 1:2 with lines title "OW-HC Dissolve" lw 3

set output 'weighted-partialsq_HW-HC.png'
# The graphic title
set title 'weighted-partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HC}(Q)"
#plot the graphic
plot  "<(sed -n '11363,11960p' weighted.sq)" using 7:4  with lines title "HW-HC scattcg" lw 3, \
      'HHH-weighted-HW-HC.sq' using 1:2 with lines title "HW-HC Dissolve" lw 3

set output 'weighted-partialsq_CT-HC.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (\305)"
set ylabel "S_{CT-HC}(Q)"
#plot the graphic
plot  "<(sed -n '11961,12558p' weighted.sq)" using 7:4  with lines title "CT-HC scattcg" lw 3, \
      'HHH-weighted-CT-HC.sq' using 1:2 with lines title "CT-HC Dissolve" lw 3

# The file we'll write to
set output 'weighted-partialsq_HC-HC.sq'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{HC-HC}(Q)"
#plot the graphic
plot  "<(sed -n '12559,13156p' weighted.sq)" using 7:4  with lines title "HC-HC scattcg" lw 3, \
      'HHH-weighted-HC-HC.sq' using 1:2 with lines title "HC-HC Dissolve" lw 3

set output 'weighted-partialsq_OW-OH.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OW-OH}(Q)"
#plot the graphic
plot  "<(sed -n '14353,14950p' weighted.sq)" using 7:4  with lines title "OW-OH scattcg" lw 3, \
      'HHH-weighted-OW-OH.sq' using 1:2 with lines title "OW-OH Dissolve" lw 3

set output 'weighted partialsq_HW-OH.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-OH}(Q)"
#plot the graphic
plot  "<(sed -n '14951,15548p' weighted.sq)" using 7:4  with lines title "HW-OH scattcg" lw 3, \
      'HHH-weighted-HW-OH.sq' using 1:2 with lines title "HW-OH Dissolve" lw 3

# The file we'll write to
set output 'weighted partialsq_CT-OH.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (q/\305)"
set ylabel " S_{CT-OH}(Q)"
#plot the graphic
plot  "<(sed -n '15549,16146p' weighted.sq)" using 7:4  with lines title "CT-OH scattcg" lw 3, \
      'HHH-weighted-CT-OH.sq' using 1:2 with lines title "CT-OH Dissolve" lw 3

set output 'weighted-partialsq_HC-OH.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HC-OH}(Q)"
#plot the graphic
plot  "<(sed -n '16147,16744p' weighted.sq)" using 7:4  with lines title "HC-OH scattcg" lw 3, \
      'HHH-weighted-HC-OH.sq' using 1:2 with lines title "HC-OH Dissolve" lw 3

set output 'weighted-partialsq_OH-OH.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OH-OH}(Q)"
#plot the graphic
plot  "<(sed -n '16745,17342p' weighted.sq)" using 7:4  with lines title "OH-OH scattcg" lw 3, \
      'HHH-weighted-OH-OH.sq' using 1:2 with lines title "OH-OH Dissolve" lw 3

# The file we'll write to
set output 'weighted-partialsq_OW-HO.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-HO}(Q)"
#plot the graphic
plot  "<(sed -n '17941,18538p' weighted.sq)" using 7:4  with lines title "OW-HO scattcg" lw 3, \
      'HHH-weighted-OW-HO.sq' using 1:2 with lines title "OW-HO Dissolve" lw 3

set output 'weighted-partialsq_HW-HO.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HO}(Q)"
#plot the graphic
plot  "<(sed -n '18539,19136p' weighted.sq)" using 7:4  with lines title "HW-HO scattcg" lw 3, \
      'HHH-weighted-HW-HO.sq' using 1:2 with lines title "HW-HO Dissolve" lw 3

set output 'weighted-partialsq_CT-HO.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{CT-HO}(Q)"
#plot the graphic
plot  "<(sed -n '19137,19734p' weighted.sq)" using 7:4  with lines title "CT-HO scattcg" lw 3, \
      'HHH-weighted-CT-HO.sq' using 1:2 with lines title "CT-HO Dissolve" lw 3

# The file we'll write to
set output 'weighted-partialsq_HC-HO.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel " S_{HC-HO}(Q)"
#plot the graphic
plot  "<(sed -n '19735,20332p' weighted.sq)" using 7:4  with lines title "HC-HO scattcg" lw 3, \
      'HHH-weighted-HC-HO.sq' using 1:2 with lines title "HC-HO Dissolve" lw 3

set output 'weighted-partialsq_OH-HO.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{OH-HO}(Q)"
#plot the graphic
plot  "<(sed -n '20333,20930p' weighted.sq)" using 7:4  with lines title "OH-HO scattcg" lw 3, \
      'HHH-weighted-OH-HO.sq' using 1:2 with lines title "OH-HO Dissolve" lw 3

set output 'weighted-partialsq_HO-HO.png'
# The graphic title
set title 'weighted partial S(Q)'
set xlabel "Q (1/\305)"
set ylabel "S_{HO-HO}(Q)"
#plot the graphic
plot  "<(sed -n '20931,21528p' weighted.sq)" using 7:4  with lines title "HO-HO scattcg" lw 3, \
      'HHH-weighted-HO-HO.sq' using 1:2 with lines title "HO-HO Dissolve" lw 3

set output 'weighted-totalsq.png'
# The graphic title
set title 'weighted total S(Q)'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'weighted_total.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'HHH-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'weighted-totalgr.png'
# The graphic title
set title 'weighted total G(r)'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'weighted_total.gr' using 1:2  with lines title "scattcg" lw 3, \
      'HHH-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3


# The file we'll write to
set output 'deut-weighted-partialsq1_OW-OW.png'
# The graphic title
set title 'weigthed-partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-OW}(Q)"
#plot the graphic
plot  "<(sed -n '1,598p' deut_weighted1.sq)" using 7:4  with lines title "OW-OW scattcg" lw 3, \
      'DHH-weighted-OW-OW.sq' using 1:2 with lines title "OW-OW Dissolve" lw 3

set output 'deut-weighted-partialsq1_OW-HW.png'
# The graphic title
set title 'weigthed-partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{OW-HW}(Q)"
#plot the graphic
plot  "<(sed -n '3589,4186p' deut_weighted1.sq)" using 7:4  with lines title "OW-HW scattcg" lw 3, \
      'DHH-weighted-OW-HW.sq' using 1:2 with lines title "OW-HW Dissolve" lw 3

set output 'deut-weighted-partialsq1_HW-HW.png'
# The graphic title
set title 'weigthed partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HW}(Q)"
#plot the graphic
plot  "<(sed -n '4187,4784p' deut_weighted1.sq)" using 7:4  with lines title "HW-HW scattcg" lw 3, \
      'DHH-weighted-HW-HW.sq' using 1:2 with lines title "HW-HW Dissolve" lw 3

# The file we'll write to
set output 'deut-weighted-partialsq1_OW-CT.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q(1/\305)"
set ylabel " S_{OW-CT}(Q)"
#plot the graphic
plot  "<(sed -n '7177,7774p' deut_weighted1.sq)" using 7:4  with lines title "OW-CT scattcg" lw 3, \
      'DHH-weighted-OW-CT.sq' using 1:2 with lines title "OW-CT Dissolve" lw 3

set output 'deut-weighted-partialsq1_HW-CT.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-CT}(Q)"
#plot the graphic
plot  "<(sed -n '7775,8372p' deut_weighted1.sq)" using 7:4  with lines title "HW-CT scattcg" lw 3, \
      'DHH-weighted-HW-CT.sq' using 1:2 with lines title "HW-CT Dissolve" lw 3

set output 'deut-weighted-partialsq1_CT-CT.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{CT-CT}(Q)"
#plot the graphic
plot  "<(sed -n '8373,8970p' deut_weighted1.sq)" using 7:4  with lines title "CT-CT scattcg" lw 3, \
      'DHH-weighted-CT-CT.sq' using 1:2 with lines title "CT-CT Dissolve" lw 3

# The file we'll write to
set output 'deut-weighted-partialsq1_OW-HC.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-HC}(Q)"
#plot the graphic
plot  "<(sed -n '10765,11362p' deut_weighted1.sq)" using 7:4  with lines title "OW-HC scattcg" lw 3, \
      'DHH-weighted-OW-HC.sq' using 1:2 with lines title "OW-HC Dissolve" lw 3

set output 'deut-weighted-partialsq1_HW-HC.png'
# The graphic title
set title 'weighted-partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HC}(Q)"
#plot the graphic
plot  "<(sed -n '11363,11960p' deut_weighted1.sq)" using 7:4  with lines title "HW-HC scattcg" lw 3, \
      'DHH-weighted-HW-HC.sq' using 1:2 with lines title "HW-HC Dissolve" lw 3

set output 'deut-weighted-partialsq_CT-HC.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (\305)"
set ylabel "S_{CT-HC}(Q)"
#plot the graphic
plot  "<(sed -n '11961,12558p' deut_weighted1.sq)" using 7:4  with lines title "CT-HC scattcg" lw 3, \
      'DHH-weighted-CT-HC.sq' using 1:2 with lines title "CT-HC Dissolve" lw 3

# The file we'll write to
set output 'deut-weighted-partialsq1_HC-HC.sq'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel " S_{HC-HC}(Q)"
#plot the graphic
plot  "<(sed -n '12559,13156p' deut_weighted1.sq)" using 7:4  with lines title "HC-HC scattcg" lw 3, \
      'DHH-weighted-HC-HC.sq' using 1:2 with lines title "HC-HC Dissolve" lw 3

set output 'deut-weighted-partialsq1_OW-OH.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{OW-OH}(Q)"
#plot the graphic
plot  "<(sed -n '14353,14950p' deut_weighted1.sq)" using 7:4  with lines title "OW-OH scattcg" lw 3, \
      'DHH-weighted-OW-OH.sq' using 1:2 with lines title "OW-OH Dissolve" lw 3

set output 'deut-weighted-partialsq1_HW-OH.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-OH}(Q)"
#plot the graphic
plot  "<(sed -n '14951,15548p' deut_weighted1.sq)" using 7:4  with lines title "HW-OH scattcg" lw 3, \
      'DHH-weighted-HW-OH.sq' using 1:2 with lines title "HW-OH Dissolve" lw 3

# The file we'll write to
set output 'deut-weighted partialsq_CT-OH.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (q/\305)"
set ylabel " S_{CT-OH}(Q)"
#plot the graphic
plot  "<(sed -n '15549,16146p' deut_weighted1.sq)" using 7:4  with lines title "CT-OH scattcg" lw 3, \
      'DHH-weighted-CT-OH.sq' using 1:2 with lines title "CT-OH Dissolve" lw 3

set output 'deut-weighted-partialsq_HC-OH.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HC-OH}(Q)"
#plot the graphic
plot  "<(sed -n '16147,16744p' deut_weighted1.sq)" using 7:4  with lines title "HC-OH scattcg" lw 3, \
      'DHH-weighted-HC-OH.sq' using 1:2 with lines title "HC-OH Dissolve" lw 3

set output 'deut-weighted-partialsq1_OH-OH.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{OH-OH}(Q)"
#plot the graphic
plot  "<(sed -n '16745,17342p' deut_weighted1.sq)" using 7:4  with lines title "OH-OH scattcg" lw 3, \
      'DHH-weighted-OH-OH.sq' using 1:2 with lines title "OH-OH Dissolve" lw 3

# The file we'll write to
set output 'deut-weighted-partialsq1_OW-HO.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel " S_{OW-HO}(Q)"
#plot the graphic
plot  "<(sed -n '17941,18538p' deut_weighted1.sq)" using 7:4  with lines title "OW-HO scattcg" lw 3, \
      'DHH-weighted-OW-HO.sq' using 1:2 with lines title "OW-HO Dissolve" lw 3

set output 'deut-weighted-partialsq1_HW-HO.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HW-HO}(Q)"
#plot the graphic
plot  "<(sed -n '18539,19136p' deut_weighted1.sq)" using 7:4  with lines title "HW-HO scattcg" lw 3, \
      'DHH-weighted-HW-HO.sq' using 1:2 with lines title "HW-HO Dissolve" lw 3

set output 'deut-weighted-partialsq1_CT-HO.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{CT-HO}(Q)"
#plot the graphic
plot  "<(sed -n '19137,19734p' deut_weighted1.sq)" using 7:4  with lines title "CT-HO scattcg" lw 3, \
      'DHH-weighted-CT-HO.sq' using 1:2 with lines title "CT-HO Dissolve" lw 3

# The file we'll write to
set output 'deut-weighted-partialsq1_HC-HO.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel " S_{HC-HO}(Q)"
#plot the graphic
plot  "<(sed -n '19735,20332p' deut_weighted1.sq)" using 7:4  with lines title "HC-HO scattcg" lw 3, \
      'DHH-weighted-HC-HO.sq' using 1:2 with lines title "HC-HO Dissolve" lw 3

set output 'deut-weighted-partialsq1_OH-HO.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{OH-HO}(Q)"
#plot the graphic
plot  "<(sed -n '20333,20930p' deut_weighted1.sq)" using 7:4  with lines title "OH-HO scattcg" lw 3, \
      'DHH-weighted-OH-HO.sq' using 1:2 with lines title "OH-HO Dissolve" lw 3

set output 'deut-weighted-partialsq1_HO-HO.png'
# The graphic title
set title 'weighted partial S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "S_{HO-HO}(Q)"
#plot the graphic
plot  "<(sed -n '20931,21528p' deut_weighted1.sq)" using 7:4  with lines title "HO-HO scattcg" lw 3, \
      'DHH-weighted-HO-HO.sq' using 1:2 with lines title "HO-HO Dissolve" lw 3

set output 'deut-weighted-totalsq1.png'
# The graphic title
set title 'weighted total S(Q)-DHH'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total1.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'DHH-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr1.png'
# The graphic title
set title 'weighted total G(r)-DHH'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total1.gr' using 1:2  with lines title "scattcg" lw 3, \
      'DHH-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3

set output 'deut-weighted-totalsq2.png'
# The graphic title
set title 'weighted total S(Q)-HDH'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total2.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'HDH-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr2.png'
# The graphic title
set title 'weighted total G(r)-HDH'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total2.gr' using 1:2  with lines title "scattcg" lw 3, \
      'HDH-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3

set output 'deut-weighted-totalsq3.png'
# The graphic title
set title 'weighted total S(Q)-HHD'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total3.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'HHD-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr3.png'
# The graphic title
set title 'weighted total G(r)-HHD'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total3.gr' using 1:2  with lines title "scattcg" lw 3, \
      'HHD-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3

set output 'deut-weighted-totalsq4.png'
# The graphic title
set title 'weighted total S(Q)-DDH'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total4.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'DDH-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr4.png'
# The graphic title
set title 'weighted total G(r)-DDH'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total4.gr' using 1:2  with lines title "scattcg" lw 3, \
      'DDH-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3

set output 'deut-weighted-totalsq5.png'
# The graphic title
set title 'weighted total S(Q)-HDD'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total5.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'HDD-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr5.png'
# The graphic title
set title 'weighted total G(r)-HDD'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total5.gr' using 1:2  with lines title "scattcg" lw 3, \
      'HDD-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3

set output 'deut-weighted-totalsq6.png'
# The graphic title
set title 'weighted total S(Q)-DDD'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total6.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'DDD-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr6.png'
# The graphic title
set title 'weighted total G(r)-DDD'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total6.gr' using 1:2  with lines title "scattcg" lw 3, \
      'DDD-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3

set output 'deut-weighted-totalsq7.png'
# The graphic title
set title 'weighted total S(Q)-H5H'
set xlabel "Q (1/\305)"
set ylabel "F_(Q)"
#plot the graphic
plot 'deut_weighted_total7.sq'  using 1:2  with lines title " scattcg" lw 3, \
      'H5H-weighted-total.sq' using 1:2 with lines title " Dissolve" lw 3

set output 'deut-weighted-totalgr7.png'
# The graphic title
set title 'weighted total G(r)-H5H'
set xlabel "r (\305)"
set ylabel "G_(r)"
#plot the graphic
plot  'deut_weighted_total7.gr' using 1:2  with lines title "scattcg" lw 3, \
      'H5H-weighted-total.gr.broad' using 1:2 with lines title "Dissolve" lw 3


