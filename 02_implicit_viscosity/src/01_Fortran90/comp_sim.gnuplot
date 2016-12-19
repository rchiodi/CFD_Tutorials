set term pngcairo enhanced dashed dl 0.5
set key bot right font " ,18pt"
set xlabel "U/U_{lid}" font " ,18pt"
set ylabel "y/L" font " ,18pt"
set xrange[-0.4:1.0]
set yrange[-0.5:0.5]
set xtics font " ,14pt"
set ytics font " ,14pt"
set title "Midline U Velocity for 50x50 Sim., Re = 100" font " ,18pt"

set linestyle 1 linetype 1 lw 5 dashtype 2 pt 3 ps 1.5 linecolor rgb "blue"
set linestyle 3 linetype 1 lw 3 dashtype 1 pt 3 ps 1.5 linecolor rgb "black"
set linestyle 2 linetype 1 lw 3 dashtype 1 pt 5 ps 1.5 linecolor rgb "red"

set output "comp.png"

plot "U_midplane_implicit" using 2:1 linestyle 3 with lines title "Sim. Implicit Visc.", \
"U_midplane_explicit" using 2:1 linestyle 1 with lines title "Sim. Explicit Visc.", \
"ghia1982" using 2:($1-0.5) linestyle 2 with points title "Ghia et al. 1982" 

