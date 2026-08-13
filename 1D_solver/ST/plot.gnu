#! /usr/bin/gnuplot/
set encoding utf8
set term pngcairo font "Times-New-Roman,18" size 1000, 800 lw 1
set tics font "Times-New-Roman,18"
set key off
set xtics nomirror
set ytics nomirror

set xlabel "{/:Italic x / L_x}"
set xrange [0:1]

set output "result.png"
  
set multiplot layout 2,2

# --- density ---
set ylabel "{/:Italic ρ / ρ_L}"
set yrange [0.1:1.1]
plot "Q.dat" using 1:2 w l lw 2 lc rgb "blue"

# --- velocity ---
set ylabel "{/:Italic u / a_L}"
set yrange [0:1]
plot "Q.dat" using 1:3 w l lw 2 lc rgb "blue"

# --- pressure ---
set ylabel "{/:Italic p / p_L}"
set yrange [0:1.1]
plot "Q.dat" using 1:4 w l lw 2 lc rgb "blue"

unset multiplot

