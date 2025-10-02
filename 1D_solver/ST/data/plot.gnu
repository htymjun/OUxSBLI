#! /usr/bin/gnuplot/
set term pngcairo size 1280, 960
set output "output.png"

set format xy "%4.1f"
set key off

set xlabel "{/=12 {/Times-New-Roman:Italic x}}"

set xrange [0:1]

set tics font "Times-New-Roman,10"
set xtics nomirror
set ytics nomirror

datafile = "Q00001.d"

set multiplot layout 2, 3

rho0 = 1.293e0
p0   = rho0 * 287.03e0 * 300.e0
u0   = sqrt(p0 / rho0)

set ylabel "{/=12 {/Symbol:Italic r}}"
set size square
set mxtics 5
set mytics 5
set yrange [0:1.1]
plot datafile using 1:($2/rho0) with lines lw 1 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic u}}"
set size square
set mxtics 5
set mytics 5
set yrange [0:1]
plot datafile using 1:($3/u0) with lines lw 1 lc rgb "blue"

set ylabel "{/=12 {/Times-New-Roman:Italic P}}"
set size square
set mxtics 5
set mytics 5
set yrange [0:1.1]
plot datafile using 1:($4/p0) with lines lw 1 lc rgb "blue"

unset multiplot

