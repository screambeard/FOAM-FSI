#set term x11 0
set yr [:]
set key bottom right
set xlabel "Simulationtime [s]"
set ylabel "Cl [-]"
set title "Cl vs time"
set grid

plot 	"./forceCoeffs/0/forceCoeffs.dat" using ($1):($3) with lines title "Cl"

set term x11 1
set yr [:]
set key bottom right
set xlabel "Simulationtime [s]"
set ylabel "Cd [-]"
set title "Cd vs time"
set grid

plot    "./forceCoeffs/0/forceCoeffs.dat" using ($1):($2) with lines title "Cd"

pause 1
reread
