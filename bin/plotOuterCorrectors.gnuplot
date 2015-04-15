set xlabel "Simulationtime [s]"
set ylabel "# outerCorrectors"
set title "Plot of nOuterCorrectors over simulationtime"
set grid

plot    "./logs/nOuterCorrectors" using ($1):($2) with lines title "#outerCorr"	

pause 1
reread