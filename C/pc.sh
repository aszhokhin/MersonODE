rm d4eq
rm D4EQ.*
rm d4eqg.ps
gcc -od4eq d4eq.c
./d4eq
gnuplot d4eq.gp
#evince d4eqg.ps