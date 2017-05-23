rm d4eq
rm D4EQ.*
gfortran -od4eq d4eq.f
./d4eq
gnuplot d4eq.gp
#evince d4eqg.ps
