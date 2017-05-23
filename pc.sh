rm d3eq
rm D3EQ.*
gfortran -od3eq d3eq.f
./d3eq
gnuplot d3eq.gp
#evince d3eqg.ps
