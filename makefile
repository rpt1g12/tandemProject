SHELL	= /bin/csh
PGM	= run.exe
FC	= mpif90
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandemAerofoil.f90 TandemTurbulence3D.f90 rpt.f90 Main3D.f90
$(PGM): $(OBJ)
	$(FC) -O3 -i_dynamic -o $(PGM) $(OBJ)
	rm *.mod
clean: 
	rm batch.*
	rm display
cleanmisc: 
	rm misc/*
cleanout: 
	rm output*
