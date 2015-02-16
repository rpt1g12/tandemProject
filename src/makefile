SHELL	= /bin/csh
PGM	= ../run.exe
PGMP	= ../post.exe
FC	= mpif90
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 Main3D.f90
OBJP	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 rptPost.f90 MainPost.f90
all:
	$(FC) -O3 -traceback -i_dynamic -o $(PGM) $(OBJ)
	rm *.mod
post:
	$(FC) -O3 -traceback -i_dynamic -o $(PGMP) $(OBJP)
	rm *.mod
clean: 
	rm batch.*
	rm display
cleanmisc: 
	rm misc/*
cleanout: 
	rm output*
	
