SHELL	= /bin/csh
PGM	= ../run.exe
PGMP	= ../post.exe
FC	= ftn
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 Main3D.f90
OBJP	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 MainPost.f90
$(PGM): $(OBJ)
	$(FC) -fast -o $(PGM) $(OBJ)
	rm *.mod
clean: 
	rm batch.*
	rm display
cleanmisc: 
	rm misc/*
cleanout: 
	rm output*
post:
	$(FC) -O3 -traceback -i_dynamic -o $(PGMP) $(OBJP)
	rm *.mod
	
