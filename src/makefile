SHELL	= /bin/csh
PGM	= ../run.exe
PGMS	= ../srun.exe
PGMP	= ../post.exe
PGMI	= ../inter.exe
FC	= mpif90
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 Main3D.f90
OBJS	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridSingle.f90 SingleTurbulence3D.f90 rpt.f90 Main3D.f90
OBJP	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 rptPost.f90 MainPost.f90
OBJI	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridSingle.f90 SingleTurbulence3D.f90 rpt.f90 rptInter.f90 Interpolate3D.f90
all:
	$(FC) -O3 -traceback -i_dynamic -o $(PGMP) $(OBJP)
	rm *.mod
run:
	$(FC) -O3 -traceback -i_dynamic -o $(PGM) $(OBJ)
	rm *.mod
post:
	$(FC) -O3 -traceback -i_dynamic -o $(PGMP) $(OBJP)
	rm *.mod
single:
	$(FC) -O3 -traceback -i_dynamic -o $(PGMS) $(OBJS)
	rm *.mod
inter:
	$(FC) -O3 -traceback -i_dynamic -o $(PGMI) $(OBJI)
	rm *.mod
clean: 
	rm batch.*
	rm display
cleanmisc: 
	rm misc/*
cleanout: 
	rm output*
	
