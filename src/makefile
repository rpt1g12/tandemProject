#SHELL	= /bin/csh
PGM	= ../run.exe
PGMS	= ../srun.exe
PGMPS	= ../spost.exe
PGMI	= ../inter.exe
#FC	= mpif90
FC	= ftn
#FLAGS	= -g -fbacktrace -fbounds-check -O4 -o
#FLAGS	= -g -traceback -O5 -o
FLAGS	= -O5 -o
OBJN	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 Grid.f90 SingleTurbulence3D.f90 rpt.f90 SubSets.f90 Main3D.f90
OBJPS	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 Grid.f90 SingleTurbulence3D.f90 rpt.f90 rptPost.f90 SubSets.f90 MainPost.f90
OBJI	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 Grid.f90 SingleTurbulence3D.f90 rpt.f90 rptInter.f90 Interpolate3D.f90
all:
	make new
post:
	$(FC) $(FLAGS) $(PGMPS) $(OBJPS)
	rm *.mod
new:
	$(FC) $(FLAGS) $(PGMS) $(OBJN)
	rm *.mod
inter:
	$(FC) $(FLAGS) $(PGMI) $(OBJI)
	rm *.mod
