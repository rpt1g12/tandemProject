#SHELL	= /bin/csh
PGM	= ../run.exe
PGMS	= ../srun.exe
PGMG	= ../grun.exe
PGMP	= ../post.exe
PGMPS	= ../spost.exe
PGMI	= ../inter.exe
FC	= mpif90
#FC	= ftn
FLAGS	= -g -fbacktrace -fbounds-check -O4 -o
OBJ	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 Main3D.f90
OBJS	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridSingle.f90 SingleTurbulence3D.f90 rpt.f90 Main3D.f90
OBJN	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 Grid.f90 SingleTurbulence3D.f90 rpt.f90 Main3D.f90
OBJG	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridGeneral.f90 GeneralTurbulence3D.f90 rpt.f90 Main3D.f90
OBJP	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridTandem.f90 TandemTurbulence3D.f90 rpt.f90 rptPost.f90 MainPost.f90
OBJPS	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridSingle.f90 SingleTurbulence3D.f90 rpt.f90 rptPost.f90 MainPost.f90
OBJI	= MainVar3D.f90 Subroutineso.f90 Subroutines3D.f90 GridSingle.f90 SingleTurbulence3D.f90 rpt.f90 rptInter.f90 Interpolate3D.f90
all:
	make new
run:
	$(FC) $(FLAGS) $(PGM) $(OBJ)
	rm *.mod
post:
	$(FC) $(FLAGS) $(PGMP) $(OBJP)
	rm *.mod
spost:
	$(FC) $(FLAGS) $(PGMPS) $(OBJPS)
	rm *.mod
single:
	$(FC) $(FLAGS) $(PGMS) $(OBJS)
	rm *.mod
new:
	$(FC) $(FLAGS) $(PGMS) $(OBJN)
	rm *.mod
inter:
	$(FC) $(FLAGS) $(PGMI) $(OBJI)
	rm *.mod
test:
	$(FC) $(FLAGS) $(PGMG) $(OBJG)
	rm *.mod
general:
	$(FC) $(FLAGS) $(PGMG) $(OBJG)
	rm *.mod
