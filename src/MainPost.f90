!*****
!***** 3D PARALLEL SOLVER
!*****

 program mainpost

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 use rpt
 use rptpost
 implicit none

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

tecplot=.false.
call setup
nread=0
do nn = 1, 3
   nread=nread+1; call postread(nread)
   ss(:,nn)=varr(:)
end do

call getMetrics
ngridv=1

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
 call average

!===== WRITE AVERAGE VALUES (MAKE SURE THIS IS THE LAST WRITTEN!)
 do nn=1,5
    varr(:)=qa(:,nn)
    nwrec=nwrec+1; call postwrite(nwrec)
 end do

!===COMPUTE FORCE COEFFICIENT
call clpost(ele=1,dir=2,nvar=1,post=.true.)

!==COMPUTE WALL DISTANCES
call getwplus(nvar=1)
if (myid==11) then
    open(7,file='out/wplus.dat')
    write(7,"('x y z x+ y+ z+')") 
 do nn = 0, lcwall;l=lwall(nn)
 write(7,"(f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)")&
       xyz(l,1),xyz(l,2),xyz(l,3),wplus(nn,1),wplus(nn,2),wplus(nn,3)
 end do
 close(7)
 if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
 wvarr=wplus(:,2)
 call wavg(dir=2,wall=.true.)
end if

!==COMPUTE SPACE AVERAGE
!if (myid==11) then
!varr(:)=qa(:,5)
!call spcavg(plane=2,dir=3,pos=0)
!end if

!==COMUPTE Q-CRITERION
call qcriterion(1)
nwrec=nwrec+1; call postwrite(nwrec) ! ADD A LINE IN POST SUBROUTINE

!==WRITE DATA SO IT CAN BE READ BY POST SUBROUTINE
if (tecplot) then
   call tdatawrite(lsta)
else
   call datawrite
end if
close(8)

!!===== WRITE TECPLOT FILE
 call post(average=.true.)

!===== END OF JOB
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
