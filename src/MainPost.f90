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

call setup
call postDat
tecplot=.false.; ispost=.true.
nread=0
! RPT-READ X,Y,Z COORDINATES
do nn = 1, 3
   if (tecplot) then
   nread=nread+1; call tpostread(nread,lsta)
   else
   nread=nread+1; call postread(nread)
   end if
   ss(:,nn)=varr(:)
end do

call getMetrics
ngridv=1
if (myid==11) then
   write(*,*) sum(area)
end if

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
call average

!===== WRITE AVERAGE VALUES (MAKE SURE THIS IS THE LAST WRITTEN!)
do nn=1,5
   varr(:)=qa(:,nn)
   nwrec=nwrec+1; call postwrite(nwrec)
end do

!===COMPUTE FORCE COEFFICIENT
ra0=0;ra1=0
do n = 0, ndata
call clpost(ele=1,dir=2,nvar=n)
call clpost(ele=1,dir=1,nvar=n)
if (myid==0) then
   write(*,*) cl(1,2),cl(1,1),n
end if
ra0=ra0+cl(1,2);ra1=ra1+cl(1,1)
end do
call clpost(ele=1,dir=2,nvar=ndata+1)
call clpost(ele=1,dir=1,nvar=ndata+1)
if (myid==0) then
   write(*,*) cl(1,2),cl(1,1),ra0/real(ndata+1),ra1/real(ndata+1)
end if

!===find location
!call findll(0.5_nr,-0.05_nr,0.0_nr,l,m)
!if (myid==m) then
!   write(*,"(f6.2,1x,f6.2,1x,f6.2,1x,i3,i10)") xyz(l,1),xyz(l,2),xyz(l,3),m,l
!    open(7,file='data/signal.dat')
!    write(7,"('t uprime')") 
!   do n = 0, ndata
!      res=0
!      do i = 2, 2
!      nread=nrec+i+(n*totVar)
!      call postread(nread)
!      res=(varr(l)-qa(l,i))
!      end do
!      write(7,"(f10.5,f10.5)") times(n),res
!   end do
!   close(7)
!end if
!
!!==COMPUTE WALL DISTANCES
!call getwplus(nvar=26)
!if (myid==11) then
!    open(7,file='out/wplus.dat')
!    write(7,"('x y z x+ y+ z+')") 
! do nn = 0, lcwall;l=lwall(nn)
! write(7,"(f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)")&
!       xyz(l,1),xyz(l,2),xyz(l,3),wplus(nn,1),wplus(nn,2),wplus(nn,3)
! end do
! close(7)
! if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
! wvarr=wplus(:,2)
! call wavg(dir=2,wall=.true.)
!end if

!==COMPUTE SPACE AVERAGE
!if (myid==11) then
!varr(:)=qa(:,5)
!call spcavg(plane=2,dir=3,pos=0)
!end if

!==COMUPTE Q-CRITERION
!do n = 0, ndata
!call qcriterion(n)
!nwrec=nwrec+1; call postwrite(nwrec) ! ADD A LINE IN POST SUBROUTINE
!end do

!==COMPUTE Cf 
if (myid==11) then
call gettw(ndata+1)
   open(7,file='data/Cf.dat')
   write(7,"('x cf')") 
   ra0=two/(amachoo**2)
   do n = 0, lcwall; l=lwall(n)
      ra1=DOT_PRODUCT(tw(n,:),wtan(n,:))
      write(7,"(f10.5,f10.5)")  xyz(l,1),-ra1*ra0
   end do
   close(7)
end if

!==COMPUTE Cp
!if (myid==11) then
!call fillqo(ndata+1)
!varr(:)=qo(:,5)
!call getatWall
!   open(7,file='data/Cp.dat')
!   write(7,"('x cp')") 
!   ra0=two/(amachoo**2)
!   do n = 0, lcwall; l=lwall(n)
!      write(7,"(f10.5,f10.5)")  xyz(l,1),(wvarr(n)-poo)*ra0
!   end do
!   close(7)
!end if

!===== WRITE TECPLOT FILE
CALL MPI_BARRIER(icom,ierr)
close(9)
!call post(average=.false.)

!===== END OF JOB
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
