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
selectcase(output)
case(1)
   call flst
case(0)
   !call postDat
   allocate(times(0:ndata))
end select
tecplot=.false.; ispost=.true.
nread=0
! RPT-READ X,Y,Z COORDINATES
selectcase(output)
case(0)
   do nn = 1, 3
      if (tecplot) then
      nread=nread+1; call tpostread(nread,lsta)
      else
      nread=nread+1; call postread(nread)
      end if
      ss(:,nn)=varr(:)
   end do
case(1)
   call p3dread(gsflag=1,nout=0)
end select

call getMetrics
ngridv=1
if (myid==7) then
   write(*,*) sum(area)
end if

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
if (favg==1) then
selectcase(output)
case(0)
!call average
case(1)
call p3daverage
end select
end if

!!===== WRITE AVERAGE VALUES 
if (fwavg==1) then
selectcase(output)
case(0)
do nn=1,5
   varr(:)=qa(:,nn)
   nwrec=nwrec+1
   call postwrite(nwrec)
end do
case(1)
   call p3dwaverage
end select
end if

!===COMPUTE FORCE COEFFICIENT
if (fcoef==1) then
do n = 0, ndata+1
call clpost(ele=1,nvar=n)
if (myid==0) then
     ra0=aoa*pi/180;ra1=cos(ra0);ra2=sin(ra0)
   write(*,*) cl(1,2)*ra1-cl(1,1)*ra2,cl(1,2)*ra2+cl(1,1)*ra1,n
end if
end do
end if

!===find location
if (floc==1) then
call findll(0.5_nr,-0.05_nr,0.0_nr,l,m)
if (myid==m) then
   write(*,"(f6.2,1x,f6.2,1x,f6.2,1x,i3,i10)") xyz(l,1),xyz(l,2),xyz(l,3),m,l
    open(7,file='data/signal.dat')
    write(7,"('t uprime')") 
   do n = 0, ndata
      res=0
      selectcase(output)
      case(0)
         do i = 2, 2
         nread=nrec+i+(n*totVar)
         call postread(nread)
         res=(varr(l)-qa(l,i))
         end do
      case(1)
         call p3dread(0,n)
         res=(qo(l,2)-qa(l,2))
      end select
      write(7,"(f10.5,f10.5)") times(n),res
   end do
   close(7)
end if
end if

!==COMPUTE WALL DISTANCES
if (fwplus==1) then
call getwplus(nvar=ndata+1)
if (myid==7) then
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
end if

!==COMUPTE Q-CRITERION
if (fqcrit==1) then
do n = ndata+1, ndata+1
call qcriterion(n)
selectcase(output)
case(0)
nwrec=nwrec+1; call postwrite(nwrec) ! ADD A LINE IN POST SUBROUTINE
case(1)
  cinput='Q'; qo(:,1)=varr(:); call wffile(cinput,n,1)
end select
end do
end if

!==WRITE WSS
if (fwss==1) then
if (output==1) then
do n = ndata+1, ndata+1
   call gettw(n)
   do i = 1, 3
      qo(:,i)=0
      if (wflag) then
         do m = 0, lcwall; l=lwall(m)
            qo(l,i)=tw(m,i)
         end do
      end if
   end do
  cinput='tw'; call wffile(cinput,n,3)
end do
end if
end if

!==COMPUTE Cf 
if (fcf==1) then
if (myid==7) then
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
end if

!==COMPUTE Cp
if (fcp==1) then
   do n = ndata+1, ndata+1
      call fillqo(n)
      ra0=two/(amachoo**2)
      qo(:,1)=(p(:)-poo)*ra0
      cinput='Cp'; call wffile(cinput,n,1)
   end do
end if

!==COMPUTE VORTICITY
if (fcurl==1) then
   do n = ndata+1, ndata+1
      call getCurl(n)
      qo(:,1:3)=ss(:,1:3)
      cinput='omega'; call wffile(cinput,n,3)
   end do
end if

!===== WRITE TECPLOT FILE
CALL MPI_BARRIER(icom,ierr)
if (output==0) then
close(9)
call post(average=.false.)
end if

!===== END OF JOB
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
