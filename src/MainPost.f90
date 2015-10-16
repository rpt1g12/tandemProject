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
call flst
ispost=.true.
nread=0
! RPT-READ X,Y,Z COORDINATES
call p3dread(gsflag=1,nout=0)

call getMetrics
ngridv=1
if (myid==7) then
   write(*,*) sum(area)
end if

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
if (favg==1) then
   call p3daverage
end if

!!===== WRITE AVERAGE VALUES 
if (fwavg==1) then
   call p3dwaverage
end if

!===== COMPUTE RMS VALUES IF NOT AVAILABLE YET
if (frms==1) then
   call p3drms
end if

!!===== WRITE RMS VALUES 
if (fwrms==1) then
   call p3dwrms
end if

!===COMPUTE FORCE COEFFICIENT
if (fcoef==1) then
   if (myid==0) then
      write(*,"(3x,'n',8x,'time',9x,'Cl',9x,'Cd',5x)")  
      write(*,"('============================================')")
   end if
   do n = 0, ndata+favgu
      call clpost(ele=1,nvar=n)
      if (myid==0) then
         ra0=aoa*pi/180;ra1=cos(ra0);ra2=sin(ra0); 
         if(n==ndata+favgu) then
            ra3=(-1)
         else
            ra3=times(n)
         end if
         write(*,"(i8,f12.5,f12.7,f12.7)") &
         n,ra3,cl(1,2)*ra1-cl(1,1)*ra2,cl(1,2)*ra2+cl(1,1)*ra1
      end if
   end do
end if

!===find location
if (floc==1) then
   call findll(0.5_k8,-0.05_k8,0.0_k8,l,m)
   if (myid==m) then
      write(*,"(f6.2,1x,f6.2,1x,f6.2,1x,i3,i10)") xyz(l,1),xyz(l,2),xyz(l,3),m,l
      open(7,file='data/signal.dat')
      write(7,"('t uprime')") 
      do n = 0, ndata
         res=0
         call p3dread(0,n)
         res=(qo(l,2)-qa(l,2))
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
      call wavg('y+',dir=2,wall=.true.)
   end if
end if

!==COMUPTE Q-CRITERION
if (fqcrit==1) then
   do n = ndata+1, ndata+1
      call qcriterion(n)
      cinput='Q'; qo(:,1)=varr(:); call wffile(cinput,n,1)
   end do
end if

!==WRITE WSS
if (fwss==1) then
   do n = ndata+1, ndata+1
      call gettw(n)
      do i = 1, 3
         qo(:,i)=0
         if (myid==7) then
            do m = 0, lcwall; l=lwall(m)
               qo(l,i)=tw(m,i)
            end do
         end if
      end do
     cinput='tw'; call wffile(cinput,n,3)
   end do
end if

!==COMPUTE Cf 
if (fcf==1) then
   if (myid==7) then
      call gettw(ndata+1)
      open(7,file='data/allCf.dat')
      write(7,"('x cf')") 
      ra0=two/(amachoo**2)
      if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
      do n = 0, lcwall; l=lwall(n)
         ra1=DOT_PRODUCT(tw(n,:),wtan(n,:))
         wvarr(n)=ra1*ra0
         write(7,"(f10.5,f10.5)")  xyz(l,1),ra1*ra0
      end do
      close(7)
      call wavg('Cf',dir=2,wall=.true.)
   end if
end if

!==COMPUTE Cp
if (fcp==1) then
   call fillqo(ndata+1)
   ra0=two/(amachoo**2)
   qo(:,1)=(p(:)-poo)*ra0
   cinput='Cp'; call wffile(cinput,ndata+1,1)
   if (myid==7) then
      varr(:)=qo(:,1)
      call wavg('Cp',dir=2,wall=.false.)
   end if
end if

!===== END OF JOB
 CALL MPI_BARRIER(icom,ierr)
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
