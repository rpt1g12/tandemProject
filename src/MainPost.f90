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
 character(20) :: cformat
 integer :: wmaster

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

call setup
call flst(fmblk)
   
!! RPT-READ X,Y,Z COORDINATES
!call rdP3dG(fmblk)
ispost=.true.
!nread=0
!call p3dread(gsflag=1,nout=0)
!
!call getMetrics
!CALL MPI_BARRIER(icom,ierr)
!!ngridv=1
!if (mb==7) then
!   if(wflag) then
!      ra1=sum(area)
!      CALL MPI_REDUCE(ra1,ra0,1,MPI_REAL8,MPI_SUM,0,bwcom,ierr)
!      call MPI_COMM_RANK(bwcom,l,ierr)
!      if(l==0) then
!         write(*,*) ra0,myid
!      end if
!   end if
!end if

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
if (favg==1) then
   call p3daverage
end if

!!===== WRITE AVERAGE VALUES 
if (fwavg==1) then
   qo(:,:)=qa(:,:)
   call wrP3dP(ndata+1,fmblk)
   !call p3dwaverage
end if

!===== COMPUTE RMS VALUES IF NOT AVAILABLE YET
if (frms==1) then
   call p3drms
end if

!!===== WRITE RMS VALUES 
if (fwrms==1) then
   qo(:,:)=qb(:,:)
   call wrP3dP(ndata+2,fmblk)
   !call p3dwrms
end if

!===COMPUTE FORCE COEFFICIENT
if (fcoef==1) then
   if (myid==0) then
      write(*,"(3x,'n',8x,'time',9x,'Cl',9x,'Cd',5x)")  
      write(*,"('============================================')")
   end if
   do n = 0, ndata+favgu
      call clpost(ele=1,nvar=n)
      if (myid==13) then
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
!if (allocated(fout)) deallocate(fout)
!if (.not.allocated(fout)) allocate(fout(0:lmx,1))
   if (mb==7) then
      do n = 0, ndata
         call qcriterion(n); fout(:,1)=varr(:)
         !call wrP3dF('Q',n,1,fmblk)
      end do
   end if
end if

!==COMUPTE VORTICITY + Q
if (fcurl==1) then
   if (allocated(fout)) deallocate(fout)
   if (.not.allocated(fout)) allocate(fout(0:lmx,1))
   if (mb==7) then
   do n = 0, ndata
      call qcriterion(n);
      call getCurl(n);
      call wrP3dP(n,fmblk,'Q+W')
      !fout(:,1:3)=qo(:,2:4)
      !call wrP3dF('Omega',n,3,fmblk)
   end do
   end if
end if

!==WRITE WSS
if (fwss==1) then
   do n = ndata+1, ndata+2
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

!==COMPUTE Cf 
if (fcf==1) then
   if (myid==7) then
      call gettw(ndata+1)
      open(7,file='data/allCfAVG.dat')
      write(7,"('x cf')") 
      ra0=two/(amachoo**2)
      if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
      do n = 0, lcwall; l=lwall(n)
         ra1=DOT_PRODUCT(tw(n,:),wtan(n,:))
         wvarr(n)=ra1*ra0
         write(7,"(f10.5,f10.5)")  xyz(l,1),ra1*ra0
      end do
      close(7)
      call wavg('CfAVG',dir=2,wall=.true.)
   end if
end if

!==COMPUTE RMS Cp
if (fcp==1) then
   call getCp(ndata+1)
   qa(:,1)=qo(:,1)
   qb(:,1)=0
    if (myid==0) then
       write(*,"('Total amout of data: ',i3)") ndata
    end if
    ! CONSTRUCT THE COEFFICIENTS ARRAY
       ns=0; ne=ndata; allocate(delt(ns:ne))
       fctr=half/(times(ne)-times(ns))
       delt(ns)=fctr*(times(ns+1)-times(ns)); 
       delt(ne)=fctr*(times(ne)-times(ne-1))
    do n=ns+1,ne-1
       delt(n)=fctr*(times(n+1)-times(n-1))
    end do
    do n=0,ndata
       if (myid==0) then
          write(*,"(f5.1,'% done')") real(n)*100.0e0/real(ndata)
       end if
       call getCp(n)
       de(:,1)=(qo(:,1)-qa(:,1))
       qb(:,1)=qb(:,1)+delt(n)*de(:,1)*de(:,1)
    end do
    qo(:,1)=sqrt(qb(:,1))
   cinput='Cp'; call wffile(cinput,ndata+2,1)
   !if (myid==7) then
   !   varr(:)=qo(:,1)
   !   call wavg('CpAVG',dir=2,wall=.false.)
   !end if
end if

!==COMPUTE RMS Cf
!qb(:,1)=0
!if (wflag) then
!   call integCoef
!   call gettw(ndata+1)
!   ra0=two/(amachoo**2)
!   do n = 0, lcwall; l=lwall(n)
!      ra1=DOT_PRODUCT(tw(n,:),wtan(n,:))
!      qa(l,1)=ra1*ra0
!   end do
!    do n=0,ndata
!       if (myid==7) then
!          write(*,"(f5.1,'% done')") real(n)*100.0e0/real(ndata)
!       end if
!       call gettw(n)
!       ra0=two/(amachoo**2)
!       do m = 0, lcwall; l=lwall(n)
!          ra1=DOT_PRODUCT(tw(m,:),wtan(m,:))
!          qo(l,1)=ra1*ra0
!       end do
!       qa(:,2)=(qo(:,1)-qa(:,1))
!       qb(:,1)=qb(:,1)+delt(n)*qa(:,2)*qa(:,2)
!    end do
!    qo(:,1)=sqrt(qb(:,1))
!end if
!cinput='Cf'; call wffile(cinput,ndata+2,1)
!!qo(:,1)=qa(:,1)
!!cinput='Cf'; call wffile(cinput,ndata+1,1)
!

!==Extract strip over time
if (fstrip) then
  cinput='out/pt'//cnzone//'.dat'
  if(mb==7) open (unit=7, file=cinput)
  do n = 0, ndata
    if (wflag) then
    call rdP3dS(n,fmblk)
       if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
       do m = 0, lcwall; l=lwall(m)
          wvarr(m)=qo(l,5)
       end do
       if(.not.allocated(svarr)) allocate(svarr(0:lze,0:ndata))
       do k = 0, lze; l=indx2(50,k,1)
          svarr(k,n)=wvarr(l)
       end do
       write(7,"(102(es15.7))") times(n),svarr(:,n)
    end if
  end do
  if(mb==7) close(7)
end if

!==== SHIFT RESTART SOLUTION
fflag=.true.
call rdRsta
m=25
do nn = 1, 5
   do k = 0, lze
      kk=k+m
      if (kk>lze) then
         kk=k+m-lze-1
      end if
      do j = 0, let; l=indx3(0,j,kk,1); ll=indx3(0,j,k,1)
         qo(ll:ll+lxi,nn)=qa(l:l+lxi,nn)
      end do
   end do
   qa(:,nn)=qo(:,nn)
end do
ndati=ndata
call wrRsta

!===== END OF JOB
 CALL MPI_BARRIER(icom,ierr)
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
