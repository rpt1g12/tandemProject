!*****
!***** 3D PARALLEL SOLVER
!*****

 program mainpost

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 use rpt
 use subsets
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

    ll=max(npro,12); allocate(ista(MPI_STATUS_SIZE,ll),ireq(ll))

    inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll
    inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll
call setup
call ssSetUp
!call ssCheck
call flst(fmblk)
   
!! RPT-READ X,Y,Z COORDINATES
call rdP3dG(fmblk)
do nss = 1, tss
   call wrP3dG_ss(fmblk,nss)
end do
ispost=.true.
!
call getMetrics
CALL MPI_BARRIER(icom,ierr)
if (intgflag) then
   write(*,"('Node ',i2,' prepared to integrate over an area of',f7.3,&
              ' using ',i4,' elements')") &
              myid,sum(aintg),lcintg
end if

!do n = 0, ndata
!   call rdP3dP(n,fmblk)
!   do nss = 1, tss
!      call wrP3dP_ss(n,fmblk,nss=nss)
!   end do
!end do

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
if (favg==1) then
   call p3dStats
end if

!===COMPUTE FORCE COEFFICIENT
if (fcoef==1) then
   if (myid==0) then
      open (unit=17, file='out/signalout0.dat')
      open (unit=18, file='out/signalout1.dat')
      write(17,"(3x,'n',8x,'time',9x,'Clp',9x,'Cdv',5x)")  
      write(18,"(3x,'n',8x,'time',9x,'Clv',9x,'Cdv',5x)")  
      write(*,"(3x,'n',8x,'time',9x,'Cdp',9x,'Cdv',5x)")  
   end if
   do n = 0, ndata
      call clPVpost(nvar=n)
      if (myid==0) then
         ra0=aoa*pi/180;ra1=cos(ra0);ra2=sin(ra0); 
         if(n>ndata) then
            ra3=(-1)
         else
            ra3=times(n)
         end if
         write(17,"(i8,f12.5,f12.7,f12.7)") &
         n,ra3,cl(1,2)*ra1-cl(1,1)*ra2,cl(1,2)*ra2+cl(1,1)*ra1
         write(18,"(i8,f12.5,f12.7,f12.7)") &
         n,ra3,cl(2,2)*ra1-cl(2,1)*ra2,cl(2,2)*ra2+cl(2,1)*ra1
         write(*,"(i8,f12.5,f12.7,f12.7)") &
         n,ra3,cl(1,2)*ra2+cl(1,1)*ra1,cl(2,2)*ra2+cl(2,1)*ra1
      end if
   end do
   if(myid==0) close(17)
   if(myid==0) close(18)
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
         call rdP3dS(n,fmblk)
         res=(qo(l,2)-qa(l,2))
         write(7,"(f10.5,f10.5)") times(n),res
      end do
      close(7)
   end if
end if


!==COMUPTE Q+W+DELTA
if (fcurl==1) then
   ! Get coefficients for time averaging
   call integCoef
   ! Store averaged quantity in qa and qb
   qb(:,:)=0
   qa(:,:)=0
   do n = 0, ndata
      !n=ndata+1
      call getAllDs(n)
      !call wrP3dP(n,fmblk,'Q+W')
      qb(:,:)=qb(:,:)+delt(n)*qo(:,:)
      do nss = 2, tss
         call wrP3dP_ss(n,fmblk,cname='Q+W+NablU',nss=nss)
      end do
      if (fwss==1) then
         ra0=two/(amachoo**2)
         qo(:,1)=0
         do i = 2, 4
            qo(:,i)=0
            if (wflag) then
               do m = 0, lcwall; l=lwall(m)
                  qo(l,i)=tw(m,i-1)
                  if (i==2) then
                  ra1=DOT_PRODUCT(tw(m,:),wtan(m,:))
                  qo(l,1)=ra1*ra0
                  end if
               end do
            end if
         end do
         qo(:,5)=(p(:)-poo)*ra0
         qa(:,:)=qa(:,:)+delt(n)*qo(:,:)
         !call wrP3dP(n,fmblk,'CftwCp')
         call wrP3dP_ss(n,fmblk,cname='Cf+tw+Cp',nss=1)
      end if
   end do
      n=ndata+1
      ! Save averaged data
      qo(:,:)=qb(:,:)
      do nss = 2, tss
         call wrP3dP_ss(n,fmblk,cname='avgQ+W+DELTA',nss=nss)
      end do
      if (fwss==1) then
         qo(:,:)=qa(:,:)
         call wrP3dP_ss(n,fmblk,cname='avgCf+tw+Cp',nss=1)
      end if
end if

!==INTEGRATION
!if (fcurl==1) then
!   !if (allocated(fout)) deallocate(fout)
!   !if (.not.allocated(fout)) allocate(fout(0:lmx,1))
!   !if (mb==7) then
!   if(myid==7.and.fintg) open (unit=7, file='out/intp1.dat')
!   !if(myid==7) open (unit=17, file='out/p1.dat')
!   !if(myid==7) open (unit=18, file='out/p2.dat')
!   do n = 0, ndata
!      !call qcriterion(n);
!      !call getCurl(n);
!      call rdP3dS(n,fmblk)
!      !call rdP3dP(n,fmblk,'Q+W')
!      if(intgflag) varr(:)=qo(:,5)
!      call integrate
!      if(myid==7.and.fintg) write(7,"(es15.7,x,es15.7)") times(n),ra0
!      !if (myid==7) then
!      !   ra0=0;ra1=0
!      !   do k = 0, 50
!      !      do i = 0, 100; ll=indx2(i,k,1); l=lwall(ll)
!      !         ra0=ra0+area(ll)
!      !         ra1=ra1+p(l)*wnor(ll,2)*area(ll)
!      !      end do
!      !   end do
!      !   write(17,"(es15.7,x,es15.7)") timo,ra1/ra0
!      !   ra0=0;ra1=0
!      !   do k = 50, 100
!      !      do i = 0, 100; ll=indx2(i,k,1); l=lwall(ll)
!      !         ra0=ra0+area(ll)
!      !         ra1=ra1+p(l)*wnor(ll,2)*area(ll)
!      !      end do
!      !   end do
!      !   write(18,"(es15.7,x,es15.7)") timo,ra1/ra0
!      !end if
!      !fout(:,1:3)=qo(:,2:4)
!      !call wrP3dF('Omega',n,3,fmblk)
!   end do
!   if(myid==7.and.fintg) close(7)
!   !if (myid==7) then
!   !   close(17);close(18)
!   !end if
!   !end if
!end if

!==Compute wall shear stress & WRITE Cf+WSS+Cp
if (fwss==1) then
   n=ndata+1
   call getAllDs(n)
   ra0=two/(amachoo**2)
   qo(:,1)=0
   do i = 2, 4
      qo(:,i)=0
      if (wflag) then
         do m = 0, lcwall; l=lwall(m)
            qo(l,i)=tw(m,i-1)
            if (i==2) then
            ra1=DOT_PRODUCT(tw(m,:),wtan(m,:))
            qo(l,1)=ra1*ra0
            end if
         end do
      end if
   end do
   qo(:,5)=(p(:)-poo)*ra0
   call wrP3dP(n+1,fmblk,'CftwCp')
   !!!!! This only works with one processor per block
   if (fcf==1) then
   !==COMPUTE Span-Average Cf 
   idum=4 ! Use only this block
      if (myid==idum) then ! use only this block
         open(7,file='out/allCfAVG.dat')
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
   !==COMPUTE Span-averaged WALL DISTANCES
   if (fwplus==1) then
      if (myid==idum) then ! Use only this block
         open(7,file='out/wplus.dat')
         write(7,"('x y z x+ y+ z+')") 
         do nn = 0, lcwall;l=lwall(nn)
         write(7,"(f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)")&
               xyz(l,1),xyz(l,2),xyz(l,3),wplus(nn,1),wplus(nn,2),wplus(nn,3)
         end do
         close(7)
         if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
         wvarr=wplus(:,1)
         call wavg('x+',dir=2,wall=.true.)
         wvarr=wplus(:,2)
         call wavg('y+',dir=2,wall=.true.)
         wvarr=wplus(:,3)
         call wavg('z+',dir=2,wall=.true.)
      end if
   end if
end if

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

if (fprobcirc==1) then
   call probCirc(ndata)
end if

if (fijkmax==1) then
   call rdP3dP(ndata+2,fmblk)
   varr=qo(:,5)
   call getijkMax(4,(/9,4/),(/(i,i=60,140,10)/),(/12,37,62,87/))
   call getvalMax(4,(/9,4/),ndata)
end if
   

!===== END OF JOB
 CALL MPI_BARRIER(icom,ierr)
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
