!*****
!***** FOURIER TRANSFORM
!*****

 program fouriertransform

 implicit none

! integer,parameter :: nr=kind(0.0d0),lmf=120,mav=1,mfilt=120,nvar=1+3*26,nvs=1
! character(13),parameter :: csgnl='signalwle.dat'
! character(14),parameter :: cspct='spectrawle.dat'

 integer,parameter :: nr=kind(0.0d0),lmf=64,mav=1,mfilt=1,nvar=3*26,nvs=0
 character(18),parameter :: csgnl='inflowsignal.dat'
 character(19),parameter :: cspct='inflowspectra.dat'

 real(nr),parameter :: pi=3.141592653589793_nr,half=0.5_nr,amach=0.24_nr
 real(nr),parameter :: turbi=1.0_nr*0.025_nr,turbl=0.04_nr,turblo=1.33898*turbl

 integer :: i,is,ie,j,k,lmt,m,mwin

 real(nr),dimension(:,:),allocatable :: vart,vdt
 real(nr),dimension(:),allocatable :: time,dt

 real(nr),dimension(0:lmf,nvar) :: varf
 real(nr),dimension(0:lmf) :: freq,ek11,ek22
 real(nr),dimension(-3:lmf+3) :: filt
 real(nr),dimension(nvar) :: vmean,varr,vari
 real(nr) :: period,tpop,fctr,fnt,cosf,sinf,ra0,ra1,ra2,ra3,res

    inquire(file=csgnl,size=m); lmt=m/(16*(1+nvar))-1
    allocate(vart(0:lmt,nvar),vdt(0:lmt,nvar),time(0:lmt),dt(0:lmt))

    open(0,file=csgnl)
 do i=0,lmt
    read(0,*) time(i),vart(i,:)
 end do
    close(0)

    mwin=min(mav-1,1); varf(:,:)=0
 do m=1,mav
       period=time(lmt)-time(0); tpop=2*pi/period; fctr=half*(time(0)+time(lmt))
    do i=0,lmt
       time(i)=time(i)-fctr
    end do
       dt(0)=half*(time(1)-time(0)); dt(lmt)=half*(time(lmt)-time(lmt-1))
    do i=1,lmt-1
       dt(i)=half*(time(i+1)-time(i-1))
    end do
    do k=1,nvar
       vdt(:,k)=dt(:)*vart(:,k)*(mwin*(sin(tpop*time(:))**2-1)+1)
       vmean(k)=sum(vdt(:,k))/period
       vdt(:,k)=vdt(:,k)-vmean(k)*dt(:)
    end do

       ra0=tpop; ra1=ra0+lmf*tpop; fctr=(ra1-ra0)/lmf
    do j=0,lmf
          freq(j)=j*fctr+ra0; varr(:)=0; vari(:)=0
       do i=0,lmt
          fnt=freq(j)*time(i); cosf=cos(fnt); sinf=sin(fnt)
          varr(:)=varr(:)+cosf*vdt(i,:)
          vari(:)=vari(:)-sinf*vdt(i,:)
       end do
          varf(j,:)=varf(j,:)+varr(:)**2+vari(:)**2
    end do
    if(m==1) then
    do k=1,nvar
       write(*,*) sum(dt(:)*(vart(:,k)-vmean(k))**2)/period,2*sum(varf(:,k))/period**2
    end do
    end if

    if(m/=mav) then
          k=m*lmt/mav
       do i=0,lmt-k
          dt(i)=time(i+k); vdt(i,:)=vart(i+k,:)
       end do
       do i=lmt-k+1,lmt
          dt(i)=time(i+k-lmt)+period; vdt(i,:)=vart(i+k-lmt,:)
       end do
          time(:)=dt(:); vart(:,:)=vdt(:,:)
    end if
 end do
    fctr=2/(mav*period); varf(:,:)=fctr*varf(:,:)

    is=0; ie=lmf
    ra0=half; ra1=9.0_nr/32; ra2=0; ra3=-1.0_nr/32
 do k=1,nvar
    if(k==nvs) then
       varf(:,k)=log10(varf(:,k)); filt(is-(/1,2,3/))=varf(is,k); filt(ie+(/1,2,3/))=varf(ie,k)
    end if
    do m=1,mfilt
          filt(is:ie)=varf(is:ie,k)
       if(k/=nvs) then
          filt(is-(/1,2,3/))=varf(is,k); filt(ie+(/1,2,3/))=varf(ie,k)
       end if
       do i=is,ie
          varf(i,k)=ra0*filt(i)+ra1*(filt(i-1)+filt(i+1))+ra2*(filt(i-2)+filt(i+2))+ra3*(filt(i-3)+filt(i+3))
       end do
    end do
    if(k==nvs) then
       varf(:,k)=10**varf(:,k)
    end if
 end do

    ra0=4*amach*turbi**2*turbl; ra1=turblo/amach
 do i=is,ie
    res=(ra1*freq(i))**2; ra2=1+res; ra3=1/ra2**(11.0_nr/6)
    ek11(i)=ra0*ra2*ra3; ek22(i)=ra0*(3+8*res)*ra3/6
 end do

    ra0=half/pi; ra1=(1.01325e5/(2e-5))**2
    open(0,file=cspct)
 do j=0,lmf
    write(0,'(es15.7)',advance='no') ra0*freq(j)
    write(0,'(es15.7)',advance='no') ek11(j)
    write(0,'(es15.7)',advance='no') ek22(j)
    if(nvs==1) then
   !    write(0,'(es15.7)',advance='no') 10*log10(ra1*abs(varf(j,1)))
       write(0,'(es15.7)',advance='no') varf(j,1)
    end if
    write(0,'(es15.7)',advance='no') 3*sum(varf(j,nvs+1:nvar-2:3))/(nvar-nvs)
    write(0,'(es15.7)',advance='no') 3*sum(varf(j,nvs+2:nvar-1:3))/(nvar-nvs)
    write(0,'(es15.7)') 3*sum(varf(j,nvs+3:nvar:3))/(nvar-nvs)
 end do
    close(0)

 end program fouriertransform

!*****
