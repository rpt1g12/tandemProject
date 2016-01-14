program test
   implicit none
   
 integer,parameter :: nr=kind(0.0d0),mav=0,mfilt=1,nvar=101
 character(*),parameter :: csgnl='out/twt007.dat'
 !character(*),parameter :: csgnl='out/sinusoid.dat'
 character(*),parameter :: cspct='out/spectra.dat'

 real(nr),parameter :: pi=3.141592653589793_nr,half=0.5_nr,amach=0.4_nr

 integer :: i,is,ie,j,k,lmt,m,mwin,stat,nslc,slc,lmf

 real(nr),dimension(:,:),allocatable :: vart,vdt
 real(nr),dimension(:),allocatable :: time,dt,tke,t,winf

 real(nr),dimension(:,:),allocatable :: varf
 real(nr),dimension(:),allocatable :: freq
 real(nr),dimension(:),allocatable :: filt
 real(nr),dimension(nvar) :: vmean,varr,vari
 real(nr) :: period,tpop,fctr,fnt,cosf,sinf,ra0,ra1,ra2,ra3,res,fmax,avgp,avgdt
 real(nr) :: ovlp,ofst,awin
 integer :: iovlp,iofst

 ovlp=0.1_nr
 ofst=1-ovlp

 open(0,file=csgnl)
 lmt=-2; stat=1
 do while (stat.ge.0)
    read(0,*,IOSTAT=stat) 
    lmt=lmt+1
 end do
 rewind(0)

    stat=0
    mwin=min(abs(mav),1); nslc=max(mav,1)
    ra0=(nslc-1)*ofst+1.0_nr
    slc=int(real(lmt,nr)/ra0)
    ra1=real(lmt,nr)/real(slc,nr)
    ofst=(ra1-1)/(nslc-1)
    if (nslc==1) ofst=0
    do while(ofst>1)
       slc=slc+1
       ra1=real(lmt,nr)/real(slc,nr)
       ofst=(ra1-1)/(nslc-1)
       ra0=(nslc-1)*ofst+1.0_nr
       stat=1
    end do
    ovlp=1-ofst
    iofst=int(ofst*slc); iovlp=slc-iofst
    if (stat==1) then
       write(*,"('The overlap has been modified to:',f5.2)") ovlp
    end if

    lmf=slc/2

 allocate(vart(0:lmt,nvar),time(0:lmt))
 allocate(vdt(0:slc,nvar),dt(0:slc),t(0:slc),winf(0:slc))
 allocate(varf(0:lmf,nvar),freq(0:lmf),filt(-3:lmf+3))

 do i = 0, lmt
    read(0,*) time(i),vart(i,:)
 end do
 close(0)

    !vart(:,1:nvar-1)=vart(:,1:nvar-1)/amach
    !time(:)=time(:)*amach

 do i = 1, nvar
   vmean(i)=avg(vart(:,i),time(:))
   vart(:,i)=vart(:,i)-vmean(i)
 end do

 ra0=(time(2)-time(1)+time(lmt)-time(lmt-1))
 do i = 1, lmt-1
   ra0=ra0+half*(time(i+1)-time(i-1))
 end do
 ra0=ra0/(lmt+1)
 ra1=1.0e0/(time(slc)-time(0))

 fmax=half*(1.0/ra0)
 write(*,"('========================================================')") 
 write(*,"(a,x,f5.2,x,a,x,f7.2,x,a,x,f5.2,x,a,f5.2,x,a,f5.2)") &
 't0=',time(0),'tN=',time(lmt),'fmax=',fmax,'fmin=',ra1,'dtAVG=',ra0
 write(*,"(a,x,i4,x,a,x,i4,x,a,i4)") &
 '#Samp=',lmt,'#Var=',nvar,'#freq=',lmf
 write(*,"(a,x,i4,x,a,x,i4,x,a,x,f5.2,x,a,x,f5.2)") &
 '#Wind=',nslc,'#SampWin=',slc,'Overlap=',ovlp,'Offset=',ofst
 write(*,"('========================================================')") 


 !vart(:,nvar)=((vart(:,1)-vmean(1))**2+&
 !       (vart(:,2)-vmean(2))**2+&
 !       (vart(:,3)-vmean(3))**2)*half
 
  varf(:,:)=0
 avgp=0.0_nr
 awin=0.5e0
 if(mwin==0) awin=1.0e0
 do m=1,nslc
       write(*,"('Averaging: ',i2,' of ',i2)") m,nslc
       is=(m-1)*iofst
       ie=is+slc
       period=time(ie)-time(is); tpop=2*pi/period; fctr=half*(time(is)+time(ie))
       avgp=avgp+period
    do i=0,slc
       t(i)=time(i+is)-time(is)-fctr
       winf(i)=window(i,slc,awin)
    end do
       dt(0)=half*(t(1)-t(0)); dt(slc)=half*(t(slc)-t(slc-1))
    do i=1,slc-1
       dt(i)=half*(t(i+1)-t(i-1))
    end do
    do k=1,nvar
       vdt(:,k)=dt(:)*vart(is:ie,k)*winf(:)
       !vmean(k)=sum(vdt(:,k))/period
       !vdt(:,k)=vdt(:,k)-vmean(k)*dt(:)
    end do
 
       ra1=fmax*2*pi-tpop; fctr=(ra1)/(lmf); ra0=tpop
    do j=0,lmf
          freq(j)=j*fctr+ra0; varr(:)=0; vari(:)=0
       do i=0,slc
          fnt=freq(j)*time(i); cosf=cos(fnt); sinf=sin(fnt)
          varr(:)=varr(:)+cosf*vdt(i,:)
          vari(:)=vari(:)-sinf*vdt(i,:)
       end do
          varf(j,:)=varf(j,:)+varr(:)**2+vari(:)**2
    end do
 end do
 
 !avgp=avgp/nslc; tpop=2*pi/avgp;
 avgp=avgp**2

 fctr=sum(winf(:)*winf(:)*dt(:))/period
 write(*,*) fctr
 fctr=(2/(avgp*fctr)); varf(:,:)=fctr*varf(:,:)

 is=0; ie=lmf
 ra0=half; ra1=9.0_nr/32; ra2=0; ra3=-1.0_nr/32
 do k=1,nvar
    do m=1,mfilt
       filt(is:ie)=varf(is:ie,k)
       filt(is-(/1,2,3/))=varf(is,k); filt(ie+(/1,2,3/))=varf(ie,k)
       do i=is,ie
          varf(i,k)=ra0*filt(i)&
          +ra1*(filt(i-1)+filt(i+1))&
          +ra2*(filt(i-2)+filt(i+2))&
          +ra3*(filt(i-3)+filt(i+3))
       end do
    end do
 end do


    ra0=half/(pi); 
    open(1,file=cspct);
    close(1,status='delete')
    open(1,file=cspct)
 do j=0,lmf
       write(1,"(102es15.7)") ra0*freq(j),varf(j,:)
 end do
    close(1)

  contains

  function avg(x,t) result(r)
    implicit none 
    real(nr) :: r
    real(nr) :: fctr
    real(nr), dimension(0:lmt), intent(in) :: x,t
    real(nr), dimension(0:lmt) :: delt
    integer :: n,i

    n=lmt
    fctr=half/(t(n)-t(0))
    delt(0)=fctr*(t(1)-t(0))
    delt(n)=fctr*(t(n)-t(n-1))
    do i = 1, n-1
       delt(i)=fctr*(t(i+1)-t(i-1))
    end do
    r=0_nr
    do i = 0, n
       r=r+delt(i)*x(i)
    end do
  end function avg

  function window(l,n,a) result(r)
     implicit none
     real(nr) :: r
     integer, intent(in) :: l,n
     real(nr) :: a,b,f
     b=1.0e0-a
     f=real(l)/real(n)
     r=a-b*cos(2*pi*f)

  end function window

end program test
