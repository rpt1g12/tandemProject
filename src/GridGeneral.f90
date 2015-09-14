!*****  
!***** 3D Generic-Grid
!*****

module gridgen

 use subroutineso
 implicit none

 integer(k4),parameter :: lnaca=1000
 real(k8),parameter :: pi4=pi/4.0_k8
 integer(k4) :: lxit,lett,lxisz,im,jm
 integer(k4) :: lxis,lxie,lxib
 integer(k4) :: lets,lete,letb

 real(k8),dimension(0:lnaca,2:3,2) :: xnaca,ynaca

 real(k8),dimension(0:5) :: xa,xb,xc,xd,xe,xf
 real(k8),dimension(0:3) :: ya,yb,yc,yd,ye,yf,yg
 real(k8),dimension(0:5,0:2,0:1) :: hslo
 real(k8),dimension(0:3,0:3,0:1) :: vslo

 real(k8),dimension(:,:),allocatable :: xx,yy,zz
 real(k8),dimension(:),allocatable :: zs
 real(k8),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(k8),dimension(:),allocatable :: pxi,qet

 real(k8) :: rs,re,rp,ts,te,shs1,she1,shs2,she2,shs,she,shswle
 real(k8) :: xo,xjct,yo,yjct,sho,pp,qq
 real(k8) :: am,err,tmp,tmpa,tmpb,gf
 real(k8) :: k01,k02,k03,k04,x0,x1
 real(k8) :: deg1,deg2

 !allocate(lxise(0:bkx-1,0:1),letse(0:bky-1,0:1),lzese(0:bkz-1,0:1))

 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,litr,smgrid,&
            domlen,span,wlew,wlea,szth1,szth2,szxt,&
            tla,tlb,cutlb,c1,delt1,ximod,etamod)

 integer(k4),intent(in) :: ngridv,litr,nthick
 integer(k4), dimension(:,:) :: lxise(0:bkx-1,0:1)
 integer(k4), dimension(:,:) :: letse(0:bky-1,0:1)
 integer(k4), dimension(:,:) :: lzese(0:bkz-1,0:1)
 real(k8),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,c1,delt1
 real(k8),intent(in) :: ximod,etamod,cutlb
 real(k8) :: lsz1,lsz2
 real(k8) :: lbl,lwle
 real(k8) :: alph
 real(k8) :: oxp,oyp
 real(k8) :: tmps,tmpe,tmpc
 real(k8) :: sha,shb,shc
 integer(k4) :: smod
 logical :: flag

    lxit=sum(lxibk(:))+(bkx-1); lett=sum(letbk(:))+(bky-1)
        lxise(0,0)=0; lxise(0,1)=lxibk(0)
    letse(0,0)=0; letse(0,1)=letbk(0)
    lzese(0,0)=0; lzese(0,1)=lzebk(0)
    do i = 1, bkx-1
          lxise(i,0)=lxise(i-1,1)+1
          lxise(i,1)=lxise(i,0)+lxibk(i)
    end do
    do j = 1, bky-1
          letse(j,0)=letse(j-1,1)+1
          letse(j,1)=letse(j,0)+letbk(j)
    end do

    shs=smgrid; she=shs
    shs1=ximod*smgrid; she1=shs1
    shs2=etamod*smgrid;
    smod=2
    tmp=(shs2+smod*shs2)*half
    lbl=max(tmp*letbk(1)/(sin(pi4)),0.10*c1/(sin(pi4)))
    lbl=min(lbl,0.65*c1/(sin(pi4)))
    if (myid==0) then
       write(*,*) lbl
    end if

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lzebk(0)))
    allocate(xp(0:lxit,0:5),yp(0:lxit,0:5))
    allocate(xq(0:lett,0:3),yq(0:lett,0:3))
    allocate(pxi(0:lxit),qet(0:lett))

if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream',form='unformatted')

 do k=0,lzebk(0)
    zs(k)=span*(real(lzebk(0)-k,k8)/lzebk(0)-half)
!---BLOCKS' BOUNDARIES
    sho=tla/litr; ll=2*litr; lsz1=ll*sho; lsz2=szth1+szxt
!---WAVY LEADING-EDGE PROFILE
    lwle=wlea*sin(2*pi*(zs(k)-zs(0))/wlew)
!---HORIZONTAL LINES
    xa(:)=-domlen;
    xb(:)=xa+lsz1
    xc(2:3)=-half*c1+lwle*cos(delt1);
    xc(0:1)=xc(2)-lbl*cos(pi4-delt1);xc(4:5)=xc(3)-lbl*cos(pi4+delt1)
    if (nthick*(nthick-3)==0) then
       xc(0:1)=xc(2);xc(4:5)=xc(3)
    end if
    xd(:)=-half*c1+c1*cos(delt1)
    xe(:)=domlen-szth1
    xf(:)=xe(:)+lsz2
!---VERTICAL LINES
    ya(:)=-domlen;
    yb=ya(:)+lsz1;
    yd(0)=zero;yd(3)=zero
    yd(1)=-lwle*sin(delt1)
    yd(2)=-c1*sin(delt1)
    ye(0:1)=yd(0:1)+lbl*sin(pi4+delt1)
    ye(2)=yd(2)+lbl*sin(pi4)*cos(delt1)
    ye(3)=ye(0)
        yc(0)=ye(0)-2*lbl*sin(pi4+delt1)
    yc(1)=yd(1)-lbl*sin(pi4-delt1)
    yc(2)=yd(2)-lbl*sin(pi4)*cos(delt1)
    yc(3)=yc(0)
    yg(:)=domlen
    yf(:)=yg(:)-lsz1

    fctr=2*pi/wlew; shswle=shs*sqrt(1+0*(fctr*wlea*cos(fctr*(zs(k)-zs(0))))**2)

!----- INITIAL AND END HORIZONTAL SLOPES
    deg1=(25_k8*pi/180_k8)
    deg2=(10_k8*pi/180_k8)
    hslo(0,:,:)=zero
    hslo(1,0,:)=(/zero,-tan(delt1)/)
    hslo(1,1,:)=(/tan(-deg1-delt1),tan(deg2-delt1)/)
    hslo(1,2,:)=(/-tan(delt1),zero/)
    hslo(4,0,:)=(/zero,-tan(delt1)/)
    hslo(4,1,:)=(/tan(deg1-delt1),tan(-deg2-delt1)/)
    hslo(4,2,:)=(/-tan(delt1),zero/)
    hslo(5,:,:)=zero
    hslo(2,0,:)=(/zero,-tan(delt1)/)
    hslo(2,1,:)=(/zero,zero/)
    hslo(2,2,:)=(/-tan(delt1),zero/)
    hslo(3,0,:)=hslo(2,0,:)
    hslo(3,1,:)=hslo(2,1,:)
    hslo(3,2,:)=hslo(2,2,:)
    select case (nthick);
    case(0);
    hslo(1,1,:)=zero
    hslo(4,1,:)=zero
    case(3);
    hslo(1,1,:)=zero
    hslo(4,1,:)=zero
    end select

!----- INITIAL AND END VERTICAL SLOPES
    vslo(0,:,:)=zero
    vslo(1,0,:)=zero
    vslo(1,1,:)= (/zero,tan(pi4+delt1)/)
    vslo(1,2,:)=(/-tan(pi4-delt1),zero/)
    vslo(1,3,:)=zero
    vslo(2,:,:)=zero
    vslo(3,:,:)=zero
    select case (nthick);
    case(0)
    vslo(1,1,:)=zero
    vslo(1,2,:)=zero
    vslo(3,1,:)=zero
    vslo(3,2,:)=zero
    case(2);
    vslo(3,1,:)=zero
    vslo(3,2,:)=zero
    case(3);
    vslo(1,1,:)=zero
    vslo(1,2,:)=zero
    end select    
    !----- AEROFOIL SURFACE GRID POINTS
    m=1
    tmp=c1-lwle;tmpa=xc(2);tmpb=yd(1);lxis=lxise(1,0);lxie=lxise(1,1);lxib=lxibk(1);alph=delt1
    if (nthick*(nthick-3)==0) then
       flag=.false.
    else
       flag=.true.
    end if
    if (flag) then
       ! READ COORDINATES FROM FILE
       open(2,file='aerofoil.dat')
       do n=2,3; do i=0,lnaca
          read(2,*) xnaca(i,n,m),ynaca(i,n,m)
          xnaca(i,n,m)=tmp*xnaca(i,n,m); ynaca(i,n,m)=tmp*ynaca(i,n,m)
       end do; end do
       close(2)
       ! DETERMINE UPPER AND LOWER SIDES
       do n=2,3
          yp(lxis,n)=zero;xp(lxis,n)=zero
          ! DETERMINE THE FIRST LL POINTS
          ll=8 ! "LL" MUST BE EQUAL TO OR LARGER THAN 4.
          do i=lxis+1,lxis+ll
             xp(i,n)=xp(i-1,n)+half*shs1; err=1
             do while(abs(err)>sml)
                yp(i,n)=naca(xp(i,n),tmp,21.0_k8,n)!ylagi(i,n,m)
                err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shs1-1;
                xp(i,n)=xp(i,n)-half**5*err*shs1
             end do
          end do
          xo=xp(lxis+ll,n); sho=sum(xp(lxis+ll-4:lxis+ll,n)*(/3,-16,36,-48,25/))/12
          ! COMPUTE THE REST OF THE POINTS
          ip=lxis+ll;im=lxib-ll 
          call gridf(xp(:,n),pxi,xo,tmp,sho,she1,lxit,im,ip)
          do i=lxis+ll+1,lxie-1
             yp(i,n)=naca(xp(i,n),tmp,21.0_k8,n)!ylagi(i,n,m)
          end do
          yp(lxie,n)=zero
          ! ROTATE AND MOVE AEROFOILS
          do i = lxis, lxie
          oxp=xp(i,n);oyp=yp(i,n)
          xp(i,n) = (oxp*cos(alph)+oyp*sin(alph))+tmpa;
          yp(i,n) = (-oxp*sin(alph)+oyp*cos(alph))+tmpb;
          end do
       end do
    else
      ip=lxis; im=lxib;
      call gridf(xp(:,2),pxi,zero,tmp,shs1,she1,lxit,im,ip)
      yp(lxis:lxie,2)=zero
      xp(lxis:lxie,3)=xp(lxis:lxie,2)
      yp(lxis:lxie,3)=yp(lxis:lxie,2)
      do n = 2, 3
      ! ROTATE AND MOVE AEROFOILS
        do i = lxis, lxie
        oxp=xp(i,n);oyp=yp(i,n)
        xp(i,n) = (oxp*cos(alph)+oyp*sin(alph))+tmpa;     
                yp(i,n) = (-oxp*sin(alph)+oyp*cos(alph))+tmpb;
        end do
      end do
    end if

!--HORIZONTAL INTERFACES
   do n = 0,4,2
   !--X-DIRECTION
   !--BLOCK0
   !--a-b
      sho=tla/litr; ll=2*litr
      ip=lxise(0,0); im=ll;
      tmpa=xa(n);sha=sho;tmpb=xb(n);shb=sho
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--b-c
      ip=ip+im; im=lxibk(0)-ll;
      tmpa=xb(n);sha=sho;tmpb=xc(n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
      if(k==0.and.n==0) then
         lxisz=lxibk(2)*(minloc(abs(xa(n)+szth1-xp(0:lxibk(0),n)),1)-1)/lxibk(0); lp=ll     
      end if
   if (n.ne.2) then
   !--BLOCK1
   !--c-d
      ip=lxise(1,0); im=lxibk(1);
      tmpa=xc(n);sha=shs1;tmpb=xd(n);shb=she1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   end if
   !--BLOCK2
   !--d-e
      ip=lxise(2,0); im=lxibk(2)-lxisz;
      tmpa=xd(n);sha=she1;tmpb=xe(n);shb=sml
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--e-f
      ip=ip+im; im=lxisz;
      tmpa=xe(n);sha=pxi(ip);tmpb=xf(n);shb=sho
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--COPY ON N+1 INTERFACE
      xp(lxise(0,0):lxise(0,1),n+1)=xp(lxise(0,0):lxise(0,1),n)
      xp(lxise(2,0):lxise(2,1),n+1)=xp(lxise(2,0):lxise(2,1),n)
      if (n.ne.2) then
      xp(lxise(1,0):lxise(1,1),n+1)=xp(lxise(1,0):lxise(1,1),n)
      end if
   end do

   !--Y-DIRECTION
   !--INTERFACE 0 & 1
   do n = 0,1
      do m = 0,2
      selectcase(n)
      case(0); k01=ya(m);k03=ya(m+1)
      case(1); k01=yc(m);k03=yc(m+1)
      end select
      selectcase(m)
      case(0); is=lxise(0,0);ie=lxise(0,1);x0=xa(n);x1=xc(n)
      case(1); is=lxise(1,0);ie=lxise(1,1);x0=xc(n);x1=xd(n)    
            case(2); is=lxise(2,0);ie=lxise(2,1);x0=xd(n);x1=xf(n)
      end select     
      k02=hslo(n,m,0);k04=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k01,k02,k03,k04,x0,x1,xp(i,n))
      end do
      end do
   end do
   !--INTERFACE 4 & 5
   do n = 4,5
      do m = 0,2
      selectcase(n)
      case(4); k01=ye(m);k03=ye(m+1)
      case(5); k01=yg(m);k03=yg(m+1)
      end select
      selectcase(m)
      case(0); is=lxise(0,0);ie=lxise(0,1);x0=xa(n);x1=xc(n)
      case(1); is=lxise(1,0);ie=lxise(1,1);x0=xc(n);x1=xd(n)
      case(2); is=lxise(2,0);ie=lxise(2,1);x0=xd(n);x1=xf(n)
      end select     
      k02=hslo(n,m,0);k04=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k01,k02,k03,k04,x0,x1,xp(i,n))
      end do
      end do
   end do
   !--INTERFACE 2 & 3
      n=2
      do m = 0,2,2
      selectcase(m)
      case(0); is=lxise(0,0);ie=lxise(0,1);x0=xa(n);x1=xc(n)
      case(2); is=lxise(2,0);ie=lxise(2,1);x0=xd(n);x1=xf(n)
      end select     
      k01=yd(m);k03=yd(m+1);k02=hslo(n,m,0);k04=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k01,k02,k03,k04,x0,x1,xp(i,n))
         yp(i,n+1)=yp(i,n)
      end do
      end do

!--VERICAL END BOUNDARIES
   !-Y-DIRECTION
   do n = 0, 3
         alph=delt1
      !-BLOCK1
      !-c-d
      ip=letse(1,0); im=letbk(1);
      tmpa=yc(n);sha=smod*shs2;tmpb=yd(n);shb=shs2*sin(pi4-alph)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      she2=qet(ip)
      !-BLOCK0
      !-a-b
      sho=tla/litr; ll=lsz1/sho
      ip=letse(0,0); im=ll;
      tmpa=ya(n);sha=sho;tmpb=yb(n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)   
            !-b-c
      ip=ip+im; im=letbk(0)-im;
      tmpa=yb(n);sha=sho;tmpb=yc(n);shb=she2
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-BLOCK2
      !-d-e
      ip=letse(2,0); im=letbk(2);
      tmpa=yd(n);sha=shs2*sin(pi4+alph);tmpb=ye(n);shb=smod*shs2
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      she2=qet(ip+im)
      !-BLOCK3
      !-e-f
      sho=tla/litr; ll=lsz1/sho
      ip=letse(3,0); im=letbk(3)-ll;
      tmpa=ye(n);sha=she2;tmpb=yf(n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-f-g
      ip=ip+im; im=ll;
      tmpa=yf(n);sha=sho;tmpb=yg(n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
   end do
   !-COPY ON N+5
   !yq(:,3)=yq(:,0)

   !-X-DIRECTION
   !-INTERFACE 0
   xq(:,0)=xa(0)
   !-INTERFACE 1
   n=1
   do m = 0,3
      selectcase(m)
      case(0); is=letse(0,0);ie=letse(0,1);x0=ya(n);x1=yc(n)
               k01=xc(m);k03=xc(m+1)
      case(1); is=letse(1,0);ie=letse(1,1);x0=yc(n);x1=yd(n)
               k01=xc(m);k03=xc(m+1)
      case(2); is=letse(2,0);ie=letse(2,1);x0=yd(n);x1=ye(n)
               k01=xc(m+1);k03=xc(m+2)
      case(3); is=letse(3,0);ie=letse(3,1);x0=ye(n);x1=yf(n)
               k01=xc(m+1);k03=xc(m+2)
      end select     
      k02=vslo(n,m,0);k04=vslo(n,m,1)
      do i = is, ie
         xq(i,n)=inter(k01,k02,k03,k04,x0,x1,yq(i,n))
      end do
   end do
   !-INTERFACE 2
   xq(:,2)=xd(0)
   !-INTERFACE 3
   xq(:,3)=xf(0)

!--GRID INTERPOLATION
   do l=0,3
      select case(l);
      case(0); js=0; je=letse(0,1); n=0; gf=0
      case(1); js=letse(1,0); je=letse(1,1); n=1; gf=1
      case(2); js=letse(2,0); je=letse(2,1); n=3; gf=0       
            case(3); js=letse(3,0); je=letse(3,1); n=4; gf=1
      end select
     do m=0,2
      select case(m)
      case(0); is=0; ie=lxise(0,1); nn=1;
               rs=abs(xa(n+1)-xa(n)); re=abs(xc(n+1)-xc(n))
      case(1); is=lxise(1,0); ie=lxise(1,1); nn=0;
               rs=abs(xc(n+1)-xc(n)); re=abs(xd(n+1)-xd(n))
      case(2); is=lxise(2,0); ie=lxise(2,1); nn=1;
               rs=abs(xd(n+1)-xd(n)); re=abs(xe(n+1)-xe(n))
      end select
               tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
      do j=js,je; do i=is,ie
         pp=real(i-is,k8)/(ie-is); qq=real(j-js,k8)/(je-js)
         tmp=sin(half*pi*pp); ra0=(gf*real(je-j,k8)/(je-js)+gf*qq); ra1=one+(2-one)*ra0**2
         pxi(i)=(1-nn)*tmp**ra1+nn*tmp**2
         ts=tmpa*tmpb*(1-pxi(i))+1-tmpb; te=1-ts
         xx(i,j)=(xp(i,n+1)-xp(i,n))*(ts*(xq(j,m)-xq(js,m)+qq*(1-tmpa))/(xq(je,m)-xq(js,m)+1-tmpa)&
                +te*(xq(j,m+1)-xq(js,m+1)+qq*(1-tmpb))/(xq(je,m+1)-xq(js,m+1)+1-tmpb))+xp(i,n)
         ts=1-pxi(i); te=1-ts
         yy(i,j)=(yp(i,n+1)-yp(i,n))*(ts*(yq(j,m)-yq(js,m))/(yq(je,m)-yq(js,m))&
                +te*(yq(j,m+1)-yq(js,m+1))/(yq(je,m+1)-yq(js,m+1)))+yp(i,n)
         zz(i,j)=zs(k)
      end do; end do
     end do
   end do

!--GRID OUTPUT
   select case(mb); 
     case(0,1,2); js=letse(0,0); je=letse(0,1); 
     case(3,4,5); js=letse(1,0); je=letse(1,1);
     case(6,7,8); js=letse(2,0); je=letse(2,1);
     case(9,10,11); js=letse(3,0); je=letse(3,1);
   end select
   select case(mb)
     case(0,3,6,9); is=lxise(0,0); ie=lxise(0,1);
     case(1,4,7,10); is=lxise(1,0); ie=lxise(1,1);
     case(2,5,8,11); is=lxise(2,0); ie=lxise(2,1)
   end select
   np=k8*(je-js+1)*(ie-is+1)
   write(1,pos=np*k+1) ((xx(i,j),i=is,ie),j=js,je)
   write(1,pos=np*(k+lzebk(0)+1)+1) ((yy(i,j),i=is,ie),j=js,je)
   write(1,pos=np*(k+2*lzebk(0)+2)+1) ((zz(i,j),i=is,ie),j=js,je)
 end do
 close(1)

!--GRID CHECKING
   if(ngridv==1) then
   open(1,file='misc/gridview'//cno(2)//cno(1)//cno(0)//'.dat')
    write(1,*) 'variables=x,y,z'; write(1,"('zone i=',i4,' j=',i4,' k=',i4)") ie-is+1,je-js+1,1
   do j=js,je; do i=is,ie
     write(1,'(2es15.7)') xx(i,j),yy(i,j),zero
   end do; end do
   close(1)
   else
   open(1,file='misc/gridview'//cno(2)//cno(1)//cno(0)//'.dat'); close(1,status='delete') 
      end if

   if(myid==0) then
      write(*,"('Grid generation is complete.')")
   end if

end if

   deallocate(xx,yy,zz,zs,xp,yp,xq,yq,pxi,qet)

 end subroutine gridaerofoil

!===== LAGRANGIAN INTERPOLATION FOR Y-COORDINATES

 function ylagi(i,n,m) result(yi)

 integer(k4),intent(in) :: i,n,m
 real(k8) :: yi

    is=0; ie=lnaca; ii=minloc(abs(xnaca(:,n,m)-xp(i,n)),1)-1; ip=ii
 if(ii-is<=1) then; ip=is+2; end if
 if(ie-ii<=1) then; ip=ie-2; end if
    yi=0; alag(:)=xp(i,n)-xnaca(ip-2:ip+2,n,m)
 do jj=-2,2
    blag(:)=xnaca(ip+jj,n,m)-xnaca(ip-2:ip+2,n,m); ao=1; bo=1
 do ii=-2,2; if(ii/=jj) then
    ao=ao*alag(ii); bo=bo*blag(ii)
 end if; end do
    yi=yi+ao*ynaca(ip+jj,n,m)/bo
 end do

 end function ylagi

!===== NACA FUNCTION
 function naca(x,c,t0,n) result(y)
    real(k8),intent(in) :: x,c,t0
    integer(k4), intent(in) :: n
    real(k8) :: y
    real(k8) :: t,k,xc

    t=t0*0.01e0
    k=0.991148635e0
    xc=x/c
    y =(t*c*k/0.2e0) &
    *(0.298222773e0*sqrt(xc)-0.127125232e0*xc- 0.357907906e0*xc**2+&
    0.291984971e0*xc**3-0.105174606e0*xc**4)
    if (n==2) then
       y=-1.0e0*y
    end if
 end function naca
 
!===== FUNCTION FOR INTERFACE POINTS
 function inter(k01,k02,k03,k04,x0,x1,x) result(y)

 real(k8) :: y
 real(k8), intent(in) :: x,x0,x1,k01,k02,k03,k04  
real(k8) :: a,b,c,d,l,invl,fx

 l=x1-x0
 invl=1.0_k8/l
 fx=(x-x0)*invl

    a=(2*fx**3-3*fx**2+1);
    b=(fx**3-2*fx**2+fx);
    c=(-2*fx**3+3*fx**2);
    d=(fx**3-fx**2);
    y=a*k01+b*l*k02+c*k03+d*l*k04;
 
 end function inter
 !=====

end module gridgen

!*****               

