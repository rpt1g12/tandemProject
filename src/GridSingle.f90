!*****
!***** 3D SINGLE AEROFOIL GRID GENERATION
!*****

module gridgen

 use subroutineso
 implicit none

 integer,parameter :: lnaca=1000
 real(nr),parameter :: pi4=pi/4.0_nr
 integer :: lxi0,lxi1,lxi2,let0,let1,lze0
 integer :: lxit,lett,lxisz,im,jm
 integer :: lxis0,lxie0,lxis1,lxie1,lxis2,lxie2
 integer :: lets0,lete0,lets1,lete1,lets2,lete2,lets3,lete3
 integer :: lxis,lxie,lxib
 integer :: lets,lete,letb

 real(nr),dimension(0:lnaca,2:3,2) :: xnaca,ynaca

 real(nr),dimension(0:5) :: xa,xb,xc,xd,xe,xf
 real(nr),dimension(0:3) :: ya,yb,yc,yd,ye,yf,yg
 real(nr),dimension(0:5,0:2,0:1) :: hslo
 real(nr),dimension(0:3,0:3,0:1) :: vslo

 real(nr),dimension(:,:),allocatable :: xx,yy,zz
 real(nr),dimension(:),allocatable :: zs
 real(nr),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(nr),dimension(:),allocatable :: pxi,qet

 real(nr) :: rs,re,rp,ts,te,shs1,she1,shs2,she2,shs,she,shswle
 real(nr) :: xo,xjct,yo,yjct,sho,pp,qq
 real(nr) :: am,err,tmp,tmpa,tmpb,gf
 real(nr) :: k1,k2,k3,k4,x0,x1
 real(nr) :: deg1,deg2

 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,litr,smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb,c1,delt1,ximod,etamod)

 integer,intent(in) :: ngridv,litr,nthick
 real(nr),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,c1,delt1
 real(nr),intent(in) :: ximod,etamod,cutlb
 real(nr) :: lsz1,lsz2
 real(nr) :: lbl,lwle
 real(nr) :: alph
 real(nr) :: oxp,oyp
 real(nr) :: tmps,tmpe,tmpc
 real(nr) :: sha,shb,shc
 integer :: smod
 logical :: flag

    lxit=lxi0+lxi1+lxi2+2; lett=2*(let0+let1)+3

    lxis0=0;lxie0=lxi0;
    lxis1=lxie0+1;lxie1=lxis1+lxi1;
    lxis2=lxie1+1;lxie2=lxit;

    lets0=0;lete0=let0;
    lets1=lete0+1;lete1=lets1+let1
    lets2=lete1+1;lete2=lets2+let1
    lets3=lete2+1;lete3=lets3+let0

    shs=smgrid; she=shs
    shs1=ximod*smgrid; she1=shs1
    shs2=etamod*smgrid;
    smod=5
    tmp=(shs2+smod*shs2)*half
    !lbl=max(tmp*let1/(sin(pi4)),0.10*c1/(sin(pi4)))
    lbl=max(tmp*80/(sin(pi4)),0.10*c1/(sin(pi4)))

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lze0))
    allocate(xp(0:lxit,0:5),yp(0:lxit,0:5))
    allocate(xq(0:lett,0:3),yq(0:lett,0:3))
    allocate(pxi(0:lxit),qet(0:lett))

if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream')

 do k=0,lze0
    zs(k)=span*(real(lze0-k,nr)/lze0-half)
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
    deg1=(25_nr*pi/180_nr)
    deg2=(10_nr*pi/180_nr)
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
    tmp=c1-lwle;tmpa=xc(2);tmpb=yd(1);lxis=lxis1;lxie=lxie1;lxib=lxi1;alph=delt1
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
                yp(i,n)=naca(xp(i,n),tmp,21.0_nr,n)!ylagi(i,n,m)
                err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shs1-1;
                xp(i,n)=xp(i,n)-half**5*err*shs1
             end do
          end do
          xo=xp(lxis+ll,n); sho=sum(xp(lxis+ll-4:lxis+ll,n)*(/3,-16,36,-48,25/))/12
          ! COMPUTE THE REST OF THE POINTS
          ip=lxis+ll;im=lxib-ll 
          call gridf(xp(:,n),pxi,xo,tmp,sho,she1,lxit,im,ip)
          do i=lxis+ll+1,lxie-1
             yp(i,n)=naca(xp(i,n),tmp,21.0_nr,n)!ylagi(i,n,m)
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
      ip=lxis0; im=ll;
      tmpa=xa(n);sha=sho;tmpb=xb(n);shb=sho
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--b-c
      ip=ip+im; im=lxi0-ll;
      tmpa=xb(n);sha=sho;tmpb=xc(n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
      if(k==0.and.n==0) then
         lxisz=lxi2*(minloc(abs(xa(n)+szth1-xp(0:lxi0,n)),1)-1)/lxi0; lp=ll     
      end if
   if (n.ne.2) then
   !--BLOCK1
   !--c-d
      ip=lxis1; im=lxi1;
      tmpa=xc(n);sha=shs1;tmpb=xd(n);shb=she1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   end if
   !--BLOCK2
   !--d-e
      ip=lxis2; im=lxi2-lxisz;
      tmpa=xd(n);sha=she1;tmpb=xe(n);shb=sml
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--e-f
      ip=ip+im; im=lxisz;
      tmpa=xe(n);sha=pxi(ip);tmpb=xf(n);shb=sho
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--COPY ON N+1 INTERFACE
      xp(lxis0:lxie0,n+1)=xp(lxis0:lxie0,n)
      xp(lxis2:lxie2,n+1)=xp(lxis2:lxie2,n)
      if (n.ne.2) then
      xp(lxis1:lxie1,n+1)=xp(lxis1:lxie1,n)
      end if
   end do

   !--Y-DIRECTION
   !--INTERFACE 0 & 1
   do n = 0,1
      do m = 0,2
      selectcase(n)
      case(0); k1=ya(m);k3=ya(m+1)
      case(1); k1=yc(m);k3=yc(m+1)
      end select
      selectcase(m)
      case(0); is=lxis0;ie=lxie0;x0=xa(n);x1=xc(n)
      case(1); is=lxis1;ie=lxie1;x0=xc(n);x1=xd(n)
      case(2); is=lxis2;ie=lxie2;x0=xd(n);x1=xf(n)
      end select     
      k2=hslo(n,m,0);k4=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k1,k2,k3,k4,x0,x1,xp(i,n))
      end do
      end do
   end do
   !--INTERFACE 4 & 5
   do n = 4,5
      do m = 0,2
      selectcase(n)
      case(4); k1=ye(m);k3=ye(m+1)
      case(5); k1=yg(m);k3=yg(m+1)
      end select
      selectcase(m)
      case(0); is=lxis0;ie=lxie0;x0=xa(n);x1=xc(n)
      case(1); is=lxis1;ie=lxie1;x0=xc(n);x1=xd(n)
      case(2); is=lxis2;ie=lxie2;x0=xd(n);x1=xf(n)
      end select     
      k2=hslo(n,m,0);k4=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k1,k2,k3,k4,x0,x1,xp(i,n))
      end do
      end do
   end do
   !--INTERFACE 2 & 3
      n=2
      do m = 0,2,2
      selectcase(m)
      case(0); is=lxis0;ie=lxie0;x0=xa(n);x1=xc(n)
      case(2); is=lxis2;ie=lxie2;x0=xd(n);x1=xf(n)
      end select     
      k1=yd(m);k3=yd(m+1);k2=hslo(n,m,0);k4=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k1,k2,k3,k4,x0,x1,xp(i,n))
         yp(i,n+1)=yp(i,n)
      end do
      end do

!--VERICAL END BOUNDARIES
   !-Y-DIRECTION
   do n = 0, 3
         alph=delt1
      !-BLOCK1
      !-c-d
      ip=lets1; im=let1;
      tmpa=yc(n);sha=smod*shs2;tmpb=yd(n);shb=shs2*sin(pi4-alph)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      she2=qet(ip)
      !-BLOCK0
      !-a-b
      sho=tla/litr; ll=lsz1/sho
      ip=lets0; im=ll;
      tmpa=ya(n);sha=sho;tmpb=yb(n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-b-c
      ip=ip+im; im=let0-im;
      tmpa=yb(n);sha=sho;tmpb=yc(n);shb=she2
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-BLOCK2
      !-d-e
      ip=lets2; im=let1;
      tmpa=yd(n);sha=shs2*sin(pi4+alph);tmpb=ye(n);shb=smod*shs2
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      she2=qet(ip+im)
      !-BLOCK3
      !-e-f
      sho=tla/litr; ll=lsz1/sho
      ip=lets3; im=let0-ll;
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
      case(0); is=lets0;ie=lete0;x0=ya(n);x1=yc(n)
               k1=xc(m);k3=xc(m+1)
      case(1); is=lets1;ie=lete1;x0=yc(n);x1=yd(n)
               k1=xc(m);k3=xc(m+1)
      case(2); is=lets2;ie=lete2;x0=yd(n);x1=ye(n)
               k1=xc(m+1);k3=xc(m+2)
      case(3); is=lets3;ie=lete3;x0=ye(n);x1=yf(n)
               k1=xc(m+1);k3=xc(m+2)
      end select     
      k2=vslo(n,m,0);k4=vslo(n,m,1)
      do i = is, ie
         xq(i,n)=inter(k1,k2,k3,k4,x0,x1,yq(i,n))
      end do
   end do
   !-INTERFACE 2
   xq(:,2)=xd(0)
   !-INTERFACE 3
   xq(:,3)=xf(0)

!--GRID INTERPOLATION
   do l=0,3
      select case(l);
      case(0); js=0; je=lete0; n=0; gf=0
      case(1); js=lets1; je=lete1; n=1; gf=1
      case(2); js=lets2; je=lete2; n=3; gf=0
      case(3); js=lets3; je=lete3; n=4; gf=1
      end select
     do m=0,2
      select case(m)
      case(0); is=0;     ie=lxie0; nn=1; rs=abs(xa(n+1)-xa(n)); re=abs(xc(n+1)-xc(n))
               tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
      case(1); is=lxis1; ie=lxie1; nn=0; rs=abs(xc(n+1)-xc(n)); re=abs(xd(n+1)-xd(n))
               tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
      case(2); is=lxis2; ie=lxie2; nn=1; rs=abs(xd(n+1)-xd(n)); re=abs(xe(n+1)-xe(n))
               tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
      end select
      do j=js,je; do i=is,ie
         pp=real(i-is,nr)/(ie-is); qq=real(j-js,nr)/(je-js)
         tmp=sin(half*pi*pp); ra0=(gf*real(je-j,nr)/(je-js)+gf*qq); ra1=one+(2-one)*ra0**2
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
     case(0,1,2); js=lets0; je=lete0; 
     case(3,4,5); js=lets1; je=lete1;
     case(6,7,8); js=lets2; je=lete2;
     case(9,10,11); js=lets3; je=lete3;
   end select
   select case(mb)
     case(0,3,6,9); is=lxis0; ie=lxie0;
     case(1,4,7,10); is=lxis1; ie=lxie1;
     case(2,5,8,11); is=lxis2; ie=lxie2
   end select
   np=nr*(je-js+1)*(ie-is+1)
   write(1,pos=np*k+1) ((xx(i,j),i=is,ie),j=js,je)
   write(1,pos=np*(k+lze0+1)+1) ((yy(i,j),i=is,ie),j=js,je)
   write(1,pos=np*(k+2*lze0+2)+1) ((zz(i,j),i=is,ie),j=js,je)
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

 integer,intent(in) :: i,n,m
 real(nr) :: yi

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
    real(nr),intent(in) :: x,c,t0
    integer, intent(in) :: n
    real(nr) :: y
    real(nr) :: t,k,xc

    t=t0*0.01e0
    k=0.991148635e0
    xc=x/c
    y =(t*c*k/0.2e0) &
    *(0.298222773e0*sqrt(xc)-0.127125232e0*xc- 0.357907906e0*xc**2+&
    0.291984971e0*xc**3-0.105174606e0*xc**4)
    if (n==2) then
       y=y*-1.0e0
    end if
 end function naca
 
!===== FUNCTION FOR INTERFACE POINTS
 function inter(k1,k2,k3,k4,x0,x1,x) result(y)

 real(nr) :: y
 real(nr), intent(in) :: x,x0,x1,k1,k2,k3,k4
 real(nr) :: a,b,c,d,l,invl,fx

 l=x1-x0
 invl=1.0_nr/l
 fx=(x-x0)*invl

    a=(2*fx**3-3*fx**2+1);
    b=(fx**3-2*fx**2+fx);
    c=(-2*fx**3+3*fx**2);
    d=(fx**3-fx**2);
    y=a*k1+b*l*k2+c*k3+d*l*k4;
 
 end function inter
 !=====

end module gridgen

!*****
