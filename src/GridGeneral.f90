!*****  
!***** 3D Generic-Grid
!*****

module gridgen

 use subroutineso
 implicit none

 real(k8),parameter :: pi4=pi/4.0_k8
 integer(k4) :: lxit,lett,lxisz,im,jm
 integer(k4) :: lxis,lxie,lxib
 integer(k4) :: lets,lete,letb

 real(k8),dimension(:,:,:),allocatable :: hslo,vslo

 real(k8),dimension(:,:),allocatable :: xx,yy,zz
 real(k8),dimension(:),allocatable :: zs
 real(k8),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(k8),dimension(:),allocatable :: pxi,qet

 real(k8) :: rs,re,rp,ts,te,shs1,she1,shs2,she2,shs,she,shswle
 real(k8) :: xo,xjct,yo,yjct,sho,pp,qq
 real(k8) :: am,err,tmp,tmpa,tmpb,gf
 real(k8) :: k01,k02,k03,k04,x0,x1
 real(k8) :: deg1,deg2

 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,litr,smgrid,&
            domlen,span,wlew,wlea,szth1,szth2,szxt,&
            tla,tlb,cutlb,c1,delt1,ximod,etamod)

 integer(k4),intent(in) :: ngridv,litr,nthick
 integer(k4), dimension(:,:) :: lxise(0:bkx-1,0:1)
 integer(k4), dimension(:,:) :: letse(0:bky-1,0:1)
 integer(k4), dimension(:,:) :: lzese(0:bkz-1,0:1)
 real(k8), dimension(:,:), allocatable :: px,py
 real(k8),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,c1,delt1
 real(k8),intent(in) :: ximod,etamod,cutlb
 real(k8) :: lsz1,lsz2
 real(k8) :: lbl,lwle
 real(k8) :: alph
 real(k8) :: oxp,oyp
 real(k8) :: tmps,tmpe,tmpc
 real(k8) :: sha,shb,shc
 integer(k4) :: smod,npx,npy
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
    lbl=0.2_k8

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lzebk(0)))
    allocate(pxi(0:lxit),qet(0:lett))
    npy=min(nthick,1)+bky;npx=bkx
    allocate(xp(0:lxit,0:npy),yp(0:lxit,0:npy))
    allocate(xq(0:lett,0:npx),yq(0:lett,0:npx))
    allocate(px(0:bkx+2,0:npy),py(0:bky+2,0:npx))
    allocate(hslo(0:npy,0:bkx-1,0:1),vslo(0:npx,0:bky-1,0:1))

if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream',form='unformatted')

 do k=0,lzebk(0)
    zs(k)=span*(real(lzebk(0)-k,k8)/lzebk(0)-half)
!---BLOCKS' BOUNDARIES
    lsz1=litr*shs1; lsz2=lsz1+szxt
    lbl=lsz1
!---VERTICAL LINES
    px(0,:)=-domlen
    px(1,:)=px(0,:)+lsz1
    px(2,:)=zero
    px(3,:)=domlen-lsz1
    px(4,:)=domlen
!---HORIZONTAL LINES
    py(0,:)=-domlen
    py(1,:)=py(0,:)+lbl
    py(2,:)=zero
    py(3,:)=domlen-lbl
    py(4,:)=domlen

!----- INITIAL AND END HORIZONTAL SLOPES
    hslo(:,:,:)=zero

!----- INITIAL AND END VERTICAL SLOPES
    vslo(:,:,:)=zero

!--HORIZONTAL INTERFACES
   do n = 0,bky
   !--X-COORDINATE
   !--BLOCK0
   !--0-1
      ll=litr
      ip=lxise(0,0); im=ll;
      tmpa=px(0,n);sha=shs1;tmpb=px(1,n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--1-2
      ip=ip+im; im=lxibk(0)-ll;
      tmpa=px(1,n);sha=shb;tmpb=px(2,n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
      if(k==0.and.n==0) then
         lxisz=lxibk(2)*(minloc(abs(px(0,n)+szth1-xp(0:lxibk(0),n)),1)-1)/lxibk(0); lp=ll     
      end if
   !--BLOCK1
   !--2-3
      ip=lxise(1,0); im=lxibk(1)-ll;
      tmpa=px(2,n);sha=shb;tmpb=px(3,n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--3-4
      ip=ip+im; im=ll;
      tmpa=px(3,n);sha=shb;tmpb=px(4,n);shb=sha
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   end do

   !--Y-COORDINATE
   do n = 0,bky
      do m = 0,bkx-1
      selectcase(n)
       case(1);     k01=py(2,m);k03=py(2,m+1)
       case(2);     k01=py(4,m);k03=py(4,m+1)
       case default;k01=py(n,m);k03=py(n,m+1) 
      end select
      selectcase(m)
       case(0);      is=lxise(0,0);ie=lxise(0,1);x0=px(0,n);x1=px(2,n)
       case(1);      is=lxise(1,0);ie=lxise(1,1);x0=px(2,n);x1=px(4,n)    
      end select     
      k02=hslo(n,m,0);k04=hslo(n,m,1)
      do i = is, ie
         yp(i,n)=inter(k01,k02,k03,k04,x0,x1,xp(i,n))
      end do
      end do
   end do

!--VERICAL END BOUNDARIES
   !-Y-COORDINATE
   do n = 0, bkx
      !-BLOCK0
      !-0-1
      ll=litr
      ip=letse(0,0); im=ll;
      tmpa=py(0,n);sha=shs1;tmpb=py(1,n);shb=shs1
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-1-2
      ip=ip+im; im=letbk(0)-ll;
      tmpa=py(1,n);sha=shb;tmpb=py(2,n);shb=shs1
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-BLOCK1
      !-2-3
      ip=letse(1,0); im=letbk(1)-ll;
      tmpa=py(2,n);sha=shs1;tmpb=py(3,n);shb=shs1
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)   
      !-3-4
      ip=ip+im; im=ll
      tmpa=py(3,n);sha=shb;tmpb=py(4,n);shb=shs1
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
   end do

   !-X-COORDINATE
   do n = 0,bkx
      do m = 0,bky-1
      selectcase(n)
       case(0);k01=px(0,m);k03=px(0,m+1) 
       case(1);k01=px(2,m);k03=px(2,m+1)
       case(2);k01=px(4,m);k03=px(4,m+1)
      end select
      selectcase(m)
       case(0);      is=letse(0,0);ie=letse(0,1);x0=py(0,n);x1=py(2,n)
       case(1);      is=letse(1,0);ie=letse(1,1);x0=py(2,n);x1=py(4,n)    
      end select     
      k02=vslo(n,m,0);k04=vslo(n,m,1)
      do i = is, ie
         xq(i,n)=inter(k01,k02,k03,k04,x0,x1,yq(i,n))
      end do
      end do
   end do

!--GRID INTERPOLATION
   do l=0,bky-1
     js=letse(l,0); je=letse(l,1); n=l;
     do m=0,bkx-1
      is=lxise(m,0); ie=lxise(m,1)
      call hermite2D(is,ie,js,je,n,m,zs(k))
     end do
   end do

!--GRID OUTPUT
   ll=mb/bkx
   js=letse(ll,0);je=letse(ll,1)
   ll=mod(mb,bkx)
   is=lxise(ll,0);ie=lxise(ll,1)
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

 subroutine hermite2D(is,ie,js,je,n,m,z)
 implicit none
 integer(k4), intent(in) :: is,ie,js,je,n,m
 real(k8), intent(in) :: z
 integer(k4) :: lx,side,lxi,let
 integer(k4) :: i,j,ii,jj
 real(k8), dimension(:,:), allocatable :: xv,xu,xuv,yv,yu,yuv
 real(k8), dimension(4,4) :: mx,my
 real(k8), dimension(1,4) :: hu,ux,uy
 real(k8), dimension(4,1) :: hv,vx,vy
 real(k8), dimension(1,1) :: val1,val2,val3
 real(k8) :: u
 
 lxi=ie-is; let=je-js
 
 xx(is:ie,js)=xp(is:ie,n); xx(is:ie,je)=xp(is:ie,n+1)
 yy(is:ie,js)=yp(is:ie,n); yy(is:ie,je)=yp(is:ie,n+1)
 xx(is,js:je)=xq(js:je,m); xx(ie,js:je)=xq(js:je,m+1)
 yy(is,js:je)=yq(js:je,m); yy(ie,js:je)=yq(js:je,m+1)
 
 allocate(xv(0:lxi,0:1),yv(0:lxi,0:1),&
          xu(0:let,0:1),yu(0:let,0:1),&
          xuv(0:1,0:1),yuv(0:1,0:1))
 
 lx=let; side=0
 xv(0  ,side)=drvbn(xx(is,js:je),lx,side)
 xv(lxi,side)=drvbn(xx(ie,js:je),lx,side)
 side=1
 xv(0  ,side)=drvbn(xx(is,js:je),lx,side)
 xv(lxi,side)=drvbn(xx(ie,js:je),lx,side)
 
 lx=let; side=0
 yv(0  ,side)=drvbn(yy(is,js:je),lx,side)
 yv(lxi,side)=drvbn(yy(ie,js:je),lx,side)
 side=1
 yv(0  ,side)=drvbn(yy(is,js:je),lx,side)
 yv(lxi,side)=drvbn(yy(ie,js:je),lx,side)
 
 lx=lxi; side=0
 xu(0  ,side)=drvbn(xx(is:ie,js),lx,side)
 xu(let,side)=drvbn(xx(is:ie,je),lx,side)
 side=1
 xu(0  ,side)=drvbn(xx(is:ie,js),lx,side)
 xu(let,side)=drvbn(xx(is:ie,je),lx,side)
 
 lx=lxi; side=0
 yu(0  ,side)=drvbn(yy(is:ie,js),lx,side)
 yu(let,side)=drvbn(yy(is:ie,je),lx,side)
 side=1
 yu(0  ,side)=drvbn(yy(is:ie,js),lx,side)
 yu(let,side)=drvbn(yy(is:ie,je),lx,side)
 
 
 lx=lxi
 do i = 0, lx
    xv(i,0) =hermite(xv(0,0),zero,xv(lxi,0),zero,lx,i)
    xv(i,1) =hermite(xv(0,1),zero,xv(lxi,1),zero,lx,i)
    yv(i,0) =hermite(xv(0,0),zero,xv(lxi,0),zero,lx,i)
    yv(i,1) =hermite(xv(0,1),zero,xv(lxi,1),zero,lx,i)
 end do
 lx=let
 do i = 0, lx
    xu(i,0) =hermite(xu(0,0),zero,xu(lxi,0),zero,lx,i)
    xu(i,1) =hermite(xu(0,1),zero,xu(lxi,1),zero,lx,i)
    yu(i,0) =hermite(xu(0,0),zero,xu(lxi,0),zero,lx,i)
    yu(i,1) =hermite(xu(0,1),zero,xu(lxi,1),zero,lx,i)
 end do
 mx(1,:)=(/xx(is,js),xx(is,je),xv ( 0 ,0),xv( 0 ,1)/)
 mx(2,:)=(/xx(ie,js),xx(ie,je),xv (lxi,0),xv(lxi,1)/)
 mx(3,:)=(/xu(0 ,0 ),xu(let,0),        0.0_k8 ,0.0_k8/)
 mx(4,:)=(/xu(0 ,1 ),xu(let,1),        0.0_k8 ,0.0_k8/)
 
 my(1,:)=(/yy(is,js),yy(is,je),yv ( 0 ,0),yv( 0 ,1)/)
 my(2,:)=(/yy(ie,js),yy(ie,je),yv (lxi,0),yv(lxi,1)/)
 my(3,:)=(/yu(0 ,0 ),yu(let,0),        0.0_k8 ,0.0_k8/)
 my(4,:)=(/yu(0 ,1 ),yu(let,1),        0.0_k8 ,0.0_k8/)
 
 do j = js, je; u=dble(j-js)/let;jj=j-js
    hv(:,1)=fillh(u)
    vx(:,1)=(/xx(is,j),xx(ie,j),xu(jj,0),xu(jj,1)/)
    vy(:,1)=(/yy(is,j),yy(ie,j),yu(jj,0),yu(jj,1)/)
    do i = is, ie; u=dble(i-is)/lxi;ii=i-is
       hu(1,:)=fillh(u)
       ux(1,:)=(/xx(i,js),xx(i,je),xv(ii,0),xv(ii,1)/)
       uy(1,:)=(/yy(i,js),yy(i,je),yv(ii,0),yv(ii,1)/)
       val1=matmul(hu,vx);val2=matmul(ux,hv);val3=matmul(hu,matmul(mx,hv))
       xx(i,j)=val1(1,1)+val2(1,1)-val3(1,1)
       val1=matmul(hu,vy);val2=matmul(uy,hv);val3=matmul(hu,matmul(my,hv))
       yy(i,j)=val1(1,1)+val2(1,1)-val3(1,1)
    end do
 end do
 
 zz(is:ie,js:je)=z
    
 end subroutine hermite2D


!===== FUNCTION FOR HERMITIAN INTERPOLATION
 function hermite(k01,k02,k03,k04,lx,x) result(y)
 real(k8) :: y
 integer(k4), intent(in) :: lx,x
 real(k8), intent(in) :: k01,k02,k03,k04  
 real(k8) :: a,b,c,d,l,invl,fx
 real(k8), dimension(0:3) :: v

 l=dble(lx)
 invl=1.0_k8/l
 fx=x*invl

 v=fillh(fx)

 y=v(0)*k01+v(2)*l*k02+v(1)*k03+v(3)*l*k04;
 
 end function hermite
 function fillh(u) result(v)
    real(k8), dimension(0:3) :: v
    real(k8) :: a0,a1,b0,b1
    real(k8), intent(in) :: u
 
    a0=2*u**3-3*u**2+1;
    a1=3*u**2-2*u**3
    b0=u**3-2*u**2+u  ;
    b1=u**3-u**2
 
    v(:)=(/a0,a1,b0,b1/)
 end function fillh

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

!===== FUNCTION FOR BOUNDARY DERIVATIVE
 function drvbn(x,lx,side) result (dx)
 implicit none
 real(k8) :: dx
 integer(k4), intent(in) :: lx,side
 real(k8), dimension(0:lx), intent(in) :: x

 select case(side)
 case(0); dx=(-1.5_k8)*x(0)+2.0_k8*x(1)-0.5_k8*x(2);
 case(1); dx=1.5_k8*x(lx)-2.0_k8*x(lx-1)+0.5_k8*x(lx-2)
 end select
 end function drvbn

!===== FUNCTION DEGREES TO RADIANS
 function tand(ad) result (r)
 implicit none
 real(k8) :: r
 real(k8), intent(in) :: ad
 real(k8) :: ar
 ar=ad*pi/180.0_k8
 r=tan(ar)
 end function tand
 !=====

end module gridgen

!*****               

