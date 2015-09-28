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
    lbl=cutlb*half/sin(pi4)
    !lbl=min(lbl,cutlb*half*0.98/sin(pi4))
    if (myid==0) then
       write(*,*) lbl
    end if

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lzebk(0)))
    allocate(pxi(0:lxit),qet(0:lett))
    npy=min(nthick,1)+bky;npx=bkx
    allocate(xp(0:lxit,0:npy),yp(0:lxit,0:npy))
    allocate(xq(0:lett,0:npx),yq(0:lett,0:npx))
    allocate(px(0:bkx+2,0:npy),py(0:bky+2,0:npx))

if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream',form='unformatted')

 do k=0,lzebk(0)
    zs(k)=span*(real(lzebk(0)-k,k8)/lzebk(0)-half)
!---BLOCKS' BOUNDARIES
    sho=tla/litr; ll=2*litr; lsz1=ll*sho; lsz2=szth1+szxt
!---WAVY LEADING-EDGE PROFILE
    lwle=wlea*sin(2*pi*(zs(k)-zs(0))/wlew)
!---VERTICAL LINES
    
    px(0,:)=-domlen
    px(1,:)=px(0,:)+lsz1
    px(2,2:3)=-half*c1+lwle*cos(delt1)
    px(2,0:1)=px(2,2)-lbl*cos(pi4-delt1); px(2,4:5)=px(2,3)-lbl*cos(pi4-delt1)
    if (nthick*(nthick-3)==0) then
       px(2,:)=-half*c1+lwle*cos(delt1)
    end if
    px(3,:)=c1*(cos(delt1)-half)
    px(bkx+1,:)=domlen-szth1
    px(bkx+2,:)=px(bkx+1,:)+lsz2
!---HORIZONTAL LINES
    py(0,:)=-domlen
    py(1,:)=py(0,:)+lsz1
    py(3,:)=zero
    py(3,1)=-lwle*sin(delt1); py(3,2)=-c1*sin(delt1)
    py(4,0)=py(3,0)+lbl;py(4,1)=py(3,1)+lbl*sin(pi4+delt1)
    py(4,2)=py(3,2)+lbl*sin(pi4)*cos(delt1)
    py(4,3)=py(4,0)
    py(2,0)=py(4,0)-2*lbl
    py(2,1)=py(3,1)-lbl*sin(pi4-delt1)
    py(2,2)=py(3,2)-lbl*sin(pi4)*cos(delt1)
    py(2,3)=py(2,0)
    py(bky+2,:)=domlen
    py(bky+1,:)=py(bky+2,:)-lsz1

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
    tmp=c1-lwle;tmpa=px(2,2);tmpb=py(3,1);lxis=lxise(1,0);lxie=lxise(1,1);lxib=lxibk(1);alph=delt1
    if (nthick*(nthick-3)==0) then
       flag=.false.
    else
       flag=.true.
    end if
    if (flag) then
       ! DETERMINE UPPER AND LOWER SIDES
       do n=2,3
          yp(lxis,n)=zero;xp(lxis,n)=zero
          ! DETERMINE THE FIRST LL POINTS
          ll=8 ! "LL" MUST BE EQUAL TO OR LARGER THAN 4.
          do i=lxis+1,lxis+ll
             xp(i,n)=xp(i-1,n)+half*shs1; err=1
             do while(abs(err)>sml)
                yp(i,n)=naca(xp(i,n),tmp,21.0_k8,n)
                err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shs1-1;
                xp(i,n)=xp(i,n)-half**5*err*shs1
             end do
          end do
          xo=xp(lxis+ll,n); sho=sum(xp(lxis+ll-4:lxis+ll,n)*(/3,-16,36,-48,25/))/12
          ! COMPUTE THE REST OF THE POINTS
          ip=lxis+ll;im=lxib-ll 
          call gridf(xp(:,n),pxi,xo,tmp,sho,she1,lxit,im,ip)
          do i=lxis+ll+1,lxie-1
             yp(i,n)=naca(xp(i,n),tmp,21.0_k8,n)
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
      tmpa=px(0,n);sha=sho;tmpb=px(1,n);shb=sho
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--b-c
      ip=ip+im; im=lxibk(0)-ll;
      tmpa=px(1,n);sha=sho;tmpb=px(2,n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
      if(k==0.and.n==0) then
         lxisz=lxibk(2)*(minloc(abs(px(0,n)+szth1-xp(0:lxibk(0),n)),1)-1)/lxibk(0); lp=ll     
      end if
   if (n.ne.2) then
   !--BLOCK1
   !--c-d
      ip=lxise(1,0); im=lxibk(1);
      tmpa=px(2,n);sha=shs1;tmpb=px(3,n);shb=she1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   end if
   !--BLOCK2
   !--d-e
      ip=lxise(2,0); im=lxibk(2)-lxisz;
      tmpa=px(3,n);sha=she1;tmpb=px(4,n);shb=sml
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--e-f
      ip=ip+im; im=lxisz;
      tmpa=px(4,n);sha=pxi(ip);tmpb=px(5,n);shb=sho
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
      case(0); k01=py(0,m);k03=py(0,m+1)
      case(1); k01=py(2,m);k03=py(2,m+1)
      end select
      selectcase(m)
      case(0); is=lxise(0,0);ie=lxise(0,1);x0=px(0,n);x1=px(2,n)
      case(1); is=lxise(1,0);ie=lxise(1,1);x0=px(2,n);x1=px(3,n)    
            case(2); is=lxise(2,0);ie=lxise(2,1);x0=px(3,n);x1=px(5,n)
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
      case(4); k01=py(4,m);k03=py(4,m+1)
      case(5); k01=py(bky+2,m);k03=py(bky+2,m+1)
      end select
      selectcase(m)
      case(0); is=lxise(0,0);ie=lxise(0,1);x0=px(0,n);x1=px(2,n)
      case(1); is=lxise(1,0);ie=lxise(1,1);x0=px(2,n);x1=px(3,n)
      case(2); is=lxise(2,0);ie=lxise(2,1);x0=px(3,n);x1=px(5,n)
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
      case(0); is=lxise(0,0);ie=lxise(0,1);x0=px(0,n);x1=px(2,n)
      case(2); is=lxise(2,0);ie=lxise(2,1);x0=px(3,n);x1=px(5,n)
      end select     
      k01=py(3,m);k03=py(3,m+1);k02=hslo(n,m,0);k04=hslo(n,m,1)
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
      tmpa=py(2,n);sha=smod*shs2;tmpb=py(3,n);shb=shs2*sin(pi4-alph)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      she2=qet(ip)
      !-BLOCK0
      !-a-b
      sho=tla/litr; ll=lsz1/sho
      ip=letse(0,0); im=ll;
      tmpa=py(0,n);sha=sho;tmpb=py(1,n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)   
            !-b-c
      ip=ip+im; im=letbk(0)-im;
      tmpa=py(1,n);sha=sho;tmpb=py(2,n);shb=she2
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-BLOCK2
      !-d-e
      ip=letse(2,0); im=letbk(2);
      tmpa=py(3,n);sha=shs2*sin(pi4+alph);tmpb=py(4,n);shb=smod*shs2
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      she2=qet(ip+im)
      !-BLOCK3
      !-e-f
      sho=tla/litr; ll=lsz1/sho
      ip=letse(3,0); im=letbk(3)-ll;
      tmpa=py(4,n);sha=she2;tmpb=py(bky+1,n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
      !-f-g
      ip=ip+im; im=ll;
      tmpa=py(bky+1,n);sha=sho;tmpb=py(bky+2,n);shb=sho
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lxit,im,ip)
   end do
   !-COPY ON N+5
   !yq(:,3)=yq(:,0)

   !-X-DIRECTION
   !-INTERFACE 0
   xq(:,0)=px(0,0)
   !-INTERFACE 1
   n=1
   do m = 0,3
      selectcase(m)
      case(0); is=letse(0,0);ie=letse(0,1);x0=py(0,n);x1=py(2,n)
               k01=px(2,m);k03=px(2,m+1)
      case(1); is=letse(1,0);ie=letse(1,1);x0=py(2,n);x1=py(3,n)
               k01=px(2,m);k03=px(2,m+1)
      case(2); is=letse(2,0);ie=letse(2,1);x0=py(3,n);x1=py(4,n)
               k01=px(2,m+1);k03=px(2,m+2)
      case(3); is=letse(3,0);ie=letse(3,1);x0=py(4,n);x1=py(bky+1,n)
               k01=px(2,m+1);k03=px(2,m+2)
      end select     
      k02=vslo(n,m,0);k04=vslo(n,m,1)
      do i = is, ie
         xq(i,n)=inter(k01,k02,k03,k04,x0,x1,yq(i,n))
      end do
   end do
   !-INTERFACE 2
   xq(:,2)=px(3,0)
   !-INTERFACE 3
   xq(:,3)=px(5,0)

!--GRID INTERPOLATION
   do l=0,3
      select case(l);
      case(0); js=letse(0,0); je=letse(0,1); n=0; gf=0
      case(1); js=letse(1,0); je=letse(1,1); n=1; gf=0
      case(2); js=letse(2,0); je=letse(2,1); n=3; gf=1       
      case(3); js=letse(3,0); je=letse(3,1); n=4; gf=1
      end select
     do m=0,2
      select case(m)
      case(0); is=lxise(0,0); ie=lxise(0,1); nn=0;
               rs=abs(px(0,n+1)-px(0,n)); re=abs(px(2,n+1)-px(2,n))
      case(1); is=lxise(1,0); ie=lxise(1,1); nn=0;
               rs=abs(px(2,n+1)-px(2,n)); re=abs(px(3,n+1)-px(3,n))
      case(2); is=lxise(2,0); ie=lxise(2,1); nn=0;
               rs=abs(px(3,n+1)-px(3,n)); re=abs(px(5,n+1)-px(5,n))
      end select
               tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
      call hermite2D(is,ie,js,je,n,m,zs(k))
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

subroutine linear2D(is,ie,js,je,n,m,z)
implicit none
integer(k4), intent(in) :: is,ie,js,je,n,m
real(k8), intent(in) :: z
real(k8) :: t,s

do j = js, je
   t=dble(j-js)/dble(je-js)
   do i = is, ie
      s=dble(i-is)/dble(ie-is)
      xx(i,j)=(1.0-s)*xq(j,m)+s*xq(j,m+1)&
             +(1.0-t)*xp(i,n)+t*xp(i,n+1)&
             -(1.0-s)*(1.0-t)*xp(js,n)&
             -(1.0-s)*t*xp(js,n+1)&
             -s*(1.0-t)*xq(je,n)&
             -s*t*xq(je,n+1)
      yy(i,j)=(1.0-s)*yq(j,m)+s*yq(j,m+1)&
             +(1.0-t)*yp(i,n)+t*yp(i,n+1)&
             -(1.0-s)*(1.0-t)*yp(js,n)&
             -(1.0-s)*t*yp(js,n+1)&
             -s*(1.0-t)*yq(je,n)&
             -s*t*yq(je,n+1)
      zz(i,j)=z
   end do
end do
end subroutine linear2D

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

 l=dble(lx)
 invl=1.0_k8/l
 fx=x*invl

    a=(2*fx**3-3*fx**2+1);
    b=(fx**3-2*fx**2+fx);
    c=(-2*fx**3+3*fx**2);
    d=(fx**3-fx**2);
    y=a*k01+b*l*k02+c*k03+d*l*k04;
 
 end function hermite
function fillh(u) result(v)
   real(k8), dimension(0:3) :: v
   real(k8) :: a0,a1,b0,b1
   real(k8), intent(in) :: u

   a0=2*u**3-3*u**2+1; a1=3*u**2-2*u**3
   b0=u**3-2*u**2+u  ; b1=u**3-u**2

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
 !=====

end module gridgen

!*****               

