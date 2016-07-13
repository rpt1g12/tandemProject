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
 real(k8),dimension(:,:),allocatable :: pxi,qet

 real(k8) :: rs,re,rp,ts,te,shs1,she1,shs2,she2,shs,she,shswle
 real(k8) :: xo,xjct,yo,yjct,sho,pp,qq
 real(k8) :: am,err,tmp,tmpa,tmpb,gf
 real(k8) :: k01,k02,k03,k04,x0,x1
 real(k8) :: deg1,deg2

 real(k8),dimension(0:1,0:1) :: szth,dlth,szsh
 real(k8), dimension(0:1) :: lhbl,lvbl
 integer(k4),dimension(0:1,0:1) :: szll,wkll
 integer(k4),dimension(0:1) :: nvbl
 integer(k4),dimension(0:1) :: nwk,nwk2
 real(k8),dimension(0:1) :: lwk,lwk2

 contains

!===== GRID GENERATION
!        <==Grid Sketch==>
!      ^
!    py+>0             1         2                        3
!     3  |----|--------|---------|--------------|---------| 5
!        |             |         |                        |
!        |- - |- - - - |- - - - -|- - - - - - - |- - - - -| 4
!        |             |         |                        |
!        |    |        |         |              |         |
!        |             |         |                        |
!        |    |        |         |              |         |
!        |             |         |                        |
!        |    |        |         |              |         |
!    1-2 |=============|<=======>|========================| 2-3
!        |    |        |         |              |         |
!        |             |         |                        |
!        |    |        |         |              |         |
!        |             |         |                        |
!        |    |        |         |              |         |
!        |             |         |                        |
!        |- - |- - - - |- - - - -|- - - - - - - |- - - - -| 4
!        |             |         |                        |
!  ^     |    |        |         |              |         |
!px+> 0  |-------------|---------|------------------------| 5
!        0    1        2         3              4         5
 subroutine gridaerofoil(ngridv,nthick,smgrid,&
            domlen,span,wlew,wlea,szth1,szth2,szxt,&
            c1,delt1,ximod,etamod)

 integer(k4),intent(in) :: ngridv,nthick
 integer(k4), dimension(:,:) :: lxise(0:bkx-1,0:1)
 integer(k4), dimension(:,:) :: letse(0:bky-1,0:1)
 integer(k4), dimension(:,:) :: lzese(0:bkz-1,0:1)
 real(k8), dimension(:,:), allocatable :: px,py
 real(k8),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,c1,delt1
 real(k8),intent(in) :: ximod,etamod
 real(k8) :: lsz1,lsz2
 real(k8) :: lwle
 real(k8) :: alph,thk=21
 real(k8) :: oxp,oyp
 real(k8) :: tmps,tmpe,tmpc
 real(k8) :: sha,shb,shc
 integer(k4) :: npx,npy,mm,nn
 integer(k4), dimension(0:1) :: smod
 integer(k4), dimension(:), allocatable :: linesy,linesx
 real(k8), dimension(:), allocatable :: degarr
 logical :: flag

    ! rpt-Total number of points in xi and eta directions
    lxit=sum(lxibk(:))+(bkx-1); lett=sum(letbk(:))+(bky-1)
    ! rpt-Start and end indices in each direction
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

    ! rpt-Smallest grid sizes
    shs=smgrid; she=shs
    shs1=ximod*smgrid; ! rpt-LE xi size 
    shs2=etamod*smgrid;! rpt-LE eta size
    she1=shs2          ! rpt-TE size both xi and eta
    smod(:)=(/4,3/) ! grid size modifiers

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lzebk(0)))

    ! rpt-Number of horizontal lines
    npy=min(nthick,1)+bky;npx=bkx

    allocate(xp(0:lxit,0:npy),yp(0:lxit,0:npy))
    allocate(xq(0:lett,0:npx),yq(0:lett,0:npx))
    allocate(pxi(0:lxit,0:npy),qet(0:lett,0:npx))
    allocate(px(0:bkx+2,0:npy),py(0:npy+2,0:npx))
    allocate(hslo(0:npy,0:bkx-1,0:1),vslo(0:npx,0:bky-1,0:1))

! rpt-Assign grid file names
if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream',form='unformatted')

!---Domain Sizes
    dlth(0,0)=0.7e0*domlen; dlth(0,1)=domlen+szxt
    dlth(1,0)=0.7e0*domlen; dlth(1,1)=domlen+szxt
!---Sponge thicknesses
    szth(0,0)=szth1; szth(0,1)=szth2+szxt ! rpt-Horizontal direction left/right boudaries
    szth(1,0)=szth1; szth(1,1)=szth2+szxt ! rpt-Vertical direction bottom/top boudaries
!----- CONSTANT ANGLES IN RADIANS
    if(.not.allocated(degarr)) allocate(degarr(3))
    degarr(:)=(/25_k8,0_k8,15_k8/); degarr(:)=degarr(:)*pi/180_k8
!---Wake Refinement
    nwk(0)=int(lxibk(2)*0.45e0) ! rpt-Wake box #xi points
    nwk(1)=int(letbk(1)*0.4e0) ! rpt-Wake box #eta points
    nwk2(0)=0.9e0*nwk(0) !rpt-wake refinement in outflow #xi points
    nwk2(1)=1.0e0*nwk(1) !rpt-wake refinement in outflow #eta points
    lwk(1)=0.5e0*c1 ! rpt-Wake box size eta direction
    lwk(0)=1.5e0*c1!real(nwk(0)/nwk(1))*lwk(1) ! rpt-Wake box size xi direction
    lwk2(0)=1.5*lwk(0) !rpt-wake refinement in outflow xi length
    lwk2(1)=1.8e0*lwk(1) !rpt-wake refinement in outflow eta length
    if (myid==0) then
       write(*,"('Wake box size: xi=',f8.4,' eta=',f8.4)")&
       lwk(0),lwk(1)
       write(*,"('Wake box Pts: xi=',i4,' eta=',i4)")&
       nwk(0),nwk(1)
       write(*,"('Wake outflow size: xi=',f8.4,' eta=',f8.4)")&
       lwk2(0),lwk2(1)
       write(*,"('Wake outflow Pts: xi=',i4,' eta=',i4)")&
       nwk2(0),nwk2(1)
    end if
!---Boundary Layer Refinement
    lhbl(0)=0.15*c1 ! rpt-LE curve bottom-horizontal lenght
    lhbl(1)=0.15*c1 ! rpt-LE curve top-horizontal lenght
    lvbl(0)=half*lwk(1) ! rpt-LE curve bottom-vertical lenght
    lvbl(1)=lwk(1) ! rpt-LE curve top-vertical lenght
    nvbl(0)=letbk(0)*0.3e0 ! rpt-#eta points bottom LE curve
    nvbl(1)=nwk(1) ! rpt-#eta points top LE curve
    if (myid==0) then
       write(*,"('BL curve horizontal size: south=',f8.4,' north=',f8.4)")&
       lhbl(0),lhbl(1)
       write(*,"('BL curve vertical size: south=',f8.4,' north=',f8.4)")&
       lvbl(0),lvbl(1)
       write(*,"('BL Pts: south=',i4,' north=',i4)")&
       nvbl(0),nvbl(1)
    end if
 do k=0,lzebk(0)
!---SPANWISE COORDINATE (UNIFORM)
    zs(k)=span*(real(lzebk(0)-k,k8)/lzebk(0)-half)
!---WAVY LEADING-EDGE PROFILE
    lwle=wlea*sin(2*pi*(zs(k)-zs(0))/wlew)
!---VERTICAL LINES
    px(0,:)=-dlth(0,0)
    px(1,:)=px(0,:)+szth(0,0)
    px(2,1:2)=-half*c1+lwle*cos(delt1); ! rpt-LE x location
    px(2,0)=px(2,2)-lhbl(0)
    px(2,3)=px(2,2)-lhbl(1)
    px(3,:)=-half*c1+c1*cos(delt1)
    px(3,3)=px(3,2)
    px(4,:)=px(3,0)+dlth(0,1)-szth(0,1)
    px(5,:)=px(3,0)+dlth(0,1)
!---HORIZONTAL LINES
    py(0,:)=-dlth(1,0)
    py(1,:)=py(0,:)+szth(1,0)
    py(2,0)=zero
    py(2,1)=-lwle*sin(delt1)
    py(2,2)=-c1*sin(delt1)
    py(2,3)=py(2,2)!-3.5e0*lhbl(1)
    py(3,:)=py(2,:)
    py(4,:)=py(2,2)+dlth(1,1)-szth(1,1)
    py(5,:)=py(2,2)+dlth(1,1)

!----- INITIAL AND END HORIZONTAL SLOPES
    !hslo(block,hline,start:end)
    hslo(0,:,:)=zero
    hslo(1,0,:)=(/zero,-tan(delt1)/)
    hslo(1,1,:)=zero
    hslo(1,2,:)=zero
    hslo(1,2,:)=(/-tan(delt1+2*degarr(2)),zero/)
    hslo(3,:,:)=zero

!----- INITIAL AND END VERTICAL SLOPES
    !hslo(block,vline,start:end)
    vslo(0,:,:)=zero
    vslo(1,0,:)=(/zero,tan(pi4+delt1)/)
    vslo(1,1,:)=(/-tan(pi4-delt1),zero/)
    vslo(2,0,:)=zero
    vslo(2,1,:)=zero
    vslo(3,:,:)=zero


!--HORIZONTAL INTERFACES
   !--X-COORDINATE
   do n = 0,npy
   !--BLOCK0
   !!--0-2 Left boundary->LE
   if (n.eq.0) then !(bottom boundary)
      ip=lxise(0,0); im=lxibk(0);
      tmpa=px(0,n);sha=220*shs1;tmpb=px(2,n);shb=shs1*smod(0);ra0=shb
      call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
         if (k==0.and.myid==0) then
            write(*,*) 'Block 0 hztnl: mesh size ratio',pxi(ip,n)/ra0
         end if
   elseif (n.eq.npy) then !(top boundary)
      ip=lxise(0,0); im=lxibk(0);
      tmpa=px(0,n);sha=220*shs1;tmpb=px(2,n);shb=shs1*smod(0)
      call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
   else !(horizontal interface)
      ip=lxise(0,0); im=lxibk(0);
      tmpa=px(0,n);sha=220*shs1;tmpb=px(2,n);shb=shs1;ra0=shb
      call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
   end if
   !--BLOCK1
   !--2-3 LE->TE (just top/bottom lines)
   select case(n)
   case(0)
      ip=lxise(1,0); im=lxibk(1);
      tmpa=px(2,n);sha=shs1*smod(0);tmpb=px(3,n);shb=she1*smod(0)
      call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
   case(3)
      ip=lxise(1,0); im=lxibk(1);
      tmpa=px(2,n);sha=shs1*smod(0);tmpb=px(3,n);shb=she1*smod(0)
      call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
   end select
   !--BLOCK2
   if (n.ne.npy) then !(all but the top boundary)
   !--3-(3+lwk(0)) Trailing edge wake box refinement
     ip=lxise(2,0); im=nwk(0)
     tmpa=px(3,n);sha=she1;tmpb=px(3,n)+lwk(0);shb=sml;ra0=sha
     call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
   !--(3+lwk(0))-5  Wake box->Right boundary
     ip=ip+im; im=lxibk(2)-im
     tmpa=tmpb;sha=pxi(ip,n);tmpb=px(5,n);shb=200*ra0!sml
     call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
         if (k==0.and.myid==0) then
            if(n==1) write(*,*) 'Horz Intfce: mesh size ratio',pxi(lxise(2,1),n)/ra0
         end if
     if (n==npy-1) ra1=shb
   else !(top boundary)
   !--3-(3+lwk(0)) Trailing edge wake box refinement
     ip=lxise(2,0); im=nwk2(0)
     tmpa=px(3,n);sha=she1*smod(0);ra0=sha
     tmpb=px(3,n)+lwk2(0);shb=sml
     call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
   !--(3+lwk(0))-5  Wake box->Right boundary
     ip=ip+im; im=lxibk(2)-im
     tmpa=tmpb;sha=pxi(ip,n);tmpb=px(5,n);shb=1.2e0*ra1!48*ra0!sml
     call gridf(xp(:,n),pxi(:,n),tmpa,tmpb,sha,shb,lxit,im,ip)
         if (k==0.and.myid==0) then
            write(*,*) 'Top boundary: mesh size ratio',pxi(lxise(2,1),n)/ra0
         end if
   end if
   end do

   !--Y-COORDINATE for horizontal lines
   !--Lower Boundarie
      yp(:,0)=py(0,0)
   !--Middle Interface 1
   !--BLOCK0
      k01=py(2,0);k03=py(2,1)
      k02=hslo(1,0,0);k04=hslo(1,0,1)
      x0=px(0,1);x1=px(2,1)
      is=lxise(0,0);ie=lxise(0,1);
      do i = is, ie
         yp(i,1)=inter(k01,k02,k03,k04,x0,x1,xp(i,1))
         yp(i,2)=yp(i,1) ! Copy interface
      end do
   !--BLOCK2
      k01=py(2,2);k03=py(2,3)
      k02=hslo(1,2,0);k04=hslo(1,2,1)
      x0=px(3,1);x1=px(5,1)
      is=lxise(2,0);ie=lxise(2,1);
      do i = is, ie
         yp(i,1)=inter2(k01,k02,k03,k04,x0,x1,xp(i,1),1)
         yp(i,2)=yp(i,1) ! Copy interface
      end do
   !--Upper Boundarie
      yp(:,3)=py(5,0)

!----- AEROFOIL SURFACE GRID POINTS
    tmp=c1-lwle;tmpa=px(2,1);tmpb=py(2,1);
    lxis=lxise(1,0);lxie=lxise(1,1);lxib=lxibk(1);
    alph=delt1;
    ! DETERMINE UPPER AND LOWER SIDES
    do n=1,2
       yp(lxis,n)=zero;xp(lxis,n)=zero
       ! DETERMINE THE FIRST LL POINTS
       ll=25 ! "LL" MUST BE EQUAL TO OR LARGER THAN 4.
       do i=lxis+1,lxis+ll
          xp(i,n)=xp(i-1,n)+half*shs1; err=1
          do while(abs(err)>sml)
             yp(i,n)=naca(xp(i,n),tmp,thk,n-1)
             err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shs1-1;
             xp(i,n)=xp(i,n)-half**8*err*shs1
          end do
       end do
       xo=xp(lxis+ll,n); sho=sum(xp(lxis+ll-4:lxis+ll,n)*(/3,-16,36,-48,25/))/12
       ! COMPUTE THE REST OF THE POINTS
       ip=lxis+ll;im=lxib-ll 
       call gridf(xp(:,n),pxi(:,n),xo,tmp,sho,she1,lxit,im,ip)
       do i=lxis+ll+1,lxie-1
          yp(i,n)=naca(xp(i,n),tmp,thk,n-1)
       end do
       yp(lxie,n)=zero
       ! ROTATE AND MOVE AEROFOILS
       do i = lxis, lxie
          oxp=xp(i,n);oyp=yp(i,n)
          xp(i,n) = (oxp*cos(alph)+oyp*sin(alph))+tmpa;
          yp(i,n) = (-oxp*sin(alph)+oyp*cos(alph))+tmpb;
       end do
    end do

!--VERICAL END BOUNDARIES
   !-Y-COORDINATE
   do n = 0, bkx
      if (n.ne.bkx) then !(all but right boundary)
         !-BLOCK0
         !-(2-0.5lwk(1))-2 LE curve-LE
         im=nvbl(0)
         ip=letse(0,1)-im;
         tmpa=py(2,n)-lvbl(0);sha=sml;tmpb=py(2,n);shb=shs2*cos(pi4+delt1)
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         !-0-1(2-0.5lwk(1)) Bottom->LE curve
         ip=letse(0,0); im=letbk(0)-im
         tmpa=py(0,n);sha=sml;tmpb=py(2,n)-lvbl(0);shb=qet(im,n)
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         if (n==bkx-1) ra2=qet(letse(0,0),n)
         !-BLOCK1
         !!-3-(3+lwk(1)) LE->LE curve
         ip=letse(1,0); im=nwk(1);
         tmpa=py(2,n);sha=shs2*cos(pi4-delt1);tmpb=py(2,n)+lvbl(1);shb=sml
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         !-(3+lwk(1))-5 LE curve->Top
         ip=ip+im; im=letbk(1)-im
         tmpa=tmpb;sha=qet(ip,n);tmpb=py(5,n);shb=sml;
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         if (n==bkx-1) ra1=qet(letse(1,1),n)
      else !(right boundary)
         !-BLOCK0
         !-(2-0.5lwk(1))-2 LE curve-LE
         im=half*nwk2(1)
         ip=letse(0,1)-im;
         tmpa=py(2,n)-lwk2(1);sha=sml;tmpb=py(2,n);shb=shs2*cos(pi4+delt1)*smod(1)
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         !-0-1(2-0.5lwk(1)) Bottom->LE curve
         ip=letse(0,0); im=letbk(0)-im
         tmpa=py(0,n);sha=ra2*0.9e0;tmpb=py(2,n)-lwk2(1);shb=qet(im,n)
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         !-BLOCK1
         !!-3-(3+lwk(1)) LE->LE curve
         ip=letse(1,0); im=nwk2(1)
         tmpa=py(2,n);sha=shs2*cos(pi4-delt1)*smod(1);ra0=sha
         tmpb=py(2,n)+lwk2(1);shb=sml
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         !-(3+lwk(1))-5 LE curve->Top
         ip=ip+im; im=lett-ip
         tmpa=tmpb;sha=qet(ip,n);tmpb=py(5,n);shb=ra1*1.2e0!ra0*32
         call gridf(yq(:,n),qet(:,n),tmpa,tmpb,sha,shb,lett,im,ip)
         if (k==0.and.myid==0) then
            write(*,*) 'Right boundary: mesh size ratio',qet(letse(1,1),n)/ra0
         end if
      end if
   end do

   !--X-COORDINATE for vertical lines
   !--Left Boundary
      xq(:,0)=px(0,0)
   !--Interface 1
      !-BLOCK0
      !-0-bl0
      is=letse(0,0);ie=letse(0,1)
      xq(is:ie,1)=px(2,0)
      !-bl0-1
      is=letbk(0)-nvbl(0);ie=letse(0,1);
      k01=px(2,0);k03=px(2,1)
      k02=vslo(1,0,0);k04=vslo(1,0,1)
      x0=yq(is,1);x1=yq(ie,1)
      do i = is, ie
         xq(i,1)=inter2(k01,k02,k03,k04,x0,x1,yq(i,1),0)
      end do
      !-BLOCK1
      !-bl1-3
      is=letse(1,0);;ie=letse(1,1)
      xq(is:ie,1)=px(2,3)
      !-2-bl1
      is=letse(1,0);ie=letse(1,0)+nvbl(1);
      k01=px(2,2);k03=px(2,3)
      k02=vslo(1,1,0);k04=vslo(1,1,1)
      x0=yq(is,1);x1=yq(ie,1)
      do i = is, ie
         xq(i,1)=inter2(k01,k02,k03,k04,x0,x1,yq(i,1),1)
      end do
   !--Interface 2
      !-BLOCK0
      is=letse(0,0);ie=letse(0,1);
      k01=px(3,0);k03=px(3,1)
      k02=vslo(2,0,0);k04=vslo(2,0,1)
      x0=yq(is,2);x1=yq(ie,2)
      do i = is, ie
         xq(i,2)=inter2(k01,k02,k03,k04,x0,x1,yq(i,2),0)
      end do
      !-BLOCK1
      is=letse(1,0);ie=letse(1,1)
      k01=px(3,2);k03=px(3,3)
      k02=vslo(2,1,0);k04=vslo(2,1,1)
      x0=yq(is,2);x1=yq(ie,2)
      do i = is, ie
         xq(i,2)=inter2(k01,k02,k03,k04,x0,x1,yq(i,2),1)
      end do
   !--Right Boundary
      xq(:,3)=px(5,0)


!--GRID INTERPOLATION
   allocate(linesy(0:bky-1)); linesy=(/0,2/)
   do l=0,bky-1
     js=letse(l,0); je=letse(l,1); n=linesy(l);
     do m=0,bkx-1
      is=lxise(m,0); ie=lxise(m,1)

      xx(is:ie,js)=xp(is:ie,n); xx(is:ie,je)=xp(is:ie,n+1)
      yy(is:ie,js)=yp(is:ie,n); yy(is:ie,je)=yp(is:ie,n+1)
      xx(is,js:je)=xq(js:je,m); xx(ie,js:je)=xq(js:je,m+1)
      yy(is,js:je)=yq(js:je,m); yy(ie,js:je)=yq(js:je,m+1)

      call hermite2D(is,ie,js,je,zs(k))
      !call interKim(is,ie,js,je,zs(k),n,m)
     end do
   end do
   deallocate(linesy)

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
       write(1,*) 'variables=x,y,z';
       write(1,"('zone i=',i4,' j=',i4,' k=',i4)") ie-is+1,je-js+1,1
      do j=js,je; do i=is,ie
        write(1,'(2es15.7)') xx(i,j),yy(i,j),zero
      end do; end do
      close(1)
   end if

   if(myid==0) then
      write(*,"('Grid generation is complete.')")
   end if

end if

   deallocate(xx,yy,zz,zs,xp,yp,xq,yq,pxi,qet,hslo,vslo)

 end subroutine gridaerofoil

 !==== 2D COONS' PATCH USING HERMITEAN INTERPOLATION ====
 subroutine hermite2D(is,ie,js,je,z)
 implicit none
 integer(k4), intent(in) :: is,ie,js,je
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
 
 
 allocate(xv(0:lxi,0:1),yv(0:lxi,0:1),&
          xu(0:let,0:1),yu(0:let,0:1),&
          xuv(0:1,0:1),yuv(0:1,0:1))
 
 ! compute extra derivatives at corners
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
 
 
 ! interpolate derivatives between corners
 lx=lxi
 do i = 0, lx
    xv(i,0) =hermite(xv(0,0),zero,xv(lx,0),zero,lx,i)
    xv(i,1) =hermite(xv(0,1),zero,xv(lx,1),zero,lx,i)
    yv(i,0) =hermite(yv(0,0),zero,yv(lx,0),zero,lx,i)
    yv(i,1) =hermite(yv(0,1),zero,yv(lx,1),zero,lx,i)
 end do

 ! interpolate derivatives between corners
 lx=let
 do i = 0, lx
    xu(i,0) =hermite(xu(0,0),zero,xu(lx,0),zero,lx,i)
    xu(i,1) =hermite(xu(0,1),zero,xu(lx,1),zero,lx,i)
    yu(i,0) =hermite(yu(0,0),zero,yu(lx,0),zero,lx,i)
    yu(i,1) =hermite(yu(0,1),zero,yu(lx,1),zero,lx,i)
 end do

 ! build X-interpolating matrix
 mx(1,:)=(/xx(is,js),xx(is,je),xv ( 0 ,0),xv( 0 ,1)/)
 mx(2,:)=(/xx(ie,js),xx(ie,je),xv (lxi,0),xv(lxi,1)/)
 mx(3,:)=(/xu(0 ,0 ),xu(let,0),0.0_k8    ,0.0_k8/)
 mx(4,:)=(/xu(0 ,1 ),xu(let,1),0.0_k8    ,0.0_k8/)
 
 ! build Y-interpolating matrix
 my(1,:)=(/yy(is,js),yy(is,je),yv ( 0 ,0),yv( 0 ,1)/)
 my(2,:)=(/yy(ie,js),yy(ie,je),yv (lxi,0),yv(lxi,1)/)
 my(3,:)=(/yu(0 ,0 ),yu(let,0),0.0_k8    ,0.0_k8/)
 my(4,:)=(/yu(0 ,1 ),yu(let,1),0.0_k8    ,0.0_k8/)
 
 ! Matrix multiplication
 do j = js, je; u=dble(j-js)/let;jj=j-js
    hv(:,1)=fillh(u) ! fill with hermitean funtions
    vx(:,1)=(/xx(is,j),xx(ie,j),xu(jj,0),xu(jj,1)/)
    vy(:,1)=(/yy(is,j),yy(ie,j),yu(jj,0),yu(jj,1)/)
    do i = is, ie; u=dble(i-is)/lxi;ii=i-is
       hu(1,:)=fillh(u) ! fill with hermitean funtions
       ux(1,:)=(/xx(i,js),xx(i,je),xv(ii,0),xv(ii,1)/)
       uy(1,:)=(/yy(i,js),yy(i,je),yv(ii,0),yv(ii,1)/)
       val1=matmul(hu,vx);val2=matmul(ux,hv);val3=matmul(hu,matmul(mx,hv))
       xx(i,j)=val1(1,1)+val2(1,1)-val3(1,1)
       val1=matmul(hu,vy);val2=matmul(uy,hv);val3=matmul(hu,matmul(my,hv))
       yy(i,j)=val1(1,1)+val2(1,1)-val3(1,1)
    end do
 end do
 
 zz(is:ie,js:je)=z ! no Z-interpolation, just copy value
    
 end subroutine hermite2D

!===== WRITE GRID LINES ====
!Write grid lines in Tecplot format
 subroutine writeLine(x,y,is,ie,cname)
   implicit none 
   integer(k4), intent(in) :: is,ie
   character(*),intent(in) :: cname
   real(k8), dimension(is:ie) :: x,y

   open(unit=10, file='out/line'//cname//'.dat')
   write(10,*) 'variables="x","y"';
   write(10,*) 'ZONE T= "'//cname//'"';
   write(10,"('i=',i4,' j=',i4)") ie-is+1,1
   do i=is,ie
     write(10,'(2es15.7)') x(i),y(i)
   end do
   close(10)
      
 end subroutine writeLine

 !==== Kim's interpolation scheme ====
 subroutine interKim(is,ie,js,je,z,n,m)
 implicit none
 integer(k4), intent(in) :: n,m,is,ie,js,je
 real(k8), intent(in) :: z
    
 ra0=abs(xp(is,n)-xp(is,n+1))-sml
 ra1=abs(xp(ie,n)-xp(ie,n+1))-sml
 if (ra0>0) then
    tmpa=1
 end if
 if (ra1>0) then
    tmpb=1
 end if
 select case(m)
 case(0); ii=0 ; nn=1; tmpa=0
 case(1); ii=0 ; nn=0
 case(2); ii=0 ; nn=1; tmpb=0
 end select
    gf=one
 do j=js,je; do i=is,ie
    pp=real(max(i-is-ii,0),kind=nr)/(ie-is-ii); qq=real(j-js,kind=nr)/(je-js)
    tmp=sin(halfpi*pp); ra0=half*((2-n)*real(je-j,kind=nr)/(je-js)+n*qq); ra1=gf+(two-gf)*ra0*ra0
    pxi(i,n)=(1-nn)*tmp**ra1+nn*tmp*tmp
    ts=tmpa*tmpb*(one-pxi(i,n))+one-tmpb; te=one-ts
    xx(i,j)=(xp(i,n+1)-xp(i,n))*(ts*(xq(j,m)-xq(js,m)+qq*(one-tmpa))/(xq(je,m)-xq(js,m)+one-tmpa)&
           +te*(xq(j,m+1)-xq(js,m+1)+qq*(one-tmpb))/(xq(je,m+1)-xq(js,m+1)+one-tmpb))+xp(i,n)
    ts=one-pxi(i,n); te=one-ts
    yy(i,j)=(yp(i,n+1)-yp(i,n))*(ts*(yq(j,m)-yq(js,m))/(yq(je,m)-yq(js,m))&
           +te*(yq(j,m+1)-yq(js,m+1))/(yq(je,m+1)-yq(js,m+1)))+yp(i,n)
    zz(i,j)=z
 end do; end do

 end subroutine interKim

!===== SUBROUTINE FOR GRID LINE GENERATION (constant size)
 subroutine cgridf(x,xxi,xo,xn,dxo,dxn,lxi,mxin,ip)

 integer(k4),intent(in) :: lxi,mxin,ip
 real(k8),dimension(0:lxi),intent(inout) :: x,xxi
 real(k8),intent(in) :: xo,dxo
 real(k8),intent(out) :: xn,dxn

 integer(k4) :: i,ii
 do i=0,mxin
    ii=i+ip;
    xn=xo+dxo*i
    x(ii)=xn
    xxi(ii)=dxo
 end do
 dxn=dxo

 end subroutine cgridf

!===== FUNCTION FOR HERMITIAN INTERPOLATION (integer entries)
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

 !==== Define Hermitean functions ====
 function fillh(u) result(v)
    real(k8), dimension(0:3) :: v
    real(k8) :: a0,a1,b0,b1
    real(k8), intent(in) :: u
 
    a0=+2*u**3 -3*u**2 +0*u**1 +1*u**0; !y0
    a1=-2*u**3 +3*u**2 +0*u**1 +0*u**0; !y1
    b0=+1*u**3 -2*u**2 +1*u**1 +0*u**0; !dydx0
    b1=+1*u**3 -1*u**2 +0*u**1 +0*u**0; !dydx1
 
    !a0=+0*u**3 +0*u**2 -1*u**1 +1*u**0; !y0
    !a1=+0*u**3 +0*u**2 +1*u**1 +0*u**0; !y1
    !b0=+0*u**3 +0*u**2 +0*u**1 +0*u**0; !dydx0
    !b1=+0*u**3 +0*u**2 +0*u**1 +0*u**0; !dydx1

    v(:)=(/a0,a1,b0,b1/)
 end function fillh

!===== NACA FUNCTION 00XX-Series
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
    if (n==0) then
       y=-1.0e0*y
    end if
 end function naca
 
!===== Hermitean interpolation FOR INTERFACE POINTS (real entries) 
 function inter(k01,k02,k03,k04,x0,x1,x) result(y)
 real(k8) :: y
 real(k8), intent(in) :: x,x0,x1,k01,k02,k03,k04  
 real(k8) :: a,b,c,d,l,invl,fx

 l=x1-x0
 invl=1.0_k8/l
 fx=(x-x0)*invl

    a=(+2*fx**3 -3*fx**2 +0*fx +1);
    b=(+1*fx**3 -2*fx**2 +1*fx +0);
    c=(-2*fx**3 +3*fx**2 +0*fx +0);
    d=(+1*fx**3 -1*fx**2 +0*fx +0);
    y=a*k01+b*l*k02+c*k03+d*l*k04;
 end function inter

 !==== Quintic hermitean interpolatin with 0-second derivative in one side
 function inter2(k01,k02,k03,k04,x0,x1,x,side) result(y)
 implicit none
 real(k8) :: y
 integer, intent(in) :: side
 real(k8), intent(in) :: x,x0,x1,k01,k02,k03,k04  
 real(k8) :: a,b,c,d,l,invl,fx

 l=x1-x0
 invl=1.0_k8/l
 fx=(x-x0)*invl

   selectcase(side)
   case(1)
    a=(-3*fx**4 +8*fx**3 -6*fx**2 +0*fx +1)
    b=(-1*fx**4 +3*fx**3 -3*fx**2 +1*fx +0);
    d=(-2*fx**4 +5*fx**3 -3*fx**2 +0*fx +0);
    c=(+3*fx**4 -8*fx**3 +6*fx**2 +0*fx +0);
   case(0)
    a=(+3*fx**4 -4*fx**3 +0*fx**2 +0*fx +1)
    b=(+2*fx**4 -3*fx**3 +0*fx**2 +1*fx +0);
    d=(+1*fx**4 -1*fx**3 +0*fx**2 +0*fx +0);
    c=(-3*fx**4 +4*fx**3 +0*fx**2 +0*fx +0);
   end select
   y=a*k01+b*l*k02+c*k03+d*l*k04;
 end function inter2

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

!===== Tangent function, entry in degrees
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

