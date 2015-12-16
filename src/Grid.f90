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

 real(k8),dimension(0:1,0:1) :: szth,wkth,blx,bly,dlth,szsh,wksh
 integer(k4),dimension(0:1,0:1) :: szll,wkll
 real(k8),dimension(0:1) :: blsh
 integer(k4),dimension(0:1) :: blll

 contains

!===== GRID GENERATION
!        <==Grid Sketch==>
!        0             1         2                        3
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
 real(k8) :: lbl,lwle
 real(k8) :: alph,thk=21
 real(k8) :: oxp,oyp
 real(k8) :: tmps,tmpe,tmpc
 real(k8) :: sha,shb,shc
 integer(k4) :: smod,npx,npy,mm,nn
 integer(k4), dimension(:), allocatable :: linesy,linesx
 real(k8), dimension(:), allocatable :: degarr
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
    shs1=ximod*smgrid; 
    shs2=etamod*smgrid; she1=shs2
    smod=2
    tmp=(shs2+smod*shs2)*half
    lbl=0.2_k8

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lzebk(0)))
    allocate(pxi(0:lxit),qet(0:lett))

    npy=min(nthick,1)+bky;npx=bkx
    allocate(xp(0:lxit,0:npy),yp(0:lxit,0:npy))
    allocate(xq(0:lett,0:npx),yq(0:lett,0:npx))
    allocate(px(0:bkx+2,0:npy),py(0:npy+2,0:npx))
    allocate(hslo(0:npy,0:bkx-1,0:1),vslo(0:npx,0:bky-1,0:1))

if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream',form='unformatted')

!---Domain Sizes
    dlth(0,0)=half*domlen; dlth(0,1)=domlen+szxt
    dlth(1,0)=half*domlen; dlth(1,1)=domlen+szxt
!---POINTS IN SPONGE
    szth(0,0)=szth1; szth(0,1)=szth2+szxt
    szth(1,0)=szth1; szth(1,1)=szth2+szxt
    szll(0,0)=min(szth(0,0)/dlth(0,0),0.15e0)*lxibk(0)
    szll(0,1)=min(szth(0,1)/dlth(0,1),0.15e0)*lxibk(bkx-1)
    szll(1,0)=min(szth(1,0)/dlth(1,0),0.15e0)*letbk(0)
    szll(1,1)=min(szth(1,1)/dlth(1,1),0.15e0)*letbk(bky-1)
    szsh(:,:)=szth(:,:)/szll(:,:)
    if (myid==0) then
       write(*,"('Sponge Pts: west=',i4,' north=',i4,' east=',i4, ' south=',i4)")&
       szll(0,0),szll(1,1),szll(0,1),szll(1,0)
    end if
!---Boundary Layer Refinement
    blx(0,0)=0.3*c1
    blx(1,0)=zero
    blx(0,1)=0.3*c1
    blx(1,1)=zero
    bly(0,0)=0.5e0*c1!thk*0.01*c1
    bly(1,0)=0.5e0*c1!thk*0.01*c1
    bly(0,1)=1.5e0*c1!thk*0.01*c1
    bly(1,1)=1.5e0*c1!thk*0.01*c1
    blll(0)=(letbk(0)-szll(1,0))*0.55e0
    blll(1)=(letbk(1)-szll(1,1))*0.55e0
    blsh(:)=bly(0,:)/blll(:)
    if (myid==0) then
       write(*,"('BL Pts: south=',i4,' north=',i4)")&
       blll(0),blll(1)
    end if
!---Wake Points
    wkth(0,0)=dlth(0,0)-szth(0,0)
    wkth(0,1)=dlth(0,1)-szth(0,1)-bly(0,1)
    wkth(1,:)=dlth(1,:)-szth(1,:)-bly(0,:)
    wkll(0,0)=lxibk(0)-szll(0,0)
    wkll(0,1)=lxibk(2)-szll(0,1)-blll(1)
    wkll(1,0)=letbk(0)-szll(1,0)-blll(0)
    wkll(1,1)=letbk(1)-szll(1,1)-blll(1)
    wksh(:,:)=wkth(:,:)/wkll(:,:)
    if (myid==0) then
       write(*,"('Wake Pts: west=',i4,' north=',i4,' east=',i4, ' south=',i4)")&
       wkll(0,0),wkll(1,1),wkll(0,1),wkll(1,0)
    end if
 do k=0,lzebk(0)
!---SPANWISE COORDINATE (UNIFORM)
    zs(k)=span*(real(lzebk(0)-k,k8)/lzebk(0)-half)
!---WAVY LEADING-EDGE PROFILE
    lwle=wlea*sin(2*pi*(zs(k)-zs(0))/wlew)
!---VERTICAL LINES
    px(0,:)=-dlth(0,0)
    px(1,:)=px(0,:)+szth(0,0)
    px(2,1:2)=-half*c1+lwle*cos(delt1);
    px(2,0)=px(2,2)-blx(0,0)
    px(2,3)=px(2,2)-blx(0,1)
    px(3,:)=-half*c1+c1*cos(delt1)
    px(4,:)=px(3,:)+dlth(0,1)-szth(0,1)
    px(5,:)=px(3,:)+dlth(0,1)
!---HORIZONTAL LINES
    py(0,:)=-dlth(1,0)
    py(1,:)=py(0,:)+szth(1,0)
    py(2,0)=zero
    py(2,1)=-lwle*sin(delt1)
    py(2,2)=-c1*sin(delt1)
    py(2,3)=py(2,2)
    py(3,:)=py(2,:)
    py(4,:)=py(2,2)+dlth(1,1)-szth(1,1)
    py(5,:)=py(2,2)+dlth(1,1)

!----- CONSTANT ANGLES IN RADIANS
    if(.not.allocated(degarr)) allocate(degarr(3))
    degarr(:)=(/25_k8,10_k8,0_k8/); degarr(:)=degarr(:)*pi/180_k8
!----- INITIAL AND END HORIZONTAL SLOPES
    hslo(0,:,:)=zero
    hslo(1,0,:)=(/zero,-tan(delt1)/)
    hslo(1,1,:)=zero
    hslo(1,2,:)=zero
    hslo(2,:,:)=hslo(1,:,:)
    hslo(3,:,:)=zero

!----- INITIAL AND END VERTICAL SLOPES
    vslo(0,:,:)=zero
    vslo(1,0,:)=(/zero,tan(pi4+delt1)/)
    vslo(1,1,:)=(/-tan(pi4-delt1),zero/)
    vslo(2,:,:)=zero
    vslo(3,:,:)=zero


!--HORIZONTAL INTERFACES
   !--X-COORDINATE
   do n = 0,npy
   !--BLOCK0
   !--0-1 Sponge
      ip=lxise(0,0); im=szll(0,0);
      tmpa=px(0,n);sha=szsh(0,0);tmpb=px(1,n);shb=wksh(0,0)
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--1-2
      ip=ip+im; im=lxibk(0)-im;
      tmpa=px(1,n);sha=shb;tmpb=px(2,n);shb=shs1
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   if ((n.ne.1).or.(n.ne.2)) then
      !--BLOCK1
      !--2-3
         ip=lxise(1,0); im=lxibk(1);
         tmpa=px(2,n);sha=shs1;tmpb=px(3,n);shb=she1
         call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   end if
   !--BLOCK2
   !--3-bl Refinement
      ip=lxise(2,0); im=blll(1)
      tmpa=px(3,n);sha=she1;tmpb=px(3,n)+bly(0,1);shb=sml
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--bl-4 Wake
      ip=lxise(2,0)+im; im=lxibk(2)-szll(0,1)-blll(1);
      tmpa=tmpb;sha=pxi(ip);tmpb=px(4,n);shb=wksh(0,1)
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   !--4-5 Sponge
      ip=ip+im; im=szll(0,1);
      tmpa=px(4,n);sha=shb;tmpb=px(5,n);shb=szsh(0,1)
      call gridf(xp(:,n),pxi,tmpa,tmpb,sha,shb,lxit,im,ip)
   end do

   !--Y-COORDINATE
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
         yp(i,1)=inter(k01,k02,k03,k04,x0,x1,xp(i,1))
         yp(i,2)=yp(i,1) ! Copy interface
      end do
   !--Upper Boundarie
      yp(:,3)=py(5,0)

!----- AEROFOIL SURFACE GRID POINTS
    tmp=c1-lwle;tmpa=px(2,1);tmpb=py(2,1);
    lxis=lxise(1,0);lxie=lxise(1,1);lxib=lxibk(1);
    alph=delt1;thk=21
    ! DETERMINE UPPER AND LOWER SIDES
    do n=1,2
       yp(lxis,n)=zero;xp(lxis,n)=zero
       ! DETERMINE THE FIRST LL POINTS
       ll=8 ! "LL" MUST BE EQUAL TO OR LARGER THAN 4.
       do i=lxis+1,lxis+ll
          xp(i,n)=xp(i-1,n)+half*shs1; err=1
          do while(abs(err)>sml)
             yp(i,n)=naca(xp(i,n),tmp,thk,n-1)
             err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shs1-1;
             xp(i,n)=xp(i,n)-half**5*err*shs1
          end do
       end do
       xo=xp(lxis+ll,n); sho=sum(xp(lxis+ll-4:lxis+ll,n)*(/3,-16,36,-48,25/))/12
       ! COMPUTE THE REST OF THE POINTS
       ip=lxis+ll;im=lxib-ll 
       call gridf(xp(:,n),pxi,xo,tmp,sho,she1,lxit,im,ip)
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

      is=lxise(0,0);ie=lxise(0,1);n=1
      call writeLine(xp(is:ie,n),yp(is:ie,n),is,ie,'0n')
      is=lxise(0,0);ie=lxise(0,1);n=2
      call writeLine(xp(is:ie,n),yp(is:ie,n),is,ie,'3s')
      is=lxise(0,0);ie=lxise(0,1);n=3
      call writeLine(xp(is:ie,n),yp(is:ie,n),is,ie,'3n')
      is=lxise(1,0);ie=lxise(1,1);n=2
      call writeLine(xp(is:ie,n),yp(is:ie,n),is,ie,'4s')
      
!--VERICAL END BOUNDARIES
   !-Y-COORDINATE
   do n = 0, bkx
      !-BLOCK0
      !-0-1
      ip=letse(0,0); im=szll(1,0);
      tmpa=py(0,n);sha=szsh(1,0);tmpb=py(1,n);shb=wksh(1,0)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lett,im,ip)
      !-1-bl0
      ip=ip+im; im=letbk(0)-(szll(1,0)+blll(0));
      tmpa=py(1,n);sha=shb;tmpb=py(2,n)-bly(0,0);shb=blsh(0)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lett,im,ip)
      !-bl0-2 refinement
      ip=ip+im; im=blll(0);
      tmpa=tmpb;sha=shb;tmpb=py(2,n);shb=shs2*cos(pi4+delt1)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lett,im,ip)
      !-BLOCK1
      !-3-bl1 refinement
      ip=letse(1,0); im=blll(1);
      tmpa=py(2,n);sha=shs2*cos(pi4-delt1);tmpb=py(2,n)+bly(0,1);shb=sml
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lett,im,ip)
      !-bl1-4
      ip=ip+im; im=letbk(1)-(szll(1,1)+blll(1));
      tmpa=tmpb;sha=qet(ip);tmpb=py(4,n);shb=wksh(1,1)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lett,im,ip)   
      !-4-5
      ip=ip+im; im=szll(1,1)
      tmpa=py(4,n);sha=shb;tmpb=py(5,n);shb=szsh(1,1)
      call gridf(yq(:,n),qet,tmpa,tmpb,sha,shb,lett,im,ip)
   end do

   !--X-COORDINATE
   !--Left Boundary
      xq(:,0)=px(0,0)
   !--Interface 1
      !-BLOCK0
      !-0-bl0
      is=letse(0,0);ie=letbk(0)-blll(0)
      xq(is:ie,1)=px(2,0)
      !-bl0-1
      k01=px(2,0);k03=px(2,1)
      k02=vslo(1,0,0);k04=vslo(1,0,1)
      x0=py(2,1)-bly(0,0);x1=py(2,1)
      is=ie;ie=letse(0,1);
      do i = is, ie
         xq(i,1)=inter2(k01,k02,k03,k04,x0,x1,yq(i,1),0)
      end do
      !-BLOCK1
      !-2-bl1
      k01=px(2,2);k03=px(2,3)
      k02=vslo(1,1,0);k04=vslo(1,1,1)
      x0=py(3,1);x1=py(3,1)+bly(0,1)
      is=letse(1,0);ie=letse(1,0)+blll(1);
      do i = is, ie
         xq(i,1)=inter2(k01,k02,k03,k04,x0,x1,yq(i,1),1)
      end do
      !-bl1-3
      is=ie+1;ie=letse(1,1)
      xq(is:ie,1)=px(2,3)
   !--Interface 2
      xq(:,2)=px(3,0)
   !--Right Boundary
      xq(:,3)=px(5,0)

      is=letse(1,0);ie=letse(1,1);n=0
      call writeLine(xq(is:ie,n),yq(is:ie,n),is,ie,'3w')
      is=letse(1,0);ie=letse(1,1);n=1
      call writeLine(xq(is:ie,n),yq(is:ie,n),is,ie,'3e')

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

   deallocate(xx,yy,zz,zs,xp,yp,xq,yq,pxi,qet)

 end subroutine gridaerofoil

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
    xv(i,0) =hermite(xv(0,0),zero,xv(lx,0),zero,lx,i)
    xv(i,1) =hermite(xv(0,1),zero,xv(lx,1),zero,lx,i)
    yv(i,0) =hermite(yv(0,0),zero,yv(lx,0),zero,lx,i)
    yv(i,1) =hermite(yv(0,1),zero,yv(lx,1),zero,lx,i)
 end do

 lx=let
 do i = 0, lx
    xu(i,0) =hermite(xu(0,0),zero,xu(lx,0),zero,lx,i)
    xu(i,1) =hermite(xu(0,1),zero,xu(lx,1),zero,lx,i)
    yu(i,0) =hermite(yu(0,0),zero,yu(lx,0),zero,lx,i)
    yu(i,1) =hermite(yu(0,1),zero,yu(lx,1),zero,lx,i)
 end do

 mx(1,:)=(/xx(is,js),xx(is,je),xv ( 0 ,0),xv( 0 ,1)/)
 mx(2,:)=(/xx(ie,js),xx(ie,je),xv (lxi,0),xv(lxi,1)/)
 mx(3,:)=(/xu(0 ,0 ),xu(let,0),0.0_k8    ,0.0_k8/)
 mx(4,:)=(/xu(0 ,1 ),xu(let,1),0.0_k8    ,0.0_k8/)
 
 my(1,:)=(/yy(is,js),yy(is,je),yv ( 0 ,0),yv( 0 ,1)/)
 my(2,:)=(/yy(ie,js),yy(ie,je),yv (lxi,0),yv(lxi,1)/)
 my(3,:)=(/yu(0 ,0 ),yu(let,0),0.0_k8    ,0.0_k8/)
 my(4,:)=(/yu(0 ,1 ),yu(let,1),0.0_k8    ,0.0_k8/)
 
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
    if (n==0) then
       y=-1.0e0*y
    end if
 end function naca
 
!===== FUNCTIONS FOR INTERFACE POINTS
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

