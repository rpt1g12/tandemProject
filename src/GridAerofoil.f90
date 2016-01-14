!*****
!***** 3D FLAT-PLATE GRID GENERATION
!*****

 module gridgen

 use subroutineso
 implicit none

 integer(kind=ni),parameter :: lnaca=1000
 integer(kind=ni) :: lxi0,lxi1,lxi2,let0,lze0
 integer(kind=ni) :: lxit,lett,lxie0,lxis1,lxie1,lxis2,lete0,lets1,lxisz,im,jm

 real(kind=nr),dimension(0:lnaca,2) :: xnaca,ynaca

 real(kind=nr),dimension(:,:),allocatable :: xx,yy,zz
 real(kind=nr),dimension(:),allocatable :: xyzmb,zs
 real(kind=nr),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(kind=nr),dimension(:),allocatable :: pxi,qet

 real(kind=nr) :: rs,re,rp,ts,te,shs,she,shswle
 real(kind=nr) :: xa,xb,xc,xd,xe,xo,ya,yb,yc,yd,yo,sho,pp,qq
 real(kind=nr) :: am,err,tmp,tmpa,tmpb,gf

 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,litr,smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb)

 integer(kind=ni),intent(in) :: ngridv,nthick,litr
 real(kind=nr),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb

    lxit=lxi0+lxi1+lxi2+2; lett=2*let0+1
    lxie0=lxi0; lxis1=lxie0+1; lxie1=lxis1+lxi1; lxis2=lxie1+1
    lete0=let0; lets1=lete0+1

    rs=nthick*domlen/(domlen+three); re=zero
    shs=smgrid; she=shs

    np=3*(lxio+1)*(leto+1)-1
    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett))
    allocate(xp(0:lxit,0:3),yp(0:lxit,0:3),xq(0:lett,0:3),yq(0:lett,0:3))
    allocate(xyzmb(0:np),zs(0:lze0),pxi(0:lxit),qet(0:lett))

!----- WAVY LEADING-EDGE PROFILE

 if(myid==mo(mb)) then
    open(9,file=cgrid); close(9,status='delete')
    open(9,file=cgrid,access='direct',form='unformatted',recl=nrecd*(np+1))

 do k=0,lze0
    zs(k)=span*(real(lze0-k,kind=nr)/lze0-half)

    xa=-domlen; xb=-half+wlea*sin(twopi*(zs(k)-zs(0))/wlew); xc=half; xd=domlen-szth1; xe=domlen+szxt
    ya=-domlen; yb=zero; yc=zero; yd=domlen

    fctr=twopi/wlew; shswle=shs*sqrt(one+0*(fctr*wlea*cos(fctr*(zs(k)-zs(0))))**two)

!----- AEROFOIL SURFACE GRID POINTS

 if(rs==zero) then
    tmp=(xc-xb)/lnaca
 do n=1,2; do i=0,lnaca
    xnaca(i,n)=i*tmp+xb; ynaca(i,n)=zero
 end do; end do
 else
    tmp=xc-xb
    open(8,file='aerofoil.dat')
 do n=1,2; do i=0,lnaca
    read(8,*) xnaca(i,n),ynaca(i,n)
    xnaca(i,n)=tmp*xnaca(i,n)+xb; ynaca(i,n)=tmp*ynaca(i,n)
 end do; end do
    close(8)
 end if

 do n=1,2
    yp(lxis1,n)=yb
 if(rs==zero) then
    ll=0; xo=xb; sho=shswle
 else
    ll=8 ! "ll" must be equal to or larger than 4.
    xp(lxis1,n)=xb
 do i=lxis1+1,lxis1+ll
    xp(i,n)=xp(i-1,n)+half*shswle; err=one
 do while(abs(err)>sml)
    yp(i,n)=ylagi(i,n)
    err=sqrt((xp(i,n)-xp(i-1,n))**two+(yp(i,n)-yp(i-1,n))**two)/shswle-one; xp(i,n)=xp(i,n)-half*err*shswle
 end do
 end do
    xo=xp(lxis1+ll,n); sho=sum(xp(lxis1+ll-4:lxis1+ll,n)*(/3.0_nr,-16.0_nr,36.0_nr,-48.0_nr,25.0_nr/))/12.0_nr
 end if
    ip=lxis1+ll; im=lxi1-ll; call gridf(xp(:,n),pxi,xo,xc,sho,she,lxit,im,ip)
 do i=lxis1+ll+1,lxie1-1
    yp(i,n)=ylagi(i,n)
 end do
    yp(lxie1,n)=yc
 end do

!----- HORIZONTAL INTERFACE POINTS IN XI-DIRECTION

    n=1
    sho=tla/litr; ll=lxi0-two*(lxi0*sho-(-half-xa))/(sho-shswle+sml); xo=xa+ll*sho
    ip=0; im=ll; call gridf(xp(:,n),pxi,xa,xo,sho,sho,lxit,im,ip)
 if(ll/=lxi0) then
    ip=ll; im=lxi0-ll; call gridf(xp(:,n),pxi,xo,xb,sho,shswle,lxit,im,ip)
 end if
 if(k==0) then
    lp=ll; lh=minloc(abs(xa+szth1-xp(0:lxi0,n)),1)-1; lxisz=lh*lxi2/lxi0
 end if
    ip=lxis2; im=lxi2-lxisz; call gridf(xp(:,n),pxi,xc,xd,she,free,lxit,im,ip)
    ip=ip+im; im=lxisz; call gridf(xp(:,n),pxi,xd,xe,pxi(ip),free,lxit,im,ip)
 do m=1,2
 select case(m)
 case(1); is=0; ie=lxie0; ii=is; xo=xb-xa; yo=yb; case(2); is=lxis2; ie=lxit; ii=ie; xo=xd-xc; yo=yc
 end select
    am=yo/sin(halfpi)**two
 do i=is,ie
    gf=(-1)**(m+1)*(xp(i,n)-xp(ii,n))/xo; yp(i,n)=am*sin(halfpi*gf)**two
 end do
 do i=is,ie
    xp(i,n+1)=xp(i,n); yp(i,n+1)=yp(i,n)
 end do
 end do

!----- TOP & BOTTOM BOUNDARY POINTS IN XI-DIRECTION

    n=0
 if(rs==zero) then
 do i=0,lxit
    xp(i,n)=xp(i,n+1)
 end do
 else
 do i=0,lh
    xp(i,n)=xp(i,n+1)
 end do
    tmpa=xb-rs; tmpb=xc+re; xo=xp(lh,n)
    ip=lh; im=lxi0-lh; call gridf(xp(:,n),pxi,xo,tmpa,sho,two*shs,lxit,im,ip)
    ip=lxis1; im=lxi1; call gridf(xp(:,n),pxi,tmpa,tmpb,two*shs,she,lxit,im,ip)
    ip=lxis2; im=lxi2-lxisz; call gridf(xp(:,n),pxi,tmpb,xd,she,free,lxit,im,ip)
    ip=ip+im; im=lxisz; call gridf(xp(:,n),pxi,xd,xe,pxi(ip),free,lxit,im,ip)
 end if
 do i=0,lxit
    yp(i,n)=ya; xp(i,n+3)=xp(i,n); yp(i,n+3)=yd
 end do

!----- VERTICAL INTERFACE POINTS IN ETA-DIRECTION

    yo=tlb*(2.5_nr-cutlb); ll=yo/(half*(shs+sho))
 do m=1,2
    if(rs==zero) then; fctr=one; else; fctr=sqrt(m*half); end if
    jp=lets1; jm=ll; call gridf(yq(:,m),qet,zero,yo,fctr*shs,sho,lett,jm,jp)
    jp=jp+jm; jm=let0-ll; call gridf(yq(:,m),qet,yo,yd,sho,free,lett,jm,jp)
 do j=0,lete0
    yq(j,m)=-yq(lett-j,m)
 end do
 end do
 if(rs==zero) then
 do j=0,lett
    xq(j,1)=xb; xq(j,2)=xc
 end do
 else
 do n=1,2
 select case(n); case(1); js=0; je=lete0; nn=0; case(2); js=lets1; je=lett; nn=3; end select
 do m=1,2
 select case(m)
 case(1); i=lxis1
    tmpa=sqrt((xp(i,n)-xp(i-2,n))**two+(yp(i,n)-yp(i-2,n))**two)
    tmpb=sqrt((xp(i,n)-xp(i+1,n))**two+(yp(i,n)-yp(i+1,n))**two)
    tmp=pi-half*acos(half*(tmpa*tmpa+tmpb*tmpb-(xp(i+1,n)-xp(i-2,n))**two-(yp(i+1,n)-yp(i-2,n))**two)/(tmpa*tmpb))
 case(2); i=lxie1
    tmpa=sqrt((xp(i,n)-xp(i-1,n))**two+(yp(i,n)-yp(i-1,n))**two)
    tmpb=sqrt((xp(i,n)-xp(i+2,n))**two+(yp(i,n)-yp(i+2,n))**two)
    tmp=half*acos(half*(tmpa*tmpa+tmpb*tmpb-(xp(i+2,n)-xp(i-1,n))**two-(yp(i+2,n)-yp(i-1,n))**two)/(tmpa*tmpb))
 end select
    res=tan(tmp); tmpa=abs(xp(i,nn)-xp(i,n)); tmpb=abs(yp(i,nn)-yp(i,n))
    qq=(half*res*tmpa)**two/abs(tmpb-res*tmpa); pp=two*sqrt(qq)/res
 do j=js,je
    xq(j,m)=(-1)**m*pp*(sqrt(abs(yq(j,m)-yp(i,n))+qq)-sqrt(qq))+xp(i,n)
 end do
 end do
 end do
 end if

!----- LEFT & RIGHT BOUNDARY POINTS IN ETA-DIRECTION

    yo=tlb*(2.5_nr-cutlb); fctr=0.8_nr; ll=yo/(half*(fctr*sho+sho))
 do m=0,3,3
    jp=lets1; jm=ll; call gridf(yq(:,m),qet,zero,yo,fctr*sho,sho,lett,jm,jp)
    jp=jp+jm; jm=let0-ll; call gridf(yq(:,m),qet,yo,yd,sho,free,lett,jm,jp)
 do j=0,lete0
    yq(j,m)=-yq(lett-j,m)
 end do
 end do
 do j=0,lett
    xq(j,0)=xa; xq(j,3)=xd
 end do
 if(k==0) then
    lq=ll
 end if

!----- GRID OUTPUT

 do n=0,2,2
 select case(n); case(0); js=0; je=lete0; case(2); js=lets1; je=lett; end select
 do m=0,2
    tmpa=int(rs/abs(rs-sml),kind=ni); tmpb=int(re/abs(re-sml),kind=ni)
 select case(m)
 case(0); is=0; ie=lxie0; ii=lh; nn=1; tmpa=0
 case(1); is=lxis1; ie=lxie1; ii=0; nn=0
 case(2); is=lxis2; ie=lxit; ii=0; nn=1; tmpb=0
 end select
    gf=one
 do j=js,je; do i=is,ie
    pp=real(max(i-is-ii,0),kind=nr)/(ie-is-ii); qq=real(j-js,kind=nr)/(je-js)
    tmp=sin(halfpi*pp); ra0=half*((2-n)*real(je-j,kind=nr)/(je-js)+n*qq); ra1=gf+(two-gf)*ra0*ra0
    pxi(i)=(1-nn)*tmp**ra1+nn*tmp*tmp
    ts=tmpa*tmpb*(one-pxi(i))+one-tmpb; te=one-ts
    xx(i,j)=(xp(i,n+1)-xp(i,n))*(ts*(xq(j,m)-xq(js,m)+qq*(one-tmpa))/(xq(je,m)-xq(js,m)+one-tmpa)&
           +te*(xq(j,m+1)-xq(js,m+1)+qq*(one-tmpb))/(xq(je,m+1)-xq(js,m+1)+one-tmpb))+xp(i,n)
    ts=one-pxi(i); te=one-ts
    yy(i,j)=(yp(i,n+1)-yp(i,n))*(ts*(yq(j,m)-yq(js,m))/(yq(je,m)-yq(js,m))&
           +te*(yq(j,m+1)-yq(js,m+1))/(yq(je,m+1)-yq(js,m+1)))+yp(i,n)
    zz(i,j)=zs(k)
 end do; end do
 end do
 end do
 select case(mb); case(0,1,2); js=0; je=lete0; case(3,4,5); js=lets1; je=lett; end select
 select case(mb)
 case(0,3); is=0; ie=lxie0; case(1,4); is=lxis1; ie=lxie1; case(2,5); is=lxis2; ie=lxit
 end select
    l=-3
 do j=js,je; do i=is,ie; l=l+3
    xyzmb(l:l+2)=(/xx(i,j),yy(i,j),zz(i,j)/)
 end do; end do
    write(9,rec=k+1) xyzmb(:)
 end do

    close(9)

!----- GRID CHECKING

 if(ngridv==1) then
    open(8,file='misc/gridview'//cnzone//'.dat')
    write(8,*) 'variables=v1,v2'; write(8,"('zone i=',i4,' j=',i4)") ie-is+1,je-js+1
 do j=js,je; do i=is,ie
    write(8,'(2es15.7)') xx(i,j),yy(i,j)
 end do; end do
    close(8)
 else
    open(8,file='misc/gridview'//cnzone//'.dat'); close(8,status='delete')
 end if
 if(myid==0) then
    write(*,"('Grid generation is complete.')")
    write(*,"('Number of cells to resolve inflow turbulence in x:',i4)") lp
    write(*,"('Number of cells to resolve inflow turbulence in y:',i4)") lq
    write(*,"('Number of cells across sponge:',i4)") minloc(abs(xa+szth1-xp(:,1)),1)-1
 end if

 end if
    deallocate(xx,yy,zz,xyzmb,xp,yp,xq,yq,pxi,qet)

 end subroutine gridaerofoil

!===== LAGRANGIAN INTERPOLATION FOR Y-COORDINATES

 function ylagi(i,n) result(yi)

 integer(kind=ni),intent(in) :: i,n
 real(kind=nr) :: yi

    is=0; ie=lnaca; ii=minloc(abs(xnaca(:,n)-xp(i,n)),1)-1; ip=ii
    if(ii-is<=1) then; ip=is+2; end if
    if(ie-ii<=1) then; ip=ie-2; end if
    yi=zero; alag(:)=xp(i,n)-xnaca(ip-2:ip+2,n)
 do jj=-2,2
    blag(:)=xnaca(ip+jj,n)-xnaca(ip-2:ip+2,n); ao=one; bo=one
 do ii=-2,2; if(ii/=jj) then
    ao=ao*alag(ii); bo=bo*blag(ii)
 end if; end do
    yi=yi+ao*ynaca(ip+jj,n)/bo
 end do

 end function ylagi

!=====

 end module gridgen

!*****