!*****
!***** 3D FLAT-PLATE GRID GENERATION
!*****

 module gridgen

 use subroutineso
 implicit none

 integer,parameter :: lnaca=1000
 integer :: lxi0,lxi1,lxi2,let0,lze0
 integer :: lxi3,lxi4
 integer :: lxit,lett,lxie0,lxis1,lxie1,lxis2,lete0,lets1,lxisz,im,jm
 integer :: lxie2,lxis3,lxie3,lxis4
 integer :: lxis,lxie,lxib
 integer :: lets,lete,letb

 real(nr),dimension(0:lnaca,2,2) :: xnaca,ynaca

 real(nr),parameter :: pi4=pi/4.0_nr

 real(nr),dimension(:,:),allocatable :: xx,yy,zz
 real(nr),dimension(:),allocatable :: zs
 real(nr),dimension(:,:),allocatable :: xp,yp,xq,yq
 real(nr),dimension(:),allocatable :: pxi,qet

 real(nr) :: rs,re,rp,ts,te,shs1,she1,shs2,she2,shs,she
 real(nr) :: xa,xb,xc,xd,xe,xo,xjct,ya,yb,yc,yd,yo,yjct,sho,pp,qq
 real(nr) :: xf,xg
 real(nr) :: am,err,tmp,tmpa,tmpb,gf
 real(nr) :: k1,k2,k3,k4,x0,x1

 contains

!===== GRID GENERATION

 subroutine gridaerofoil(ngridv,nthick,litr,smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb,gap,c1,c2,delt1,delt2,ximod,etamod)

 integer,intent(in) :: ngridv,litr,nthick,cutlb
 real(nr),intent(in) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,gap,c1,c2,delt1,delt2
 real(nr),intent(in) :: ximod,etamod
 real(nr) :: alph
 real(nr) :: oxp,oyp
 real(nr) :: tmps,tmpe,tmpc

    lxit=lxi0+lxi1+lxi2+lxi3+lxi4+4; lett=2*let0+1
    lxie0=lxi0; lxis1=lxie0+1; lxie1=lxis1+lxi1; lxis2=lxie1+1
    lxie2=lxis2+lxi2; lxis3=lxie2+1; lxie3=lxis3+lxi3; lxis4=lxie3+1
    lete0=let0; lets1=lete0+1

    rs=1*domlen/(domlen+3);
    re=0
    shs=smgrid; she=shs
    shs1=ximod*smgrid; she1=shs1
    shs2=etamod*smgrid; she2=shs2

    allocate(xx(0:lxit,0:lett),yy(0:lxit,0:lett),zz(0:lxit,0:lett),zs(0:lze0))
    allocate(xp(0:lxit,0:3),yp(0:lxit,0:3))
    allocate(xq(0:lett,0:5),yq(0:lett,0:5))
    allocate(pxi(0:lxit),qet(0:lett))

!----- WAVY LEADING-EDGE PROFILE

 if(myid==mo(mb)) then
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    open(1,file='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat',access='stream')

 do k=0,lze0
    zs(k)=span*(real(lze0-k,nr)/lze0-half)

!----- BLOCKS' BOUNDARIES
    xa=-domlen; xb=-half*c1; xc=half*c1;
    xd=xc+gap+wlea*sin(2*pi*(zs(k)-zs(0))/wlew); xe=xc+gap+c2
    xf=xe-half*c2+(domlen-szth1); xg=xf+szth1+szxt
    ya=-domlen; yb=0; yc=0; yd=domlen

!----- AEROFOIL SURFACE GRID POINTS

do m=1,2
select case(m)
case(1);tmp=c1;tmpa=xb;tmpb=xc;tmpc=0_nr;lxis=lxis1;lxie=lxie1;lxib=lxi1;alph=delt1
case(2);tmp=xe-xd;tmpa=xc+gap;tmpb=xe;tmpc=xd-tmpa;lxis=lxis3;lxie=lxie3;lxib=lxi3;alph=delt2
end select
 if(rs==0) then
    tmp=(tmpb-tmpa)/lnaca
 do n=1,2; do i=0,lnaca
    xnaca(i,n,m)=i*tmp+tmpa; ynaca(i,n,m)=0   ! create naca x coordinates
 end do; end do
 else
    open(2,file='aerofoil.dat')
 do n=1,2; do i=0,lnaca
    read(2,*) xnaca(i,n,m),ynaca(i,n,m)   ! read coordinates from file
    !xnaca(i,n,m)=tmp*xnaca(i,n,m); ynaca(i,n,m)=tmp*ynaca(i,n,m)
 end do; end do
    close(2)
 end if

 do n=1,2
    yp(lxis,n)=yb   ! yb = 0

    if(rs==0) then
       ll=0; xo=0.0e0; sho=shs1
    else
       ll=8 ! "ll" must be equal to or larger than 4.
       xp(lxis,n)=0.0e0
       do i=lxis+1,lxis+ll
          xp(i,n)=xp(i-1,n)+half*shs1; err=1
          do while(abs(err)>sml)
             yp(i,n)=ylagi(i,n,m)
             err=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)/shs1-1;
             xp(i,n)=xp(i,n)-half*err*shs1
          end do
       end do
       xo=xp(lxis+ll,n); sho=sum(xp(lxis+ll-4:lxis+ll,n)*(/3,-16,36,-48,25/))/12
    end if
   
    am=2; xjct=half*(xo+1.0e0); err=1
    
    do while(abs(err)>sml)
       ip=lxis+ll; im=(lxib-ll)/2; call ogridf(xp(:,n),pxi,xo,xjct,sho,am,0,lxit,im/2,im,ip)
       ip=ip+im; im=lxib-ll-im; call ogridf(xp(:,n),pxi,xjct,1.0_nr,pxi(ip),am,0,lxit,im/2,im,ip)
       err=pxi(ip+im)/she1-1; xjct=xjct+half*err*she1
    end do
    do i=lxis+ll+1,lxie-1
       yp(i,n)=ylagi(i,n,m)
    end do
    yp(lxie,n)=yc

    do i = lxis, lxie
    oxp=xp(i,n);oyp=yp(i,n)
    xp(i,n) = tmp*(oxp*cos(alph)+oyp*sin(alph))+tmpa+tmpc*cos(alph);
    yp(i,n) = tmp*(-oxp*sin(alph)+oyp*cos(alph))-tmpc*sin(alph);
    end do
    
 end do
end do

!----- HORIZONTAL INTERFACE POINTS IN XI-DIRECTION

    n=1
    sho=tla/litr; ll=2*litr!0!lxi0-2*(lxi0*sho-(-half-xa))/(sho-shs)
    am=2; xo=xa+ll*sho
    ip=0; im=ll; call ogridf(xp(:,n),pxi,xa,xo,sho,am,0,lxit,im/2,im,ip)
    xjct=half*(xo+xb); err=1
 do while(abs(err)>sml)
    ip=ll; im=(lxi0-ll)/2; call ogridf(xp(:,n),pxi,xo,xjct,sho,am,0,lxit,im/2,im,ip)
    ip=ip+im; im=lxi0-ll-im; call ogridf(xp(:,n),pxi,xjct,xb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/shs1-1; xjct=xjct+half*err*shs1
 end do
 if(k==0) then
    lxisz=lxi4*(minloc(abs(xa+szth1-xp(0:lxi0,n)),1)-1)/lxi0; lp=ll     ! lxi2 was
                                                                        !used before
 end if
 tmpa=xp(lxie1,n); tmpb=xp(lxis3,n)
 xjct=(tmpa+tmpb)*half; err=1
 do while(abs(err)>sml)
    ip=lxis2; im=lxi2/2; call ogridf(xp(:,n),pxi,tmpa,xjct,she1,am,0,lxit,im/2,im,ip)
    ip=ip+im; im=lxi2-im; call ogridf(xp(:,n),pxi,xjct,tmpb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/she1-1; xjct=xjct+half*err*shs1
 end do
    tmpa=xp(lxie3,n); tmpb=xf
    ip=lxis4; im=lxi4-lxisz; call ogridf(xp(:,n),pxi,tmpa,tmpb,she1,am,0,lxit,im/2,im,ip)
    tmpa=xf; tmpb=xg
    xjct=(tmpa+tmpb)*half; err=1
 do while(abs(err)>sml)
    ip=lxis4+(lxi4-lxisz); im=lxisz/2; call ogridf(xp(:,n),pxi,tmpa,xjct,pxi(ip),am,0,lxit,im,im,ip)
    ip=ip+im; im=lxisz-im; call ogridf(xp(:,n),pxi,xjct,tmpb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/sho-1; xjct=xjct+half*err*sho
 end do
 do m=1,3
 select case(m)
 case(1); is=0; ie=lxie0; k1=yc; k2=0.0_nr; k3=yp(lxis1,n); k4=-tan(delt1);
          x0=xa; x1=xp(lxis1,n)
 case(2); is=lxis2; ie=lxie2; k1=yp(lxie1,n); k2=-tan(delt1); k3=yp(lxis3,n); k4=-tan(delt2);
          x0=xp(lxie1,n); x1=xp(lxis3,n)
 case(3); is=lxis4; ie=lxit; k1=yp(lxie3,n); k2=-tan(delt2); k3=yc; k4=0.0_nr
          x0=xp(lxie3,n); x1=xg
 !case(1); is=0; ie=lxie0; ii=is; xo=xb-xa; yo=yb
 !case(2); is=lxis2; ie=lxie2; ii=is; xo=xp(lxis3,n)-xp(lxie1,n); yo=yp(lxis3,n)-yp(lxie1,n)
 !case(3); is=lxis4; ie=lxit; ii=ie; xo=xg-xp(lxie3,n); yo=yp(lxie3,n)
 end select
    !am=yo/sin(half*pi)**2   ! am = yo/1 = 0
 do i=is,ie
    !gf=(-1)**(m+1)*(xp(i,n)-xp(ii,n))/xo; yp(i,n)=am*sin(half*pi*gf)**2
    yp(i,n)=hinter(k1,k2,k3,k4,x0,x1,xp(i,n))
 end do
 do i=is,ie
    xp(i,n+1)=xp(i,n); yp(i,n+1)=yp(i,n)
 end do
 end do

!----- TOP & BOTTOM BOUNDARY POINTS IN XI-DIRECTION

    n=0
 if(rs==0) then
 do i=0,lxit
    xp(i,n)=xp(i,n+1)
 end do
 else
    sho=tla/litr; 
    am=2; xo=xa+ll*sho;
    ip=0; im=ll; call ogridf(xp(:,n),pxi,xa,xo,sho,am,0,lxit,im/2,im,ip)
    tmpa=xb-rs; tmpb=xp(lxie1,1)+re
    xjct=half*(xo+tmpa); err=1
 do while(abs(err)>sml)
    ip=ll; im=(lxi0-ll)/2; call ogridf(xp(:,n),pxi,xo,xjct,sho,am,0,lxit,im/2,im,ip)
    ip=ip+im; im=lxi0-ll-im; call ogridf(xp(:,n),pxi,xjct,tmpa,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/(2*shs1)-1; xjct=xjct+half*err*(2*shs1)
 end do
    xjct=half*(tmpa+tmpb); err=1
 do while(abs(err)>sml)
    ip=lxis1; im=lxi1/2; call ogridf(xp(:,n),pxi,tmpa,xjct,2*shs1,am,0,lxit,0,im,ip)
    ip=ip+im; im=lxi1-im; call ogridf(xp(:,n),pxi,xjct,tmpb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/she1-1; xjct=xjct+half*err*she1
 end do
    tmpa=tmpb; tmpb=xp(lxis3,1)-rs
    xjct=half*(tmpa+tmpb); err=1
 do while(abs(err)>sml)
    ip=lxis2; im=lxi2/2; call ogridf(xp(:,n),pxi,tmpa,xjct,she1,am,0,lxit,0,im,ip)
    ip=ip+im; im=lxi2-im; call ogridf(xp(:,n),pxi,xjct,tmpb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/(2*shs1)-1; xjct=xjct+half*err*she1
 end do
    tmpa=tmpb; tmpb=xp(lxie3,1)+re
    xjct=half*(tmpa+tmpb); err=1
 do while(abs(err)>sml)
    ip=lxis3; im=lxi3/2; call ogridf(xp(:,n),pxi,tmpa,xjct,2*shs1,am,0,lxit,0,im,ip)
    ip=ip+im; im=lxi3-im; call ogridf(xp(:,n),pxi,xjct,tmpb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/she1-1; xjct=xjct+half*err*she1
 end do
    ip=lxis4; im=lxi4-lxisz; call ogridf(xp(:,n),pxi,tmpb,xf,she1,am,0,lxit,im/2,im,ip)
    tmpa=xf; tmpb=xg
    xjct=(tmpa+tmpb)*half; err=1
 do while(abs(err)>sml)
    ip=lxis4+(lxi4-lxisz); im=lxisz/2; call ogridf(xp(:,n),pxi,tmpa,xjct,pxi(ip),am,0,lxit,im,im,ip)
    ip=ip+im; im=lxisz-im; call ogridf(xp(:,n),pxi,xjct,tmpb,pxi(ip),am,0,lxit,im/2,im,ip)
    err=pxi(ip+im)/sho-1; xjct=xjct+half*err*sho
 end do
 end if
 do i=0,lxit
    yp(i,n)=ya; xp(i,n+3)=xp(i,n); yp(i,n+3)=yd
 end do

!----- LEFT & RIGHT BOUNDARY POINTS IN ETA-DIRECTION

    ll=let0!4*tlb/(shs+sho)
 do m=0,5,5
 select case(m); case(0); xo=xa; case(5); xo=xf; end select
    am=2; yjct=half*yd; err=1
 do while(abs(err)>sml)
    jp=lets1; jm=ll/2; call ogridf(yq(:,m),qet,0.0_nr,yjct,shs2,am,0,lett,jm/2,jm,jp) !etamod will reduce the normal spacing
    jp=jp+jm; jm=ll-jm; call ogridf(yq(:,m),qet,yjct,yd,qet(jp),am,0,lett,jm/2,jm,jp)
    err=qet(jp+jm)/sho-1; yjct=yjct+half*err*sho
 end do
 do j=0,lete0
    yq(j,m)=-yq(lett-j,m)
 end do
 do j=0,lett
    xq(j,m)=xo
 end do
 end do
 if(k==0) then
    lq=ll
 end if

!----- VERTICAL INTERFACE POINTS IN ETA-DIRECTION

 if(rs==0) then
 do j=0,lett
    xq(j,1)=xb; xq(j,2)=xp(lxie1,1); yq(j,1)=yq(j,0); yq(j,2)=yq(j,3)
 end do
 else
 do n=1,2
     do m=1,4
     select case(n)
     case(1); tmpb=yd;
     tmpe=tla/litr; lets=lets1; !tmps=sqrt(half)*shs2; !etamod will reduce the normal spacing
       select case(m)
       case(1); tmpa=yp(lxis1,n); yjct=(tmpa+tmpb)*half;tmps=cos(pi4-delt1)*shs2; 
       case(2); tmpa=yp(lxie1,n); yjct=(tmpa+tmpb)*half;tmps=cos(pi4-delt1)*shs2;
       case(3); tmpa=yp(lxis3,n); yjct=(tmpa+tmpb)*half;tmps=cos(pi4-delt2)*shs2;
       case(4); tmpa=yp(lxie3,n); yjct=(tmpa+tmpb)*half;tmps=cos(pi4-delt2)*shs2;
       end select
     case(2); tmpa=ya; 
     tmps=tla/litr; lets=0; !tmpe=sqrt(half)*shs2; !etamod will reduce the normal spacing
       select case(m)
       case(1); tmpb=yp(lxis1,n); yjct=(tmpa+tmpb)*half;tmpe=cos(pi4+delt1)*shs2;
       case(2); tmpb=yp(lxie1,n); yjct=(tmpa+tmpb)*half;tmpe=cos(pi4+delt1)*shs2;
       case(3); tmpb=yp(lxis3,n); yjct=(tmpa+tmpb)*half;tmpe=cos(pi4+delt2)*shs2;
       case(4); tmpb=yp(lxie3,n); yjct=(tmpa+tmpb)*half;tmpe=cos(pi4+delt2)*shs2;
       end select
     end select
     am=2; err=1; 
     do while(abs(err)>sml)
        jp=lets; jm=ll/2; call ogridf(yq(:,m),qet,tmpa,yjct,tmps,am,0,lett,jm/2,jm,jp)
        jp=jp+jm; jm=ll-jm; call ogridf(yq(:,m),qet,yjct,tmpb,qet(jp),am,0,lett,jm/2,jm,jp)
        err=qet(jp+jm)/tmpe-1; yjct=yjct+half*err*tmpe
     end do
        !jp=jp+jm; jm=let0-ll; call ogridf(yq(:,m),qet,2*tlb,yd,qet(jp),am,0,lett,jm/2,jm,jp)
     !do j=0,lete0
     !   yq(j,m)=-yq(lett-j,m)
     !end do
     end do
 end do
 do n=1,2
 select case(n); case(1); js=0; je=lete0; nn=0; case(2); js=lets1; je=lett; nn=3; end select
 do m=1,4
 select case(m)
 case(1); i=lxis1; alph=delt1
    tmpa=sqrt((xp(i,n)-xp(i-2,n))**2+(yp(i,n)-yp(i-2,n))**2)
    tmpb=sqrt((xp(i,n)-xp(i+1,n))**2+(yp(i,n)-yp(i+1,n))**2)
    tmp=half*acos(half*(tmpa**2+tmpb**2-(xp(i+1,n)-xp(i-2,n))**2-(yp(i+1,n)-yp(i-2,n))**2)/(tmpa*tmpb))
 case(2); i=lxie1; alph=delt1
    tmpa=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)
    tmpb=sqrt((xp(i,n)-xp(i+2,n))**2+(yp(i,n)-yp(i+2,n))**2)
    tmp=half*acos(half*(tmpa**2+tmpb**2-(xp(i+2,n)-xp(i-1,n))**2-(yp(i+2,n)-yp(i-1,n))**2)/(tmpa*tmpb))
 case(3); i=lxis3; alph=delt2
    tmpa=sqrt((xp(i,n)-xp(i-2,n))**2+(yp(i,n)-yp(i-2,n))**2)
    tmpb=sqrt((xp(i,n)-xp(i+1,n))**2+(yp(i,n)-yp(i+1,n))**2)
    tmp=half*acos(half*(tmpa**2+tmpb**2-(xp(i+1,n)-xp(i-2,n))**2-(yp(i+1,n)-yp(i-2,n))**2)/(tmpa*tmpb))
 case(4); i=lxie3; alph=delt2
    tmpa=sqrt((xp(i,n)-xp(i-1,n))**2+(yp(i,n)-yp(i-1,n))**2)
    tmpb=sqrt((xp(i,n)-xp(i+2,n))**2+(yp(i,n)-yp(i+2,n))**2)
    tmp=half*acos(half*(tmpa**2+tmpb**2-(xp(i+2,n)-xp(i-1,n))**2-(yp(i+2,n)-yp(i-1,n))**2)/(tmpa*tmpb))
 end select
    sho=tan(tmp-(-1)**(m+n)*alph)
    tmpa=abs(xp(i,nn)-xp(i,n)); tmpb=abs(yp(i,nn)-yp(i,n))
    qq=(half*sho*tmpa)**2/abs(tmpb-sho*tmpa); pp=2*sqrt(qq)/sho
 do j=js,je
    xq(j,m)=(-1)**m*pp*(sqrt(abs(yq(j,m)-yp(i,n))+qq)-sqrt(qq))+xp(i,n)
 end do
 end do
 end do
 end if

!----- GRID OUTPUT

 if(rs==0) then
 do j=0,lett; do i=0,lxit
    xx(i,j)=half*(xp(i,1)+xp(i,2)); yy(i,j)=half*(yq(j,1)+yq(j,2)); zz(i,j)=zs(k)
 end do; end do
 else
 do n=0,2,2
 select case(n); case(0); js=0; je=lete0; case(2); js=lets1; je=lett; end select
 do m=0,4
 select case(m)
 case(0); is=0; ie=lxie0; nn=1; tmpa=0; tmpb=int(rs/abs(rs-sml))
 case(1); is=lxis1; ie=lxie1; nn=0; tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
 case(2); is=lxis2; ie=lxie2; nn=0; tmpa=int(re/abs(re-sml)); tmpb=int(rs/abs(rs-sml))
 case(3); is=lxis3; ie=lxie3; nn=0; tmpa=int(rs/abs(rs-sml)); tmpb=int(re/abs(re-sml))
 case(4); is=lxis4; ie=lxit; nn=1; tmpa=int(re/abs(re-sml)); tmpb=0
 end select
    gf=1
 do j=js,je; do i=is,ie
    pp=real(i-is,nr)/(ie-is); qq=real(j-js,nr)/(je-js)
    tmp=sin(half*pi*pp); ra0=((2-n)*real(je-j,nr)/(je-js)+n*qq)/2; ra1=gf+(2-gf)*ra0**2
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
 end if
 select case(mb); case(0,1,2,3,4); js=0; je=lete0; case(5,6,7,8,9); js=lets1; je=lett; end select
 select case(mb)
 case(0,5); is=0; ie=lxie0; case(1,6); is=lxis1; ie=lxie1; case(2,7); is=lxis2; ie=lxie2
 case(3,8); is=lxis3; ie=lxie3; case(4,9); is=lxis4; ie=lxit;
 end select
    np=nr*(je-js+1)*(ie-is+1)
    write(1,pos=np*k+1) ((xx(i,j),i=is,ie),j=js,je)
    write(1,pos=np*(k+lze0+1)+1) ((yy(i,j),i=is,ie),j=js,je)
    write(1,pos=np*(k+2*lze0+2)+1) ((zz(i,j),i=is,ie),j=js,je)
 end do

    close(1)

!----- GRID CHECKING

 if(ngridv==1) then
    open(1,file='misc/gridview'//cno(2)//cno(1)//cno(0)//'.dat')
    write(1,*) 'variables=v1,v2'; write(1,"('zone i=',i4,' j=',i4)") ie-is+1,je-js+1
 do j=js,je; do i=is,ie
    write(1,'(2es15.7)') xx(i,j),yy(i,j)
 end do; end do
    close(1)
 else
    open(1,file='misc/gridview'//cno(2)//cno(1)//cno(0)//'.dat'); close(1,status='delete')
 end if
 if(myid==0) then
    write(*,"('Grid generation is complete.')")
    write(*,"('Number of cells to resolve inflow turbulence in x:',i4)") lp
    write(*,"('Number of cells to resolve inflow turbulence in y:',i4)") lq
    write(*,"('Number of cells across sponge:',i4)") minloc(abs(xa+szth1-xp(:,1)),1)-1
 end if

 end if
    deallocate(xx,yy,zz,xp,yp,xq,yq,pxi,qet)

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

!===== FUNCTION FOR HORIZONTAL INTERFACE

 function hinter(k1,k2,k3,k4,x0,x1,x) result(y)

 real(nr) :: y
 real(nr), intent(in) :: x,x0,x1,k1,k2,k3,k4
 real(nr) :: a,b,c,d,l,invl,fx

 l = x1 - x0
 invl = 1.0_nr/l
 fx = (x - x0)*invl
 a = k1; b = l*k2
 d = l*k4 + 2*k3 + b + 2*a
 c = k3 - d - (a+b)
 y = a + b*fx + c*fx**2 + d*fx**3

 end function hinter

 !=====

 end module gridgen

!*****
