!!*****
!!***** BASIC SUBROUTINES
!!*****

 module subroutineso

! use mainvar2d
 use mainvar3d
 implicit none

 contains

!===== SUBROUTINE FOR CHOLESKY DECOMPOSITION OF PENTADIAGONAL MATRICES

 subroutine penta(xu,xl,albes,albee,alpha,beta,is,ie)

 integer,intent(in) :: is,ie
 real(nr),dimension(0:lim,3),intent(inout) :: xu
 real(nr),dimension(0:lim,2),intent(inout) :: xl
 real(nr),dimension(-2:2,0:2),intent(in) :: albes,albee
 real(nr),intent(in) :: alpha,beta

 do i=is,ie
    xl(i,:)=one; xu(i,:)=one
 end do
    i=is
    xu(i,1)=one
    xu(i,2)=albes(1,0)
    xu(i,3)=albes(2,0)
    i=is+1
    xl(i,2)=albes(-1,1)*xu(i-1,1)
    xu(i,1)=one/(one-xu(i-1,2)*xl(i,2))
    xu(i,2)=albes(1,1)-xu(i-1,3)*xl(i,2)
    xu(i,3)=albes(2,1)
    i=is+2
    xl(i,1)=albes(-2,2)*xu(i-2,1)
    xl(i,2)=(albes(-1,2)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=albes(1,2)-xu(i-1,3)*xl(i,2)
    xu(i,3)=albes(2,2)
 do i=is+3,ie-3
    xl(i,1)=beta*xu(i-2,1)
    xl(i,2)=(alpha-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=alpha-xu(i-1,3)*xl(i,2)
    xu(i,3)=beta
 end do
    i=ie-2
    xl(i,1)=albee(2,2)*xu(i-2,1)
    xl(i,2)=(albee(1,2)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=albee(-1,2)-xu(i-1,3)*xl(i,2)
    xu(i,3)=albee(-2,2)
    i=ie-1
    xl(i,1)=albee(2,1)*xu(i-2,1)
    xl(i,2)=(albee(1,1)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=albee(-1,1)-xu(i-1,3)*xl(i,2)
    i=ie
    xl(i,1)=albee(2,0)*xu(i-2,1)
    xl(i,2)=(albee(1,0)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=one/(one-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
 do i=is,ie
    xu(i,2:3)=xu(i,2:3)*xu(i,1)
 end do

 end subroutine penta

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS: ORIGINAL VERSION

 subroutine fcbco(fltk,fltkbc,albef,fbc)
 
 real(nr),intent(in) :: fltk,fltkbc
 real(nr),dimension(-2:2,0:2),intent(inout) :: albef
 real(nr),dimension(0:4),intent(inout) :: fbc
 real(nr) :: alphf,betf,af,bf,cf

    res=(fltk-fltkbc)/3; ra0=fltkbc; ra1=fltkbc+res; ra2=fltkbc+2*res; albef(:,:)=0; albef(0,:)=1

    call fcint(ra0,half,alphf,betf,af,bf,cf)
    fctr=1/(1+5*alphf**3+alphf**2*(8+22*betf)+betf*(5+4*betf+60*betf**2)+5*alphf*(1+betf*(3+10*betf)))
    albef(1,0)=(alphf*(1+alphf)*(1+4*alphf)+2*alphf*(7+3*alphf)*betf+24*(1-alphf)*betf**2-80*betf**3)*fctr
    albef(2,0)=(alphf**3+betf+14*alphf**2*betf+60*betf**3+alphf*betf*(3+46*betf))*fctr

    call fcint(ra1,half,alphf,betf,af,bf,cf)
    fctr=1/(1+alphf**2*(9-5*betf)+alphf*(5+(35-29*betf)*betf)+betf*(5+6*(9-10*betf)*betf))
    albef(-1,1)=(9*alphf**3+10*betf**2*(8*betf-1)+5*alphf**2*(1+8*betf)+alphf*(1+betf*(4+81*betf)))*fctr
    albef(1,1)=(alphf*(1+alphf*(5+9*alphf))+alphf*(5+36*alphf)*betf+(55*alphf-1)*betf**2+10*betf**3)*fctr
    albef(2,1)=betf*(1+alphf*(5+9*alphf)+5*betf+35*alphf*betf+50*betf**2)*fctr

    call fcint(ra2,half,alphf,betf,af,bf,cf)
    albef(:,2)=(/betf,alphf,one,alphf,betf/); fbc(:)=(/bf+5*cf,af-10*cf,af-5*cf,bf+cf,cf/)
    
 end subroutine fcbco

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS: MODIFIED VERSION

 subroutine fcbcm(fltk,fltkbc,albef,fa,fb,fc)
 
 real(nr),intent(in) :: fltk,fltkbc
 real(nr),dimension(-2:2,0:2),intent(inout) :: albef
 real(nr),dimension(0:2),intent(inout) :: fa,fb,fc
 real(nr) :: alphf,betf,af,bf,cf

    res=(fltk-fltkbc)/3; ra0=fltkbc; ra1=fltkbc+res; ra2=fltkbc+2*res

    call fcint(ra0,half,alphf,betf,af,bf,cf)
    albef(:,0)=(/zero,zero,one,alphf,betf/); fa(0)=af; fb(0)=bf; fc(0)=cf
    call fcint(ra1,half,alphf,betf,af,bf,cf)
    albef(:,1)=(/zero,alphf,one,alphf,betf/); fa(1)=af; fb(1)=bf; fc(1)=cf
    call fcint(ra2,half,alphf,betf,af,bf,cf)
    albef(:,2)=(/betf,alphf,one,alphf,betf/); fa(2)=af; fb(2)=bf; fc(2)=cf

 end subroutine fcbcm

!===== SUBROUTINE FOR INTERIOR FILTER COEFFICIENTS

 subroutine fcint(fltk,fltr,alphf,betf,af,bf,cf)
 
 real(nr),intent(in) :: fltk,fltr
 real(nr),intent(inout) :: alphf,betf,af,bf,cf
 real(nr),dimension(3) :: cosf

    cosf(1)=cos(fltk); cosf(2)=cos(2*fltk); cosf(3)=cos(3*fltk)
    fctr=1/(30+5*(7-16*fltr)*cosf(1)+2*(1+8*fltr)*cosf(2)-3*cosf(3))
    alphf=(20*(2*fltr-1)-30*cosf(1)+12*(2*fltr-1)*cosf(2)-2*cosf(3))*fctr
    betf=(2*(13-8*fltr)+(33-48*fltr)*cosf(1)+6*cosf(2)-cosf(3))*half*fctr
    af=60*(1-fltr)*cos(half*fltk)**4*fctr; bf=-0.4_nr*af; cf=af/15

 end subroutine fcint

!===== SUBROUTINE FOR SUBDOMAIN-BOUNDARY COEFFICIENTS

 subroutine sbcco(np,nq)

 integer,intent(in) :: np,nq
 integer,dimension(1) :: imax
 integer,dimension(:),allocatable :: ipvt
 real(nr),dimension(:),allocatable :: temp
 real(nr),dimension(:,:),allocatable :: ax,bx,rx,sx

 do nt=0,1; lp=2*nt-1
 if(nt==0) then; ll=lmd; is=1; ie=2*(ll+1)
    allocate(ipvt(ie),temp(ie),ax(ie,ie),bx(ie,ie),rx(ie,ie),sx(ie,ie)); ax(:,:)=0; bx(:,:)=0
    ax(is,is:is+2)=albed(0:2,0,np); ax(is+1,is:is+3)=albed(-1:2,1,np); ax(is+2,is:is+4)=albed(-2:2,2,np)
    ra0=sum(abc(:,0)); bx(is,is:is+6)=(/-ra0,a01,a02,a03,a04,a05,a06/)
    ra1=sum(abc(:,1)); bx(is+1,is:is+6)=(/a10,-ra1,a12,a13,a14,a15,a16/)
    ra2=sum(abc(:,2)); bx(is+2,is:is+6)=(/a20,a21,-ra2,a23,a24,a25,a26/)
 do i=is+3,ie-3
    ax(i,i-2:i+2)=(/beta,alpha,one,alpha,beta/); bx(i,i-3:i+3)=(/-ac,-ab,-aa,zero,aa,ab,ac/)
 end do
 end if
 if(nt==1) then; ll=lmf; is=1; ie=2*(ll+1)
    allocate(ipvt(ie),temp(ie),ax(ie,ie),bx(ie,ie),rx(ie,ie),sx(ie,ie)); ax(:,:)=0; bx(:,:)=0
    ax(is,is:is+2)=albef(0:2,0,nq); ax(is+1,is:is+3)=albef(-1:2,1,nq); ax(is+2,is:is+4)=albef(-2:2,2,nq)
 if(nq==0) then
    ra2=sum(fbc(:)); bx(is+2,is:is+5)=(/fbc(0),fbc(1),-ra2,fbc(2),fbc(3),fbc(4)/)
 else
    bx(is,is:is+3)=(/-2*(fam(0)+fbm(0)+fcm(0)),fam(0),fbm(0),fcm(0)/)
    bx(is+1,is:is+4)=(/fam(1),-2*(fam(1)+fbm(1)+fcm(1)),fam(1),fbm(1),fcm(1)/)
    bx(is+2,is:is+5)=(/fbm(2),fam(2),-2*(fam(2)+fbm(2)+fcm(2)),fam(2),fbm(2),fcm(2)/)
    res=fex(1)+fex(2)+fex(3)
    bx(is,is)=bx(is,is)+fam(0)*(res+1)+fbm(0)*(2*res+1)+fcm(0)*(3*res+1)
    bx(is+1,is)=bx(is+1,is)+fbm(1)*(res+1)+fcm(1)*(2*res+1)
    bx(is+2,is)=bx(is+2,is)+fcm(2)*(res+1)
 do ip=1,3; ii=is+ip*mfbi
    bx(is,ii)=bx(is,ii)-fex(ip)*(fam(0)+2*fbm(0)+3*fcm(0))
    bx(is+1,ii)=bx(is+1,ii)-fex(ip)*(fbm(1)+2*fcm(1))
    bx(is+2,ii)=bx(is+2,ii)-fex(ip)*fcm(2)
 end do
 end if
 do i=is+3,ie-3
    ax(i,i-2:i+2)=(/betf,alphf,one,alphf,betf/)
    bx(i,i-3:i+3)=(/fc,fb,fa,-2*(fa+fb+fc),fa,fb,fc/)
 end do
 end if
    ax(ie-2,ie:is:-1)=ax(is+2,:); bx(ie-2,ie:is:-1)=lp*bx(is+2,:)
    ax(ie-1,ie:is:-1)=ax(is+1,:); bx(ie-1,ie:is:-1)=lp*bx(is+1,:)
    ax(ie,ie:is:-1)=ax(is,:); bx(ie,ie:is:-1)=lp*bx(is,:)

    rx(:,:)=ax(:,:); ipvt(:)=(/(i,i=is,ie)/)
 do i=is,ie
    imax(:)=maxloc(abs(rx(i:ie,i))); m=i-1+imax(1)
 if(m/=i) then
    ipvt((/m,i/))=ipvt((/i,m/)); rx((/m,i/),:)=rx((/i,m/),:)
 end if
    ra0=1/rx(i,i); temp(:)=rx(:,i)
 do j=is,ie
    ra1=ra0*rx(i,j); rx(:,j)=rx(:,j)-ra1*temp(:); rx(i,j)=ra1
 end do
    rx(:,i)=-ra0*temp(:); rx(i,i)=ra0
 end do
    sx(:,ipvt(:))=rx(:,:); rx(:,:)=ax(:,:)
    i=ie/2-1; rx(i,i+2)=0; rx(i+1,i+2)=0; rx(i+1,i+3)=0
    i=ie/2+2; rx(i,i-2)=0; rx(i-1,i-2)=0; rx(i-1,i-3)=0
    ax(:,:)=matmul(rx(:,:),matmul(sx(:,:),bx(:,:)))
    i=ie/2+1; pbco(ll:0:-1,0,nt)=ax(i,is:is+ll); pbci(0:ll,0,nt)=ax(i,is+ll+1:ie)
    i=ie/2+2; pbco(ll:0:-1,1,nt)=ax(i,is:is+ll); pbci(0:ll,1,nt)=ax(i,is+ll+1:ie)
    deallocate(ipvt,temp,ax,bx,rx,sx)
 end do

 end subroutine sbcco

!===== SUBROUTINE FOR MOVING FRAME VELOCITIES

 subroutine movef(dtko,dtk)

 real(nr),intent(in) :: dtko,dtk

 if(nsmf==0) then
    ra0=pi/timf; ra1=ra0*min(timo,timf); ra2=ra0*min(timo+dtko,timf)

    fctr=1-cos(ra1)
    dfdt=ra0*sin(ra2)
    progmf=half*(fctr+dtk*dfdt)
    umf(:)=progmf*uoo(:)

    fctr=sin(ra1)
    dfdt=ra0*cos(ra2)
    progmf=half*ra0*(fctr+dtk*dfdt)
    dudtmf(:)=progmf*uoo(:)
 else
    umf(:)=uoo(:); dudtmf(:)=0
 end if

 end subroutine movef

!===== SUBROUTINE FOR BLASIUS LAMINAR BOUNDARY LAYER

 subroutine lambl(x,y,blu,blv,blm,lbl)

 integer,intent(in) :: lbl
 real(nr),dimension(0:lbl),intent(in) :: x,y
 real(nr),dimension(0:lbl),intent(inout) :: blu,blv,blm

 real(nr) :: blas,spr,eta,etb

    blas=0.3320573362151963_nr; spr=sqrt(prndtl)
    ra0=half*pi; ra1=3*blas*blas/560; ra2=11*blas/420
 do i=0,lbl
    fctr=1/sqrt(x(i))
    eta=fctr*sqrtrema*y(i); etb=spr*eta
    res=0.000001_nr*eta**4*exp(eta**1.3625_nr)
    blu(i)=(blas*eta+ra1*eta**4+res)/(1+ra2*eta**3+res)
    blv(i)=0.8604_nr*fctr*sqrtremai*(sin(ra0*tanh(exp(0.177_nr*eta)-1)))**1.96_nr
    res=0.000001_nr*etb**4*exp(etb**1.3625_nr)
    blm(i)=(blas*etb+ra1*etb**4+res)/(1+ra2*etb**3+res)
 end do

 end subroutine lambl

!===== SUBROUTINE FOR GRID LINE GENERATION

 subroutine gridf(x,xxi,xo,xn,dxo,dxn,lxi,mxin,ip)

 integer,intent(in) :: lxi,mxin,ip
 real(nr),dimension(0:lxi),intent(inout) :: x,xxi
 real(nr),intent(in) :: xo,xn,dxo,dxn

 integer :: i,ii
 real(nr) :: dxoo,dxnn,aa,bb,cc,ee,dd,xi,fctr

    dxoo=dxo; dxnn=dxn
 if(dxo==sml) then
    dxoo=2*(xn-xo)/mxin-dxn
 end if
 if(dxn==sml) then
    dxnn=2*(xn-xo)/mxin-dxo
 end if
    aa=6*(xn-xo)-3*mxin*(dxoo+dxnn)
    bb=15*(xo-xn)+mxin*(8*dxoo+7*dxnn)
    cc=10*(xn-xo)-mxin*(6*dxoo+4*dxnn)
    dd=mxin*dxoo; fctr=one/mxin
 do i=0,mxin
    ii=i+ip; xi=i*fctr
    x(ii)=aa*xi**5+bb*xi**4+cc*xi**3+dd*xi+xo
    xxi(ii)=fctr*(5*aa*xi**4+4*bb*xi**3+3*cc*xi**2+dd)
 end do

 end subroutine gridf

!===== SUBROUTINE FOR GRID LINE GENERATION : OLD VERSION

 subroutine ogridf(x,xxi,xo,xn,dxs,am,ns,lxi,mxic,mxin,ip)

 integer,intent(in) :: ns,lxi,mxic,mxin,ip
 real(nr),dimension(0:lxi),intent(inout) :: x,xxi
 real(nr),intent(in) :: xo,xn,dxs,am

 integer :: i,ii
 real(nr) :: amp0,amm0,alp0,alp1,s0,s1,c0,c1,aa,bb,xi,xic,xin,xii

    amp0=am+1; amm0=am-1
    xic=mxic; xin=mxin
    alp0=xin+amm0*xic; alp1=am*xin-amm0*xic
 if(ns==0) then
    s0=dxs; s1=(amp0*(xn-xo)-alp0*s0)/alp1
 else
    s1=dxs; s0=(amp0*(xn-xo)-alp1*s1)/alp0
 end if
    c0=xo; c1=xn-s1*xin
    aa=(s1-s0)*xic/xin; bb=(s1-s0)*(xin-xic)/xin
 if(mxic/=0) then
 do i=0,mxic
    ii=i+ip; xi=i; xii=xi/xic
    x(ii)=s0*xi+aa*xic*xii**amp0/amp0+c0; xxi(ii)=s0+aa*xii**am
 end do
 end if
 if(mxic/=mxin) then
 do i=mxic,mxin
    ii=i+ip; xi=i; xii=(xin-xi)/(xin-xic)
    x(ii)=s1*xi+bb*(xin-xic)*xii**amp0/amp0+c1; xxi(ii)=s1-bb*xii**am
 end do
 end if

 end subroutine ogridf

!===== SUBROUTINE FOR CHARACTER STRING CONVERSION

 subroutine strio(nfile,lh,cinput)

 integer,intent(in) :: nfile
 integer,intent(inout) :: lh
 character(16),intent(in) :: cinput

 do ll=1,len_trim(cinput)
    write(nfile,pos=4*lh+1) ichar(cinput(ll:ll)); lh=lh+1
 end do
    write(nfile,pos=4*lh+1) 0; lh=lh+1

 end subroutine strio

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 3D

 function indx3(i,j,k,nn) result(lm)

 integer,intent(in) :: i,j,k,nn
 integer :: lm

 select case(nn)
 case(1); lm=(k*(let+1)+j)*(lxi+1)+i
 case(2); lm=(j*(let+1)+i)*(lxi+1)+k
 case(3); lm=(i*(let+1)+k)*(lxi+1)+j
 end select

 end function indx3

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 2D

 function indx2(i,j,nn) result(lm)

 integer,intent(in) :: i,j,nn
 integer :: lm

 select case(nn)
 case(1); lm=j*(lxi+1)+i
 case(2); lm=i*(lxi+1)+j
 end select

 end function indx2

!=====

 end module subroutineso

!*****
