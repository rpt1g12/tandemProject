!*****
!***** BASIC SUBROUTINES
!*****

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
 real(nr),dimension(0:8),intent(in) :: albes,albee
 real(nr),intent(in) :: alpha,beta

 do i=is,ie
    xl(i,:)=1; xu(i,:)=1
 end do
    i=is
    xu(i,1)=1
    xu(i,2)=albes(0)
    xu(i,3)=albes(1)
    i=is+1
    xl(i,2)=albes(2)*xu(i-1,1)
    xu(i,1)=1/(1-xu(i-1,2)*xl(i,2))
    xu(i,2)=albes(3)-xu(i-1,3)*xl(i,2)
    xu(i,3)=albes(4)
    i=is+2
    xl(i,1)=albes(5)*xu(i-2,1)
    xl(i,2)=(albes(6)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=1/(1-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=albes(7)-xu(i-1,3)*xl(i,2)
    xu(i,3)=albes(8)
 do i=is+3,ie-3
    xl(i,1)=beta*xu(i-2,1)
    xl(i,2)=(alpha-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=1/(1-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=alpha-xu(i-1,3)*xl(i,2)
    xu(i,3)=beta
 end do
    i=ie-2
    xl(i,1)=albee(8)*xu(i-2,1)
    xl(i,2)=(albee(7)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=1/(1-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=albee(6)-xu(i-1,3)*xl(i,2)
    xu(i,3)=albee(5)
    i=ie-1
    xl(i,1)=albee(4)*xu(i-2,1)
    xl(i,2)=(albee(3)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=1/(1-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
    xu(i,2)=albee(2)-xu(i-1,3)*xl(i,2)
    i=ie
    xl(i,1)=albee(1)*xu(i-2,1)
    xl(i,2)=(albee(0)-xu(i-2,2)*xl(i,1))*xu(i-1,1)
    xu(i,1)=1/(1-xu(i-2,3)*xl(i,1)-xu(i-1,2)*xl(i,2))
 do i=is,ie
    xu(i,2:3)=xu(i,2:3)*xu(i,1)
 end do

 end subroutine penta

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS: ORIGINAL VERSION

 subroutine fcbco(fltk,fltkbc)
 
 real(nr),intent(in) :: fltk,fltkbc

    res=(fltk-fltkbc)/3; ra0=fltkbc; ra1=fltkbc+res; ra2=fltkbc+2*res

    call fcint(ra0,half)
    fctr=1/(1+5*alphf**3+alphf**2*(8+22*betf)+betf*(5+4*betf+60*betf**2)+5*alphf*(1+betf*(3+10*betf)))
    alphf01=(alphf*(1+alphf)*(1+4*alphf)+2*alphf*(7+3*alphf)*betf+24*(1-alphf)*betf**2-80*betf**3)*fctr
    betf02=(alphf**3+betf+14*alphf**2*betf+60*betf**3+alphf*betf*(3+46*betf))*fctr

    call fcint(ra1,half)
    fctr=1/(1+alphf**2*(9-5*betf)+alphf*(5+(35-29*betf)*betf)+betf*(5+6*(9-10*betf)*betf))
    alphf10=(9*alphf**3+10*betf**2*(8*betf-1)+5*alphf**2*(1+8*betf)+alphf*(1+betf*(4+81*betf)))*fctr
    alphf12=(alphf*(1+alphf*(5+9*alphf))+alphf*(5+36*alphf)*betf+(55*alphf-1)*betf**2+10*betf**3)*fctr
    betf13=betf*(1+alphf*(5+9*alphf)+5*betf+35*alphf*betf+50*betf**2)*fctr

    call fcint(ra2,half)
    betf20=betf; alphf21=alphf; alphf23=alphf; betf24=betf
    f20=fb+5*fc; f21=fa-10*fc; f23=fa-5*fc; f24=fb+fc; f25=fc
    
 end subroutine fcbco

!===== SUBROUTINE FOR BOUNDARY FILTER COEFFICIENTS: MODIFIED VERSION

 subroutine fcbcm(fltk)
 
 real(nr),intent(in) :: fltk

    res=(fltk-pi)/3; ra0=pi; ra1=pi+res; ra2=pi+2*res

    call fcint(ra0,half); alphf01=alphf; betf02=betf
    call fcint(ra1,half); alphf10=alphf; alphf12=alphf; betf13=betf
    call fcint(ra2,half); betf20=betf; alphf21=alphf; alphf23=alphf; betf24=betf
    
 end subroutine fcbcm

!===== SUBROUTINE FOR INTERIOR FILTER COEFFICIENTS

 subroutine fcint(fltk,fltr)
 
 real(nr),intent(in) :: fltk,fltr
 real(nr),dimension(3) :: cf

    cf(1)=cos(fltk); cf(2)=cos(2*fltk); cf(3)=cos(3*fltk)
    fctr=1/(30+5*(7-16*fltr)*cf(1)+2*(1+8*fltr)*cf(2)-3*cf(3))
    alphf=(20*(2*fltr-1)-30*cf(1)+12*(2*fltr-1)*cf(2)-2*cf(3))*fctr
    betf=(2*(13-8*fltr)+(33-48*fltr)*cf(1)+6*cf(2)-cf(3))*half*fctr
    fa=60*(1-fltr)*cos(half*fltk)**4*fctr; fb=-0.4_nr*fa; fc=fa/15

 end subroutine fcint

!===== SUBROUTINE FOR SUBDOMAIN-BOUNDARY COEFFICIENTS

 subroutine sbcco

 integer,dimension(1) :: imax
 integer,dimension(:),allocatable :: ipvt
 real(nr),dimension(:),allocatable :: temp
 real(nr),dimension(:,:),allocatable :: ax,bx,rx,sx

 do nt=0,1; lp=2*nt-1
 if(nt==0) then; ll=lmd; is=1; ie=2*(ll+1)
    allocate(ipvt(ie),temp(ie),ax(ie,ie),bx(ie,ie),rx(ie,ie),sx(ie,ie)); ax(:,:)=0; bx(:,:)=0
    ra0=sum(abc(:,0)); ra1=sum(abc(:,1)); ra2=sum(abc(:,2))
    ax(is,is:is+2)=(/one,alpha01,beta02/); bx(is,is:is+6)=(/-ra0,a01,a02,a03,a04,a05,a06/)
    ax(is+1,is:is+3)=(/alpha10,one,alpha12,beta13/); bx(is+1,is:is+6)=(/a10,-ra1,a12,a13,a14,a15,a16/)
    ax(is+2,is:is+4)=(/beta20,alpha21,one,alpha23,beta24/); bx(is+2,is:is+6)=(/a20,a21,-ra2,a23,a24,a25,a26/)
 do i=is+3,ie-3
    ax(i,i-2:i+2)=(/beta,alpha,one,alpha,beta/); bx(i,i-3:i+3)=(/-ac,-ab,-aa,zero,aa,ab,ac/)
 end do
 end if
 if(nt==1) then; ll=lmf; is=1; ie=2*(ll+1)
    allocate(ipvt(ie),temp(ie),ax(ie,ie),bx(ie,ie),rx(ie,ie),sx(ie,ie)); ax(:,:)=0; bx(:,:)=0
    ax(is,is:is+2)=(/one,alphf01,betf02/)
    ax(is+1,is:is+3)=(/alphf10,one,alphf12,betf13/)
    ax(is+2,is:is+4)=(/betf20,alphf21,one,alphf23,betf24/)
 if(nfbco==0) then
    bx(is+2,is:is+5)=(/f20,f21,-(f20+f21+f23+f24+f25),f23,f24,f25/)
 else
    m=nfbcn; res=fa+fb+fc; fctr=fex1+fex2+fex3; ra0=fa+2*fb+3*fc; ra1=fb+2*fc; ra2=fc
    bx(is,is:is+3)=(/ra0*fctr-res,fa,fb,fc/)
    bx(is,is+m:is+3*m:m)=bx(is,is+m:is+3*m:m)-ra0*(/fex1,fex2,fex3/)
    bx(is+1,is:is+4)=(/ra1*fctr+res,-2*res,fa,fb,fc/)
    bx(is+1,is+m:is+3*m:m)=bx(is+1,is+m:is+3*m:m)-ra1*(/fex1,fex2,fex3/)
    bx(is+2,is:is+5)=(/ra2*fctr+fb+fc,fa,-2*res,fa,fb,fc/)
    bx(is+2,is+m:is+3*m:m)=bx(is+2,is+m:is+3*m:m)-ra2*(/fex1,fex2,fex3/)
 end if
 do i=is+3,ie-3
    ax(i,i-2:i+2)=(/betf,alphf,one,alphf,betf/); bx(i,i-3:i+3)=(/fc,fb,fa,-2*(fa+fb+fc),fa,fb,fc/)
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

    ra0=pi/timf; ra1=ra0*min(timo,timf); ra2=ra0*min(timo+dtko,timf)

    fctr=1-cos(ra1)
    dfdt=ra0*sin(ra2)
    progmf=half*(fctr+dtk*dfdt)
    umf(:)=uoo(:)*((1-nsmf)*progmf+nsmf)

    fctr=sin(ra1)
    dfdt=ra0*cos(ra2)
    progmf=half*ra0*(fctr+dtk*dfdt)
    dudtmf(:)=uoo(:)*(1-nsmf)*progmf

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

 subroutine gridf(x,xxi,xo,xn,dxs,am,ns,lxi,mxic,mxin,ip)

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

 end subroutine gridf

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