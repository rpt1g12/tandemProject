!***************************
!***** RPT INTERPOLATION MODULE *****
!***************************

module  rptinter

use mainvar3d
use subroutineso
use subroutines3d
use problemcase
use mpi
use rpt
implicit none
 real(nr),dimension(:,:),allocatable :: xyz2,ixis
 integer :: lxii,leti,lzei
 integer :: lxiio,letio,lzeio
 integer :: l2,ncom,ilze0,olze0

real(nr),dimension(:),allocatable :: xxi,xet,xze
real(nr),dimension(:),allocatable :: yxi,yet,yze
real(nr),dimension(:),allocatable :: zxi,zet,zze
real(nr),dimension(:),allocatable :: f
real(nr),dimension(:),allocatable :: fxi,fet,fze
real(nr),dimension(:),allocatable :: fetxi
real(nr),dimension(:),allocatable :: fzexi,fzeet,fzeetxi
real(nr),dimension(:),allocatable :: lvarr,lvarr2

real(nr), dimension(0:1,1:3) :: bounds
real(nr), dimension(1:5) :: outside

integer, dimension(:), allocatable :: ilximb,iletmb,ilzemb
integer, dimension(:), allocatable :: ilxibk,iletbk,ilzebk
integer, dimension(:), allocatable :: olximb,oletmb,olzemb
integer, dimension(:), allocatable :: olxibk,oletbk,olzebk
logical :: iflag,igflag
 character(20) :: icrestart
 integer(ni) :: niter,ilmx
contains

subroutine interSetUp()
    integer :: i1,i2

    open(9,file='inputi.dat')
    read(9,*) cinput !Flags
    read(9,*) i1,i2
    read(9,*) cinput !Blocks
    read(9,*) mbk,bkx,bky,bkz

    iflag=(i1==1);igflag=(i2==1)
    allocate(ilxibk(0:bkx-1),iletbk(0:bky-1),ilzebk(0:bkz-1))
    allocate(olxibk(0:bkx-1),oletbk(0:bky-1),olzebk(0:bkz-1))

    read(9,*) cinput !Input
    read(9,*) cinput,ilxibk(0:bkx-1) ! # points in xi per columns
    read(9,*) cinput,iletbk(0:bky-1) ! # points in eta per row
    read(9,*) cinput,ilzebk(0:bkz-1) ! # points in zeta per plane
    read(9,*) cinput !Output
    read(9,*) cinput,olxibk(0:bkx-1) ! # points in xi per columns
    read(9,*) cinput,oletbk(0:bky-1) ! # points in eta per row
    read(9,*) cinput,olzebk(0:bkz-1) ! # points in zeta per plane
    close(9)


end subroutine interSetUp
!====================================================================================
!=====PROBLEM SETUP
!====================================================================================
 subroutine setup(io)
 integer, intent(in) :: io
!===== INPUT PARAMETERS

    open(9,file='inputo.dat')
    read(9,*) cinput!,mbk,bkx,bky,bkz
    read(9,*) cinput,nts,nto
    read(9,*) cinput,nscrn,nsgnl
    read(9,*) cinput,ndata
    read(9,*) cinput,nkrk
    read(9,*) cinput,nviscous
    read(9,*) cinput,nsmf
    read(9,*) cinput,nfskp
    read(9,*) cinput,nrestart
    read(9,*) cinput,nextrabc,nextgcic
    read(9,*) cinput,reoo,tempoo
    read(9,*) cinput,amach1,amach2,amach3
    read(9,*) cinput,wtemp
    read(9,*) cinput,cfl
    read(9,*) cinput,tmax,timf,tsam
    read(9,*) cinput,fltk,fltkbc
    read(9,*) cinput,dto
    read(9,*) cinput,forcing,amfor
    read(9,*) cinput,aoa
    read(9,*) cinput,LES,smago1,smago2
    read(9,*) cinput,output,ogrid,osol,oblock
    close(9)

    cinput=cinput; fltk=pi*fltk; fltkbc=pi*fltkbc
    rhooo=one; poo=one/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1*amach1+amach2*amach2+amach3*amach3)
    srefoo=111/tempoo; srefp1dre=(srefoo+one)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=one/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo

    abc(:,0)=(/a01,a02,a03,a04,a05,a06/)
    abc(:,1)=(/a10,a12,a13,a14,a15,a16/)
    abc(:,2)=(/a20,a21,a23,a24,a25,a26/)

    ll=3+5*(ndata+1)
    !allocate(times(0:ndata),cfilet(-1:ndata),ctecplt(-1:ndata),varm(0:1,0:mpro),varmin(ll),varmax(ll))
    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),npc(0:mbk,3))
    if (.not.allocated(lxibk)) allocate(lxibk(0:bkx-1),letbk(0:bky-1),lzebk(0:bkz-1))

    call inputext

    selectcase(io)
    case(0); 
       lxibk(:)=ilxibk(:)
       letbk(:)=iletbk(:)
       lzebk(:)=ilzebk(:)
    case(1);
       lxibk(:)=olxibk(:)
       letbk(:)=oletbk(:)
       lzebk(:)=olzebk(:)
    end select

    ! Use input data to fill old arrays
    do k = 0, bkz-1
       do j = 0, bky-1
          do i = 0, bkx-1; l=k*(bkx*bky)+j*bkx+i
             lximb(l)=lxibk(i)
             letmb(l)=letbk(j)
             lzemb(l)=lzebk(k)
             
             npc(l,1)=1
             npc(l,2)=1
             npc(l,3)=1
          end do
       end do
    end do
    lze0 = lzebk(0)

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

    mo(0)=0
 do mm=1,mbk
    mo(mm)=mo(mm-1)+npc(mm-1,1)*npc(mm-1,2)*npc(mm-1,3)
 end do
 do mm=0,mbk
    if(myid>=mo(mm)) mb=mm 
 end do
    lxio=lximb(mb); leto=letmb(mb); lzeo=lzemb(mb)
     ! rpt- Create communicator per block
     CALL MPI_COMM_SPLIT(icom,mb,myid,bcom,ierr)   

    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10)
    cno=achar(no+48); cnzone=cno(2)//cno(1)//cno(0)
    czone='zone'//cnzone;
    cgrid='misc/grid'//cnzone//'.dat';
    crestart='rsta/restart'//cnzone//'.dat'
    icrestart='irsta/restart'//cnzone//'.dat'

    no(4)=myid/10000; no(3)=mod(myid,10000)/1000;
    no(2)=mod(myid,1000)/100; no(1)=mod(myid,100)/10; no(0)=mod(myid,10)
    cno=achar(no+48); cnnode=cno(4)//cno(3)//cno(2)//cno(1)//cno(0)

    call domdcomp

    ip=mod(myid-mo(mb),npc(mb,1))
    jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

    ! rpt- Store processors coordinates
    mpc=(/ip,jp,kp/)

    ncds(1)=mo(ms(1))+kp*npc(ms(1),2)*npc(ms(1),1)+jp*npc(ms(1),1)+npc(ms(1),1)-1
    ncde(1)=mo(me(1))+kp*npc(me(1),2)*npc(me(1),1)+jp*npc(me(1),1)

    ncds(2)=mo(ms(2))+kp*npc(ms(2),2)*npc(ms(2),1)+(npc(ms(2),2)-1)*npc(ms(2),1)+ip
    ncde(2)=mo(me(2))+kp*npc(me(2),2)*npc(me(2),1)+ip

    ncds(3)=mo(ms(3))+(npc(ms(3),3)-1)*npc(ms(3),2)*npc(ms(3),1)+jp*npc(ms(3),1)+ip
    ncde(3)=mo(me(3))+jp*npc(me(3),1)+ip

 do nn=1,3
 select case(nn)
 case (1); ll=lxio; lp=ip; mp=1
 case (2); ll=leto; lp=jp; mp=npc(mb,1)
 case (3); ll=lzeo; lp=kp; mp=npc(mb,1)*npc(mb,2)
 end select
    ma=npc(mb,nn)
 if(ma==1) then
    l=ll;
    nbc(0,nn)=nbcs(nn); nbc(1,nn)=nbce(nn);
    ncd(0,nn)=ncds(nn); ncd(1,nn)=ncde(nn)
 end if
 if(ma>=2) then
 if(lp==0) then
       l=ll-((ll+1)/ma)*(ma-1);
       nbc(0,nn)=nbcs(nn); nbc(1,nn)=40;
       ncd(0,nn)=ncds(nn); ncd(1,nn)=myid+mp
 end if
 if(lp>0.and.lp<ma-1) then
       l=(ll+1)/ma-1;
       nbc(0,nn)=40; nbc(1,nn)=40;
       ncd(0,nn)=myid-mp; ncd(1,nn)=myid+mp
 end if
 if(lp==ma-1) then
       l=(ll+1)/ma-1;
       nbc(0,nn)=40; nbc(1,nn)=nbce(nn);
       ncd(0,nn)=myid-mp; ncd(1,nn)=ncde(nn)
 end if
 end if
 select case(nn); case (1); lxi=l; case (2); let=l; case (3); lze=l; end select
 end do

!===== SUBDOMAIN SIZES & WRITING START POSITIONS IN OUTPUT FILE

 if(myid==0) then
    lxim(0)=lxi; letm(0)=let; lzem(0)=lze
 do mp=1,mpro
    itag=1; call MPI_RECV(lxim(mp),1,MPI_INTEGER4,mp,itag,icom,ista,ierr)
    itag=2; call MPI_RECV(letm(mp),1,MPI_INTEGER4,mp,itag,icom,ista,ierr)
    itag=3; call MPI_RECV(lzem(mp),1,MPI_INTEGER4,mp,itag,icom,ista,ierr)
 end do
 else; itag=myid
    itag=1; call MPI_SEND(lxi,1,MPI_INTEGER4,0,itag,icom,ierr)
    itag=2; call MPI_SEND(let,1,MPI_INTEGER4,0,itag,icom,ierr)
    itag=3; call MPI_SEND(lze,1,MPI_INTEGER4,0,itag,icom,ierr)
 end if
    call MPI_BCAST(lxim(:),npro,MPI_INTEGER4,0,icom,ierr)
    call MPI_BCAST(letm(:),npro,MPI_INTEGER4,0,icom,ierr)
    call MPI_BCAST(lzem(:),npro,MPI_INTEGER4,0,icom,ierr)

    ltomb=(lxio+1)*(leto+1)*(lzeo+1)

    lmx=(lxi+1)*(let+1)*(lze+1)-1
    lim=(lxi+1)+(let+1)+(lze+1)-1

    ijk(1,1)=lxi; ijk(2,1)=let; ijk(3,1)=lze
    ijk(1,2)=let; ijk(2,2)=lze; ijk(3,2)=lxi
    ijk(1,3)=lze; ijk(2,3)=lxi; ijk(3,3)=let

    nbsize(:)=(ijk(2,:)+1)*(ijk(3,:)+1)

 do mm=0,mbk
    lpos(mo(mm))=0
 do i=1,npc(mm,1)-1
    mp=mo(mm)+i
    lpos(mp)=lpos(mp-1)+lxim(mp-1)+1
 end do
    jp=npc(mm,1)
       do j=1,npc(mm,2)-1;
          do i=0,npc(mm,1)-1
    mp=mo(mm)+j*jp+i
    lpos(mp)=lpos(mp-jp)+(lximb(mm)+1)*(letm(mp-jp)+1)
           end do;
       end do
    kp=npc(mm,1)*npc(mm,2)
       do k=1,npc(mm,3)-1;
          do j=0,npc(mm,2)-1;
             do i=0,npc(mm,1)-1
    mp=mo(mm)+k*kp+j*jp+i
    lpos(mp)=lpos(mp-kp)+(lximb(mm)+1)*(letmb(mm)+1)*(lzem(mp-kp)+1)
             end do;
          end do;
 end do
 end do

    ! rpt- Find start indices depending on proc coordinates
    allocate(ibegin(0:npc(mb,1)))
    allocate(jbegin(0:npc(mb,2)))
    allocate(kbegin(0:npc(mb,3)))
    ! setup first process
    ibegin(0)=0
    jbegin(0)=0
    kbegin(0)=0
    ! setup i-start indices
    do i=1,npc(mb,1)
    mp=mo(mb)+i
    ibegin(i)=ibegin(i-1)+lxim(mp-1)+1
    end do
    ! setup j-start indices
    do j=1,npc(mb,2)
    mp=mo(mb)+j*npc(mb,1)
    jbegin(j)=jbegin(j-1)+letm(mp-1)+1
    end do
    ! setup k-start indices
    do k=1,npc(mb,3)
    mp=mo(mb)+k*npc(mb,1)*npc(mb,2)
    kbegin(k)=kbegin(k-1)+lzem(mp-1)+1
    end do

    ! rpt- #Points in block per direction
    mbijkl=(/lxio,leto,lzeo/)+1
    ! rpt- #Points in proccessor per direction
    mpijkl=(/lxi,let,lze/)+1
    ! rpt- Starts in proccessor per direction
    mpijks=(/ibegin(mpc(1)),jbegin(mpc(2)),kbegin(mpc(3))/)
    ! rpt- Ends in proccessor per direction
    mpijke=mpijks+(/lxi,let,lze/)

 end subroutine setup

!====================================================================================
!=====PREPARE ARRAYS
!====================================================================================
 subroutine deallocateArrays

 if(allocated(iit)) deallocate(iit,idsgnl,lsgnl)
 if(allocated(qo)) deallocate(qo,qa,de)
 if(allocated(xim)) deallocate(xim,etm,zem,rr,ss)
 if(allocated(p)) deallocate(p,yaco,varr)
 if(allocated(lximb)) deallocate(lximb,letmb,lzemb,lhmb,mo,npc)
 if(allocated(mxc)) deallocate(mxc,ran,sit,ait,xit,yit,zit)
 if(allocated(drva1)) deallocate(drva1,drva2,drva3)
 if(allocated(drvb1)) deallocate(drvb1,drvb2,drvb3)
 if(allocated(send1)) deallocate(send1,send2,send3)
 if(allocated(recv1)) deallocate(recv1,recv2,recv3)
 if(allocated(cm1)) deallocate(cm1,cm2,cm3)
 if(allocated(xu)) deallocate(xu,yu,xl,yl,li,sa,sb)
 if(allocated(lio)) deallocate (lio)
 if(allocated(ibegin)) deallocate(ibegin,jbegin,kbegin)
 if(allocated(rpex))deallocate(rpex,sbcc)
 if(allocated(q8)) deallocate(q8)
 if(allocated(xyz4)) deallocate(xyz4)
 if(allocated(fout)) deallocate(fout)
    
 end subroutine deallocateArrays

!====================================================================================
!=====PREPARE ARRAYS
!====================================================================================
 subroutine prepareArrays
!===== ALLOCATION OF MAIN ARRAYS

    allocate(qo(0:lmx,5),qa(0:lmx,5),de(0:lmx,5))
    allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
    allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx))

    ii=nbsize(1)-1; jj=nbsize(2)-1; kk=nbsize(3)-1
    allocate(drva1(0:ii,5,0:1),drva2(0:jj,5,0:1),drva3(0:kk,5,0:1))
    allocate(drvb1(0:ii,5,0:1),drvb2(0:jj,5,0:1),drvb3(0:kk,5,0:1))
    allocate(send1(0:ii,0:2,0:1),send2(0:jj,0:2,0:1),send3(0:kk,0:2,0:1))
    allocate(recv1(0:ii,0:2,0:1),recv2(0:jj,0:2,0:1),recv3(0:kk,0:2,0:1))
    allocate(cm1(0:ii,3,0:1),cm2(0:jj,3,0:1),cm3(0:kk,3,0:1))

    allocate(xu(0:lim,3),yu(0:lim,3),xl(0:lim,2),yl(0:lim,2),li(0:lim),sa(0:lim),sb(0:lim))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

    albed(-2:2,0,0)=(/zero,zero,one,alpha01,beta02/)
    albed(-2:2,1,0)=(/zero,alpha10,one,alpha12,beta13/)
    albed(-2:2,2,0)=(/beta20,alpha21,one,alpha23,beta24/)

    albed(-2:2,0,1)=(/zero,zero,one,alpha,beta/)
    albed(-2:2,1,1)=(/zero,alpha,one,alpha,beta/)
    albed(-2:2,2,1)=(/beta,alpha,one,alpha,beta/)

    call fcbcm(fltk,fltkbc,albef(:,:,0),fam(:),fbm(:),fcm(:))
    call fcint(fltk,half,alphf,betf,fa,fb,fc)
    albef(-2:2,0,1)=(/zero,zero,one,alphf,betf/)
    albef(-2:2,1,1)=(/zero,alphf,one,alphf,betf/)
    albef(-2:2,2,1)=(/betf,alphf,one,alphf,betf/)

    pbco(:,:,:)=0; pbci(:,:,:)=0; call sbcco
 do nt=0,1; do j=0,1; ii=lmd+nt*(lmf-lmd)
    pbcot(j,nt)=sum(pbco(0:ii,j,nt))
 end do; end do

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

 do nn=1,3
 select case(nn)
       case(1); is=0; ie=is+lxi;
       case(2); is=lxi+1; ie=is+let;
       case(3); is=lxi+let+2; ie=is+lze
 end select
 do ip=0,1; np=nbc(ip,nn)
 select case(np)
 case(10,20,25,30); ndf(ip,0,nn)=0; ndf(ip,1,nn)=0
 case(35,40,45); ndf(ip,0,nn)=1; ndf(ip,1,nn)=1
 end select
 end do
    ns=ndf(0,0,nn); ne=ndf(1,0,nn)
    call penta(xu(:,:),xl(:,:),albed(:,:,ns),albed(:,:,ne),alpha,beta,is,ie)
    ns=ndf(0,1,nn); ne=ndf(1,1,nn)
    call penta(yu(:,:),yl(:,:),albef(:,:,ns),albef(:,:,ne),alphf,betf,is,ie)
 end do

    allocate(lio(0:let,0:lze))
 do k=0,lze; kp=k*(leto+1)*(lxio+1)
 do j=0,let; jp=j*(lxio+1)
    lio(j,k)=jp+kp
 end do
 end do

  end subroutine prepareArrays

!====================================================================================
!=====CREATE GRID
!====================================================================================
  subroutine getGrid

    call makegrid
    call MPI_BARRIER(icom,ierr)

    call rdGrid
    allocate(xyz4(0:lmx,3))
    xyz4(:,:)=ss(:,:)

    call MPI_BARRIER(icom,ierr)
 if(myid==mo(mb)) then
    open(9,file=cgrid); close(9,status='delete')
 end if
  end subroutine getGrid


!====================================================================================
!=====COMPUTE METRIC TERMS
!====================================================================================
 subroutine getMetrics
 !===== COMPUTE INVERSE METRICS
     rr(:,1)=ss(:,1)
    m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)
 
     rr(:,1)=ss(:,2)
    m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)
 
     rr(:,1)=ss(:,3)
    m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     de(:,1)=rr(:,1); de(:,2)=rr(:,2); de(:,3)=rr(:,3)
 
     allocate(xxi(0:ltomb-1),xet(0:ltomb-1),xze(0:ltomb-1))
     allocate(yxi(0:ltomb-1),yet(0:ltomb-1),yze(0:ltomb-1))
     allocate(zxi(0:ltomb-1),zet(0:ltomb-1),zze(0:ltomb-1))

 !===== COMPUTE METRICS
     xim(:,1)=qa(:,2)*de(:,3)-de(:,2)*qa(:,3)
     xim(:,2)=de(:,2)*qo(:,3)-qo(:,2)*de(:,3)
     xim(:,3)=qo(:,2)*qa(:,3)-qa(:,2)*qo(:,3)
     etm(:,1)=qa(:,3)*de(:,1)-de(:,3)*qa(:,1)
     etm(:,2)=de(:,3)*qo(:,1)-qo(:,3)*de(:,1)
     etm(:,3)=qo(:,3)*qa(:,1)-qa(:,3)*qo(:,1)
     zem(:,1)=qa(:,1)*de(:,2)-de(:,1)*qa(:,2)
     zem(:,2)=de(:,1)*qo(:,2)-qo(:,1)*de(:,2)
     zem(:,3)=qo(:,1)*qa(:,2)-qa(:,1)*qo(:,2)
    

    yaco(:)=three/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
               +qa(:,1)*xim(:,2)+qa(:,2)*etm(:,2)+qa(:,3)*zem(:,2)&
               +de(:,1)*xim(:,3)+de(:,2)*etm(:,3)+de(:,3)*zem(:,3))

     varr=xim(:,1)*yaco(:); call joinBlock; xxi=lvarr
     varr=etm(:,1)*yaco(:); call joinBlock; xet=lvarr
     varr=zem(:,1)*yaco(:); call joinBlock; xze=lvarr
     varr=xim(:,2)*yaco(:); call joinBlock; yxi=lvarr
     varr=etm(:,2)*yaco(:); call joinBlock; yet=lvarr
     varr=zem(:,2)*yaco(:); call joinBlock; yze=lvarr
     varr=xim(:,3)*yaco(:); call joinBlock; zxi=lvarr
     varr=etm(:,3)*yaco(:); call joinBlock; zet=lvarr
     varr=zem(:,3)*yaco(:); call joinBlock; zze=lvarr

 do nn=1,3; do ip=0,1; i=ip*ijk(1,nn)
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
 select case(nn)
 case(1); rv(:)=yaco(l)*xim(l,:); fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3)); cm1(jk,:,ip)=fctr*rv(:)
 case(2); rv(:)=yaco(l)*etm(l,:); fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3)); cm2(jk,:,ip)=fctr*rv(:)
 case(3); rv(:)=yaco(l)*zem(l,:); fctr=one/sqrt(rv(1)*rv(1)+rv(2)*rv(2)+rv(3)*rv(3)); cm3(jk,:,ip)=fctr*rv(:)
 end select
 end do
 end do
 end do; end do

!===== EXTRA COEFFICIENTS FOR GCBC/GCIC

    cbca(:,:)=zero; cbca(1,1:2)=albed(1:2,0,0);
    cbca(2,1:3)=albed(0:2,1,0); cbca(3,1:3)=albed(-1:1,2,0)
 if(mbci>=4) then
    cbca(3,4)=albed(2,2,0)
    do i=4,mbci
       cbca(i,i-3:i)=(/beta,alpha,one,alpha/);
       if(i<mbci) then; cbca(i,i+1)=beta; end if
    end do
 end if
    rbci(:)=zero; rbci(1:3)=(/one,albed(-1,1,0),albed(-2,2,0)/)
    call mtrxi(cbca,cbcs,1,mbci); sbci(:)=-matmul(cbcs(:,:),rbci(:))
    ! rpt- New added
 !???????????????????
    fctr=pi/(mbci+1); res=zero
 do i=1,mbci; res=res+one
    sbci(i)=half*sbci(i)*(one+cos(res*fctr))
 end do
    lp=-1; ll=-1; rr(:,1)=zero
 do nn=1,3; do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn); iq=1-2*ip
    if((np-10)*(np-20)*(np-25)*(np-30)==0) then
       do k=0,ijk(3,nn); do j=0,ijk(2,nn); l=indx3(i,j,k,nn)
       if((np-20)*(np-25)==0) then
          lp=lp+1; call extrabcc(de(lp,1))
       end if
          ll=ll+1; res=one/yaco(l); rr(l,1)=rr(l,1)+one; rr(ll,2)=res; rr(ll,3)=l+sml
       do ii=1,mbci; l=indx3(i+iq*ii,j,k,nn)
          ll=ll+1; rr(l,1)=rr(l,1)+one; rr(ll,2)=res*sbci(ii); rr(ll,3)=l+sml
       end do
       end do; end do
    end if
 end do; end do
    lq=ll; allocate(rpex(0:lp),sbcc(0:lq))
 do ll=0,lp
    rpex(ll)=de(ll,1)
 end do
 do ll=0,lq; l=rr(ll,3)
    sbcc(ll)=rr(ll,2)/rr(l,1)
 end do
 !???????????????????
 end subroutine getMetrics

!====================================================================================
!=====READ RESTART FILE
!====================================================================================
 subroutine readRestart

    open(9,file=crestart,access='stream'); lh=0
    read(9,pos=nr*lh+1) niter; lh=lh+1
    read(9,pos=nr*lh+1) ndt; lh=lh+1
    read(9,pos=nr*lh+1) dt; lh=lh+1
    read(9,pos=nr*lh+1) dts; lh=lh+1
    read(9,pos=nr*lh+1) dte; lh=lh+1
    read(9,pos=nr*lh+1) timo; lh=lh+1
    lp=lpos(myid)+lh
    if ((tsam-timo)/tsam<0.05e0) then
       tsam=timo
    end if
 do m=1,5; lq=(m-1)*ltomb
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    read(9,pos=nr*(lp+lq+lio(j,k))+1) qo(l:l+lxi,m)
 end do; end do
 end do
    close(9)
 end subroutine readRestart

!====================================================================================
!=====  READ RAW RESTART
!====================================================================================
  subroutine rdIRsta()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp,offset
     integer :: fh,amode,qarr,iolen
     integer, dimension (4) :: gsizes,lsizes,starts
     real(k8) :: rbuf
     integer(k4) :: ibuf

      if (myid==0) then
         write(*,"('Reading restart file..')") 
      end if

     wrlen=5*(lmx+1)
     amode=MPI_MODE_RDONLY
     CALL MPI_TYPE_EXTENT(MPI_REAL8,iolen,ierr)

     gsizes(:)=(/mbijkl(:),5/)
     lsizes(:)=(/mpijkl(:),5/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,qarr,ierr) 
     CALL MPI_TYPE_COMMIT(qarr,ierr)
     

     CALL MPI_FILE_OPEN(bcom,crestart,amode,info,fh,ierr)
     lh=0

         offset=lh*iolen ! Iteration Number
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1; niter=ibuf
         offset=lh*iolen ! 10*(n/10)+1
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1; ndt=ibuf
         offset=lh*iolen ! Timestep
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; dt=rbuf
         offset=lh*iolen ! ?
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; dts=rbuf
         offset=lh*iolen ! ?
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; dte=rbuf
         offset=lh*iolen ! time
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; timo=rbuf

     disp=lh*iolen
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL8,qarr,'native',info,ierr)
     CALL MPI_FILE_READ_ALL(fh,qo,wrlen,MPI_REAL8,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(qarr,ierr)


  end subroutine rdIRsta
!====================================================================================
!=====READ GRID
!====================================================================================
 subroutine readGrid

         open(9,file='data/grid'//cnzone,access='stream'); lh=0
         lp=lpos(myid)
      do m=1,3; lq=(m-1)*ltomb
      do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         read(9,pos=nr*(lp+lq+lio(j,k))+1) xyz2(l:l+lxi,m)
      end do; end do
      end do
         close(9)
 end subroutine readGrid
!====================================================================================
!=====WRITE GRID
!====================================================================================
 subroutine writeGrid

         open(9,file='data/grid'//cnzone,access='stream'); lh=0
         lp=lpos(myid)
      do m=1,3; lq=(m-1)*ltomb
      do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         write(9,pos=nr*(lp+lq+lio(j,k))+1) xyz2(l:l+lxi,m)
      end do; end do
      end do
         close(9)
 end subroutine writeGrid
!====================================================================================
!=====  READ RAW GRID
!====================================================================================
  subroutine rdIGrid()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp
     integer :: fh,amode,garr
     integer, dimension (4) :: gsizes,lsizes,starts
     character(16) :: cout

     cout='data/grid'//cnzone

     wrlen=3*(lmx+1)
     amode=MPI_MODE_RDONLY

     gsizes(:)=(/mbijkl(:),3/)
     lsizes(:)=(/mpijkl(:),3/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,garr,ierr) 
     CALL MPI_TYPE_COMMIT(garr,ierr)
     
     disp=0

     CALL MPI_FILE_OPEN(bcom,cout,amode,info,fh,ierr)
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL8,garr,'native',info,ierr)
     CALL MPI_FILE_READ_ALL(fh,xyz2,wrlen,MPI_REAL8,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(garr,ierr)

     if(myid==mo(mb)) CALL MPI_FILE_DELETE(cout,info,ierr)
  end subroutine rdIGrid

!====================================================================================
!=====  WRITE RAW GRID
!====================================================================================
  subroutine wrIGrid()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp
     integer :: fh,amode,garr
     integer, dimension (4) :: gsizes,lsizes,starts
     character(16) :: cout

     cout='data/grid'//cnzone

     wrlen=3*(lmx+1)
     amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

     gsizes(:)=(/mbijkl(:),3/)
     lsizes(:)=(/mpijkl(:),3/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,garr,ierr) 
     CALL MPI_TYPE_COMMIT(garr,ierr)
     
     disp=0

     CALL MPI_FILE_OPEN(bcom,cout,amode,info,fh,ierr)
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL8,garr,'native',info,ierr)
     CALL MPI_FILE_WRITE_ALL(fh,xyz2,wrlen,MPI_REAL8,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(garr,ierr)

  end subroutine wrIGrid
!====================================================================================
!=====WRITE RESTART FILE
!====================================================================================
 subroutine writeRestart

         open(9,file='i'//crestart,access='stream'); lh=0
      if(myid==mo(mb)) then
         write(9,pos=nr*lh+1) niter; lh=lh+1
         write(9,pos=nr*lh+1) ndt; lh=lh+1
         write(9,pos=nr*lh+1) dt; lh=lh+1
         write(9,pos=nr*lh+1) dts; lh=lh+1
         write(9,pos=nr*lh+1) dte; lh=lh+1
         write(9,pos=nr*lh+1) timo; lh=lh+1
      else
         lh=lh+6
      end if
         lp=lpos(myid)+lh
      do m=1,5; lq=(m-1)*ltomb
      do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         write(9,pos=nr*(lp+lq+lio(j,k))+1) qb(l:l+lxi,m)
      end do; end do
      end do
         close(9)
 end subroutine writeRestart
!====================================================================================
!=====COMPUTE VARIABLE DERIVATIVES
!====================================================================================
 subroutine getDeri(n)
 implicit none
 integer, intent(in) :: n

     if(.not.allocated(fetxi))allocate(f(0:ltomb-1))
     if(.not.allocated(fxi))allocate(fxi(0:ltomb-1),fet(0:ltomb-1),fze(0:ltomb-1))
     if(.not.allocated(fetxi))allocate(fetxi(0:ltomb-1))
     if(.not.allocated(fzexi))allocate(fzexi(0:ltomb-1),fzeet(0:ltomb-1))
     if(.not.allocated(fzeetxi))allocate(fzeetxi(0:ltomb-1))

     varr=qo(:,n); call joinBlock; f=lvarr
     rr(:,1)=qo(:,n)
     m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     varr=rr(:,1); call joinBlock; fxi=lvarr
     varr=rr(:,2); call joinBlock; fet=lvarr
     varr=rr(:,3); call joinBlock; fze=lvarr
 
     rr(:,1)=rr(:,2)
     m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(1,1,m)
     varr=rr(:,1); call joinBlock; fetxi=lvarr
 
     rr(:,1)=rr(:,3)
     m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(2,1,m); call deriv(1,1,m)
     varr=rr(:,1); call joinBlock; fzexi=lvarr
     varr=rr(:,2); call joinBlock; fzeet=lvarr

     rr(:,1)=rr(:,2)
     m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(1,1,m)
     varr=rr(:,1); call joinBlock; fzeetxi=lvarr

 end subroutine getDeri
    

!====================================================================================
!=====COMPUTE VARIABLE DERIVATIVES
!====================================================================================
 subroutine interpolate(n,tol)
 implicit none
 integer, intent(in) :: n
 real(nr), intent(in) :: tol
 integer :: i,j,k
 real(nr),dimension(12) :: r
 real(nr) :: res
 integer :: m
 integer :: xi0,xi1,et0,et1,ze0,ze1
 integer :: l000,l010,l100,l110,l001,l011,l101,l111
 real(nr) :: x000,x010,x100,x110,x001,x011,x101,x111
 real(nr) :: xn00,xn10,xn0,xn01,xn11,xn1
 real(nr),dimension(3) :: xs,start,xin,hxi,hxn
 real(nr),dimension(3,3) :: jaco
 real(nr) :: err1,thisxi,thiset,thisze,xn,yn,zn
 real(nr) :: xxin,xetn,xzen
 real(nr) :: yxin,yetn,yzen
 real(nr) :: zxin,zetn,zzen
 if ((myid==0).and.(n==1)) then
    write(*,"('Interpolation procedure in progress')")
 end if

do k = 0, lzei
   write(*,*) myid,k
   do j = 0, leti
      do i = 0, lxii;l2=indx4(i,j,k,1)
         xs(:)=(/xyz2(l2,1),xyz2(l2,2),xyz2(l2,3)/)
         if (xs(1)<bounds(0,1)) then
             qb(l2,n)=outside(n)
         elseif (xs(1)>bounds(1,1)) then
             qb(l2,n)=outside(n)
         elseif (xs(2)<bounds(0,2)) then
             qb(l2,n)=outside(n)
         elseif (xs(2)>bounds(1,2)) then
             qb(l2,n)=outside(n)
         else
           if (n==1) then
                 if (i==0) then
                    start(1)=0
                 else
                    l=indx4(i-1,j,k,1)
                    start(1)=(ixis(l,1))
                 end if
                 start(1)=max(start(1),0.0_nr);
                 start(1)=min(start(1),real(lxio,nr))
                 if (j==0) then
                    start(2)=(nint(real(leto*j/letio,nr)))
                 else
                    l=indx4(i,j-1,k,1)
                    start(2)=(ixis(l,2))
                 end if
                 start(2)=max(start(2),0.0_nr);
                 start(2)=min(start(2),real(leto,nr))
                 if (k==0) then
                    start(3)=(nint(real(lzeo*k/lzeio,nr)))
                 else
                    l=indx4(i,j,k-1,1)
                    start(3)=(ixis(l,3))
                 end if
                 start(3)=max(start(3),0.0_nr);
                 start(3)=min(start(3),real(lzeo,nr))
              xin(:)=(/start(1),start(2),start(3)/)
              hxi(:)=(/0,0,0/)
              err1=1
              do while(err1>tol)
              xin=xin+hxi
              xin(:)=max(xin(:),(/0.0_nr,0.0_nr,0.0_nr/));
              xin(:)=min(xin(:),(/real(lxio,nr),real(leto,nr),real(lzeo,nr)/))
              thisxi=xin(1);thiset=xin(2);thisze=xin(3)
              if (mod(thisxi,1.0_nr)==0) then
                 if (thisxi==lxio) then
                    xi0=thisxi-1;xi1=thisxi;
                    else
                    xi0=thisxi;xi1=thisxi+1
                 end if
                 else
                 xi0=floor(thisxi);xi1=ceiling(thisxi)
              end if
              if (mod(thiset,1.0_nr)==0) then
                 if (thiset==leto) then
                    et0=thiset-1;et1=thiset;
                    else
                    et0=thiset;et1=thiset+1
                 end if
                 else
                 et0=floor(thiset);et1=ceiling(thiset)
              end if
              if (mod(thisze,1.0_nr)==0) then
                 if (thisze==lzeo) then
                    ze0=thisze-1;ze1=thisze;
                    else
                    ze0=thisze;ze1=thisze+1
                 end if
                 else
                 ze0=floor(thisze);ze1=ceiling(thisze)
              end if
              l000=indx5(xi0,et0,ze0,1);
              l010=indx5(xi0,et1,ze0,1);
              l100=indx5(xi1,et0,ze0,1);
              l110=indx5(xi1,et1,ze0,1);
              l001=indx5(xi0,et0,ze1,1);
              l101=indx5(xi1,et0,ze1,1);
              l011=indx5(xi0,et1,ze1,1);
              l111=indx5(xi1,et1,ze1,1);
              do m = 1, 12
                 selectcase(m)
                 case(1,2,3)
                 x000=xyz(l000,m);x010=xyz(l010,m);x100=xyz(l100,m);x110=xyz(l110,m);
                 x001=xyz(l001,m);x101=xyz(l101,m);x011=xyz(l011,m);x111=xyz(l111,m)
                 case(4)
                 x000=xxi(l000);x010=xxi(l010);x100=xxi(l100);x110=xxi(l110);
                 x001=xxi(l001);x101=xxi(l101);x011=xxi(l011);x111=xxi(l111)
                 case(5)
                 x000=xet(l000);x010=xet(l010);x100=xet(l100);x110=xet(l110);
                 x001=xet(l001);x101=xet(l101);x011=xet(l011);x111=xet(l111)
                 case(6)
                 x000=xze(l000);x010=xze(l010);x100=xze(l100);x110=xze(l110);
                 x001=xze(l001);x101=xze(l101);x011=xze(l011);x111=xze(l111)
                 case(7)
                 x000=yxi(l000);x010=yxi(l010);x100=yxi(l100);x110=yxi(l110);
                 x001=yxi(l001);x101=yxi(l101);x011=yxi(l011);x111=yxi(l111)
                 case(8)
                 x000=yet(l000);x010=yet(l010);x100=yet(l100);x110=yet(l110);
                 x001=yet(l001);x101=yet(l101);x011=yet(l011);x111=yet(l111)
                 case(9)
                 x000=yze(l000);x010=yze(l010);x100=yze(l100);x110=yze(l110);
                 x001=yze(l001);x101=yze(l101);x011=yze(l011);x111=yze(l111)
                 case(10)
                 x000=zxi(l000);x010=zxi(l010);x100=zxi(l100);x110=zxi(l110);
                 x001=zxi(l001);x101=zxi(l101);x011=zxi(l011);x111=zxi(l111)
                 case(11)
                 x000=zet(l000);x010=zet(l010);x100=zet(l100);x110=zet(l110);
                 x001=zet(l001);x101=zet(l101);x011=zet(l011);x111=zet(l111)
                 case(12)
                 x000=zze(l000);x010=zze(l010);x100=zze(l100);x110=zze(l110);
                 x001=zze(l001);x101=zze(l101);x011=zze(l011);x111=zze(l111)
                 end select
                 xn00=(xi1-thisxi)*x000+(thisxi-xi0)*x100;
                 xn10=(xi1-thisxi)*x010+(thisxi-xi0)*x110;
                 xn0=(et1-thiset)*xn00+(thiset-et0)*xn10;
                 xn01=(xi1-thisxi)*x001+(thisxi-xi0)*x101;
                 xn11=(xi1-thisxi)*x011+(thisxi-xi0)*x111;
                 xn1=(et1-thiset)*xn01+(thiset-et0)*xn11;
                 res=(ze1-thisze)*xn0+(thisze-ze0)*xn1;
                 select case (m)
                 case(1);xn=res;case(2);yn=res;case(3);zn=res;
                 case(4);xxin=res;case(5);xetn=res;case(6);xzen=res;
                 case(7);yxin=res;case(8);yetn=res;case(9);yzen=res;
                 case(10);zxin=res;case(11);zetn=res;case(12);zzen=res;
                 end select
              end do
              hxn=xs(:)-(/xn,yn,zn/)
              jaco(1,:)=(/xxin,yxin,zxin/)
              jaco(2,:)=(/xetn,yetn,zetn/)
              jaco(3,:)=(/xzen,yzen,zzen/)
              hxi(1)=sum(jaco(1,:)*hxn(:))
              hxi(2)=sum(jaco(2,:)*hxn(:))
              hxi(3)=sum(jaco(3,:)*hxn(:))
              err1=sqrt(hxn(1)**2+hxn(2)**2+hxn(3)**2)
              end do
           ixis(l2,:)=(/thisxi,thiset,thisze/)
           end if
         qb(l2,n)=htrilinr(ixis(l2,1),ixis(l2,2),ixis(l2,3))
         end if
      end do
   end do
end do
 if ((myid==mo(mb)).and.(n==5)) then
    write(*,"('Block ',i2,' done!')") mb
 end if
 end subroutine interpolate

!====================================================================================
!=====  WRITE RAW INTERPOLATED RESTART 
!====================================================================================
  subroutine wrIRsta()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp,offset
     integer :: amode,iolen
     integer, dimension (4) :: gsizes,lsizes,starts
     real(k8) :: rbuf
     integer(k4) :: ibuf

        wrlen=5*(lmx+1)
        if(.not.allocated(q8)) allocate(q8(0:lmx,5))
        q8(:,:)=qb(:,:)
        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
        CALL MPI_TYPE_EXTENT(MPI_REAL8,iolen,ierr)

        if (.not.qflag) then
           gsizes(:)=(/mbijkl(:),5/)
           lsizes(:)=(/mpijkl(:),5/)
           starts(:)=(/mpijks(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                            MPI_ORDER_FORTRAN,MPI_REAL8,qarr,ierr) 
           CALL MPI_TYPE_COMMIT(qarr,ierr)
           qflag=.true.
        end if
        

        if(myid==mo(mb)) CALL MPI_FILE_DELETE(icrestart,info,ierr)
        CALL MPI_FILE_OPEN(bcom,icrestart,amode,info,qfh,ierr)
        lh=0
        if (myid==mo(mb)) then
            ibuf=niter; offset=lh*iolen ! Iteration Number
            CALL MPI_FILE_WRITE_AT(qfh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ndt; offset=lh*iolen ! ?
            CALL MPI_FILE_WRITE_AT(qfh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            rbuf=dt; offset=lh*iolen ! Timestep
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
            rbuf=dts; offset=lh*iolen ! ?
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
            rbuf=dte; offset=lh*iolen ! ?
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
            rbuf=timo; offset=lh*iolen ! time
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
        else
            lh=lh+6
        end if
        disp=lh*iolen
        CALL MPI_FILE_SET_VIEW(qfh,disp,MPI_REAL8,qarr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL(qfh,q8,wrlen,MPI_REAL8,ista,ierr)
        CALL MPI_FILE_CLOSE(qfh,ierr)
        CALL MPI_TYPE_FREE(qarr,ierr)
        if (myid==0) then
           write(*,"('Restart file written!')") 
        end if

  end subroutine wrIRsta
!====================================================================================
! ====JOIN BLOCK DATA
!====================================================================================
 subroutine joinBlock
 implicit none
 integer :: i,j,k,l,lp

 if(.not.allocated(lvarr)) allocate(lvarr(0:ltomb-1),lvarr2(0:ltomb-1))
 lvarr=0_nr;lvarr2=0_nr

 do k = 0, lze
    do j = 0, let
       do i = 0, lxi; l=indx5(i,j,k,1) !RPT CHECK THIS!!!!!
          lp=i+lio(j,k)+lpos(myid)
          lvarr2(lp)=varr(l)
       end do
    end do
 end do

 CALL MPI_ALLREDUCE(lvarr2,lvarr,ltomb,MPI_REAL8,MPI_SUM,ncom,ierr)
    
 end subroutine joinBlock

!====================================================================================
!=====COPY OVER SPAN
!====================================================================================
 subroutine spanCopy()
 integer :: i,j,k,kk,l,l2

 do k = 0, lzei
    kk=k-(k/lze)*lze
    if(myid==0)write(*,*) k/lze,kk
    do j = 0, leti
       do i = 0, lxii; l=indx4(i,j,k,1); l2=indx3(i,j,kk,1)
          qb(l,:)=qo(l2,:)
       end do
    end do
 end do
    
 end subroutine spanCopy
!====================================================================================
!=====COPY OVER SPAN ITS SPAN-AVERAGE
!====================================================================================
 subroutine spanAvg()
 integer :: i,j,k,kk,l,l2,plmx

 plmx=nbsize(3)-1
 qa(0:plmx,:)=zero
 do k = 0, lze; l=k*(plmx+1)
    qa(0:plmx,:)=qa(0:plmx,:)+qo(l:l+plmx,:)
 end do
 ra1=1.0_k8/real(lze+1,k8)
 qa(0:plmx,:)=qa(0:plmx,:)*ra1
 do k = 0, lzei
    if(myid==0)write(*,*) k
    do j = 0, leti
       do i = 0, lxii; l=indx4(i,j,k,1); l2=indx3(i,j,0,1)
          qb(l,:)=qa(l2,:)
       end do
    end do
 end do
    
 end subroutine spanAvg

!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 3D

 function indx5(i,j,k,nn) result(lm)

 integer,intent(in) :: i,j,k,nn
 integer :: lm

 lm=i+lio(j,k)+lpos(myid)


 end function indx5
!===== FUNCTION FOR MAIN INDEX TRANSFORMATION IN 3D

 function indx4(i,j,k,nn) result(lm)

 integer,intent(in) :: i,j,k,nn
 integer :: lm

 select case(nn)
 case(1); lm=(k*(leti+1)+j)*(lxii+1)+i
 case(2); lm=(j*(leti+1)+i)*(lxii+1)+k
 case(3); lm=(i*(leti+1)+k)*(lxii+1)+j
 end select

 end function indx4

 
!===== FUNCTION FOR 3D HERMITEAN INTERPOLATION
 function htrilinr(xi,et,ze) result(fs)
 implicit none
 real(nr) :: fs
 real(nr), intent(in) :: xi,et,ze
 integer :: xi0,xi1,et0,et1,ze0,ze1
 integer ::  l000,l010,l100,l110,&
             l001,l011,l101,l111
 real(nr) :: f000,f010,f100,f110,&
             f001,f011,f101,f111
 real(nr) :: fxi000,fxi010,fxi100,fxi110,&
             fxi001,fxi011,fxi101,fxi111
 real(nr) :: fet000,fet010,fet100,fet110,&
             fet001,fet011,fet101,fet111
 real(nr) :: fze000,fze010,fze100,fze110,&
             fze001,fze011,fze101,fze111
 real(nr) :: fetxi000,fetxi010,fetxi100,fetxi110,&
             fetxi001,fetxi011,fetxi101,fetxi111
 real(nr) :: fzexi000,fzexi010,fzexi100,fzexi110,&
             fzexi001,fzexi011,fzexi101,fzexi111
 real(nr) :: fzeet000,fzeet010,fzeet100,fzeet110,&
             fzeet001,fzeet011,fzeet101,fzeet111
 real(nr) :: fzeetxi000,fzeetxi010,fzeetxi100,fzeetxi110,&
             fzeetxi001,fzeetxi011,fzeetxi101,fzeetxi111
 real(nr) :: f00,f10,fet00,fet10,f0,fzeet00,fzeet10
 real(nr) :: f01,f11,fet01,fet11,f1,fzeet01,fzeet11
 real(nr) :: fze00,fze10,fze01,fze11,fze0,fze1

 if (mod(xi,1.0_nr)==0) then
    if (xi==lxio) then
       xi0=xi-1;xi1=xi;
       else
       xi0=xi;xi1=xi+1
    end if
    else
    xi0=floor(xi);xi1=ceiling(xi)
 end if
 if (mod(et,1.0_nr)==0) then
    if (et==leto) then
       et0=et-1;et1=et;
       else
       et0=et;et1=et+1
    end if
    else
    et0=floor(et);et1=ceiling(et)
 end if
 if (mod(ze,1.0_nr)==0) then
    if (ze==lzeo) then
       ze0=ze-1;ze1=ze;
       else
       ze0=ze;ze1=ze+1
    end if
    else
    ze0=floor(ze);ze1=ceiling(ze)
 end if
 l000=indx5(xi0,et0,ze0,1);f000=f(l000)
 l010=indx5(xi0,et1,ze0,1);f010=f(l010)
 l100=indx5(xi1,et0,ze0,1);f100=f(l100)
 l110=indx5(xi1,et1,ze0,1);f110=f(l110)
 l001=indx5(xi0,et0,ze1,1);f001=f(l001)
 l101=indx5(xi1,et0,ze1,1);f101=f(l101)
 l011=indx5(xi0,et1,ze1,1);f011=f(l011)
 l111=indx5(xi1,et1,ze1,1);f111=f(l111)
 fxi000=fxi(l000);fet000=fet(l000);fze000=fze(l000);
 fxi010=fxi(l010);fet010=fet(l010);fze010=fze(l010);
 fxi100=fxi(l100);fet100=fet(l100);fze100=fze(l100);
 fxi110=fxi(l110);fet110=fet(l110);fze110=fze(l110);
 fxi001=fxi(l001);fet001=fet(l001);fze001=fze(l001);
 fxi011=fxi(l011);fet011=fet(l011);fze011=fze(l011);
 fxi101=fxi(l101);fet101=fet(l101);fze101=fze(l101);
 fxi111=fxi(l111);fet111=fet(l111);fze111=fze(l111);
 fetxi000=fetxi(l000);fzexi000=fzexi(l000);fzeet000=fzeet(l000);fzeetxi000=fzeetxi(l000);
 fetxi010=fetxi(l010);fzexi010=fzexi(l010);fzeet010=fzeet(l010);fzeetxi010=fzeetxi(l010);
 fetxi100=fetxi(l100);fzexi100=fzexi(l100);fzeet100=fzeet(l100);fzeetxi100=fzeetxi(l100);
 fetxi110=fetxi(l110);fzexi110=fzexi(l110);fzeet110=fzeet(l110);fzeetxi110=fzeetxi(l110);
 fetxi001=fetxi(l001);fzexi001=fzexi(l001);fzeet001=fzeet(l001);fzeetxi001=fzeetxi(l001);
 fetxi011=fetxi(l011);fzexi011=fzexi(l011);fzeet011=fzeet(l011);fzeetxi011=fzeetxi(l011);
 fetxi101=fetxi(l101);fzexi101=fzexi(l101);fzeet101=fzeet(l101);fzeetxi101=fzeetxi(l101);
 fetxi111=fetxi(l111);fzexi111=fzexi(l111);fzeet111=fzeet(l111);fzeetxi111=fzeetxi(l111);

 f00=fhermite(f000,fxi000,f100,fxi100,xi0,xi1,xi);
 f10=fhermite(f010,fxi010,f110,fxi110,xi0,xi1,xi);
 fet00=fhermite(fet000,fetxi000,fet100,fetxi100,xi0,xi1,xi);
 fet10=fhermite(fet010,fetxi010,fet110,fetxi110,xi0,xi1,xi);
 f0=fhermite(f00,fet00,f10,fet10,et0,et1,et);
 fzeet00=fhermite(fzeet000,fzeetxi000,fzeet100,fzeetxi100,xi0,xi1,xi);
 fzeet10=fhermite(fzeet010,fzeetxi010,fzeet110,fzeetxi110,xi0,xi1,xi);

 f01=fhermite(f001,fxi001,f101,fxi101,xi0,xi1,xi);
 f11=fhermite(f011,fxi011,f111,fxi111,xi0,xi1,xi);
 fet01=fhermite(fet001,fetxi001,fet101,fetxi101,xi0,xi1,xi);
 fet11=fhermite(fet011,fetxi011,fet111,fetxi111,xi0,xi1,xi);
 f1=fhermite(f01,fet01,f11,fet11,et0,et1,et);
 fzeet01=fhermite(fzeet001,fzeetxi001,fzeet101,fzeetxi101,xi0,xi1,xi);
 fzeet11=fhermite(fzeet011,fzeetxi011,fzeet111,fzeetxi111,xi0,xi1,xi);

 fze00=fhermite(fze000,fzexi000,fze100,fzexi100,xi0,xi1,xi);
 fze10=fhermite(fze010,fzexi010,fze110,fzexi110,xi0,xi1,xi);
 fze01=fhermite(fze001,fzexi001,fze101,fzexi101,xi0,xi1,xi);
 fze11=fhermite(fze011,fzexi011,fze111,fzexi111,xi0,xi1,xi);
 fze0=fhermite(fze00,fzeet00,fze10,fzeet10,et0,et1,et);
 fze1=fhermite(fze01,fzeet01,fze11,fzeet11,et0,et1,et);
 fs=fhermite(f0,fze0,f1,fze1,ze0,ze1,ze);

 end function htrilinr

!===== FUNCTION FOR INTERFACE POINTS
 function fhermite(k1,k2,k3,k4,x0,x1,x) result(y)

 real(nr) :: y
 integer, intent(in) :: x0,x1
 real(nr), intent(in) :: x,k1,k2,k3,k4
 real(nr) :: a,b,c,d,l,invl,fx

 l=x1-x0
 invl=1.0_nr/l
 fx=(x-x0)*invl

    a=(2*fx**3-3*fx**2+1);
    b=(fx**3-2*fx**2+fx);
    c=(-2*fx**3+3*fx**2);
    d=(fx**3-fx**2);
    y=a*k1+b*l*k2+c*k3+d*l*k4;
 
 end function fhermite

end module rptinter
