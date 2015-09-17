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
 integer :: l2,color,ncom

real(nr),dimension(:),allocatable :: xxi,xet,xze
real(nr),dimension(:),allocatable :: yxi,yet,yze
real(nr),dimension(:),allocatable :: zxi,zet,zze
real(nr),dimension(:),allocatable :: f
real(nr),dimension(:),allocatable :: fxi,fet,fze
real(nr),dimension(:),allocatable :: fetxi
real(nr),dimension(:),allocatable :: fzexi,fzeet,fzeetxi
real(nr),dimension(:),allocatable :: lvarr,lvarr2
real(nr) :: dti,dtsi,dtei,timoi,wtsi,wtei
integer :: ndti,ni

real(nr), dimension(0:1,1:3) :: bounds
real(nr), dimension(1:5) :: outside
contains

!====================================================================================
!=====PROBLEM SETUP
!====================================================================================
 subroutine setup(ilxi0,ilxi1,ilxi2,ilet0,ilet1,ilze0)
 integer, intent(in) :: ilxi0,ilxi1,ilxi2,ilet0,ilet1,ilze0 
!===== INPUT PARAMETERS

    open(9,file='inputo.dat',shared)
    read(9,*) cinput,mbk
    read(9,*) cinput,nts
    read(9,*) cinput,nscrn,nsgnl
    read(9,*) cinput,ndata
    read(9,*) cinput,nkrk
    read(9,*) cinput,nviscous
    read(9,*) cinput,nsmf
    read(9,*) cinput,nfskp
    read(9,*) cinput,nrestart
    read(9,*) cinput,reoo,tempoo
    read(9,*) cinput,amach1,amach2,amach3
    read(9,*) cinput,wtemp
    read(9,*) cinput,cfl
    read(9,*) cinput,tmax,timf,tsam
    read(9,*) cinput,fltk,fltkbc
    read(9,*) cinput,dto
    close(9)

    cinput=cinput; fltk=pi*fltk; fltkbc=pi*fltkbc
    rhooo=1; poo=1/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1**2+amach2**2+amach3**2)
    srefoo=111.0_nr/tempoo; srefp1dre=(srefoo+1)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=1/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo
    ! rpt-Initialising the record count 
    nwrec=0

    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),npc(0:mbk,3))

    call inputext

    lxi0=ilxi0
    lxi1=ilxi1
    lxi2=ilxi2
    let0=ilet0
    let1=ilet1
    lze0=ilze0

    lximb(0:5)=(/lxi0,lxi1,lxi2,lxi0,lxi1,lxi2/)
    lximb(6:11)=(/lxi0,lxi1,lxi2,lxi0,lxi1,lxi2/)
    letmb(0:5)=(/let0,let0,let0,let1,let1,let1/)
    letmb(6:11)=(/let1,let1,let1,let0,let0,let0/)
    lzemb(0:5)=(/lze0,lze0,lze0,lze0,lze0,lze0/)
    lzemb(6:11)=(/lze0,lze0,lze0,lze0,lze0,lze0/)
!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

    mo(0)=0
 do mm=1,mbk
    mo(mm)=mo(mm-1)+npc(mm-1,1)*npc(mm-1,2)*npc(mm-1,3)
 end do
 do mm=0,mbk
 if(myid>=mo(mm)) then; mb=mm; end if
 end do
    lxio=lximb(mb); leto=letmb(mb); lzeo=lzemb(mb)

    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10)
    cno=achar(no+48); cnzone=cno(2)//cno(1)//cno(0)
    czone='zone'//cnzone;
    coutput='out/output'//cnzone//'.plt'
    ctecout='out/tecout'//cnzone//'.plt'
    cgrid='misc/grid'//cnzone//'.dat';
    crestart='rsta/restart'//cnzone//'.dat'

    no(4)=myid/10000; no(3)=mod(myid,10000)/1000;
    no(2)=mod(myid,1000)/100; no(1)=mod(myid,100)/10; no(0)=mod(myid,10)
    cno=achar(no+48); cnnode=cno(4)//cno(3)//cno(2)//cno(1)//cno(0)
    cdata='misc/data'//cnnode//'.dat';
    cturb='misc/turb'//cnnode//'.dat'

    call domdcomp

    ip=mod(myid-mo(mb),npc(mb,1))
    jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

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
    l=ll; nbc(0,nn)=nbcs(nn); nbc(1,nn)=nbce(nn); ncd(0,nn)=ncds(nn); ncd(1,nn)=ncde(nn)
 end if
 if(ma>=2) then
 if(lp==0) then
    l=ll-((ll+1)/ma)*(ma-1); nbc(0,nn)=nbcs(nn); nbc(1,nn)=40; ncd(0,nn)=ncds(nn); ncd(1,nn)=myid+mp
 end if
 if(lp>0.and.lp<ma-1) then
    l=(ll+1)/ma-1; nbc(0,nn)=40; nbc(1,nn)=40; ncd(0,nn)=myid-mp; ncd(1,nn)=myid+mp
 end if
 if(lp==ma-1) then
    l=(ll+1)/ma-1; nbc(0,nn)=40; nbc(1,nn)=nbce(nn); ncd(0,nn)=myid-mp; ncd(1,nn)=ncde(nn)
 end if
 end if
 select case(nn); case (1); lxi=l; case (2); let=l; case (3); lze=l; end select
 end do

!===== SUBDOMAIN SIZES & WRITING START POSITIONS IN OUTPUT FILE

    lxim(myid)=lxi; letm(myid)=let; lzem(myid)=lze
 do mp=0,mpro
    call MPI_BCAST(lxim(mp),1,MPI_INTEGER,mp,icom,ierr)
    call MPI_BCAST(letm(mp),1,MPI_INTEGER,mp,icom,ierr)
    call MPI_BCAST(lzem(mp),1,MPI_INTEGER,mp,icom,ierr)
 end do

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
 do j=1,npc(mm,2)-1; do i=0,npc(mm,1)-1
    mp=mo(mm)+j*jp+i
    lpos(mp)=lpos(mp-jp)+(lximb(mm)+1)*(letm(mp-jp)+1)
 end do; end do
    kp=npc(mm,1)*npc(mm,2)
 do k=1,npc(mm,3)-1; do j=0,npc(mm,2)-1; do i=0,npc(mm,1)-1
    mp=mo(mm)+k*kp+j*jp+i
    lpos(mp)=lpos(mp-kp)+(lximb(mm)+1)*(letmb(mm)+1)*(lzem(mp-kp)+1)
 end do; end do; end do
 end do

    allocate(lio(0:let,0:lze))
 do k=0,lze; kp=k*(leto+1)*(lxio+1)
 do j=0,let; jp=j*(lxio+1)
    lio(j,k)=jp+kp
 end do
 end do

 end subroutine setup

!====================================================================================
!=====PREPARE ARRAYS
!====================================================================================
 subroutine deallocateArrays

 deallocate(iit,idsgnl,lsgnl)
 deallocate(qo,qa,de,xim,etm,zem,rr,ss,p,yaco,varr)
 deallocate(lximb,letmb,lzemb,lhmb,mo,npc)
 deallocate(mxc,ran,sit,ait,xit,yit,zit)
 deallocate(drva1,drva2,drva3)
 deallocate(drvb1,drvb2,drvb3)
 deallocate(send1,send2,send3)
 deallocate(recv1,recv2,recv3)
 deallocate(cm1,cm2,cm3)
 deallocate(xu,yu,xl,yl,li,sa,sb,lio)
    
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

!===== EXTRA COEFFICIENTS FOR GCBC/GCIC

    cbca(:,:)=0; cbca(1,1:2)=albed(1:2,0,0); cbca(2,1:3)=albed(0:2,1,0); cbca(3,1:3)=albed(-1:1,2,0)
 if(mbci>=4) then
    cbca(3,4)=albed(2,2,0)
 do i=4,mbci
    cbca(i,i-3:i)=(/beta,alpha,one,alpha/); if(i<mbci) then; cbca(i,i+1)=beta; end if
 end do
 end if
    rbci(:)=0; rbci(1:3)=(/one,albed(-1,1,0),albed(-2,2,0)/)
    call mtrxi(cbca,cbcs,1,mbci); sbci(:)=-matmul(cbcs(:,:),rbci(:))

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

 do nn=1,3
 select case(nn)
 case(1); is=0; ie=is+lxi; case(2); is=lxi+1; ie=is+let; case(3); is=lxi+let+2; ie=is+lze
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

  end subroutine prepareArrays

!====================================================================================
!=====CREATE GRID
!====================================================================================
  subroutine getGrid

    call makegrid
    call MPI_BARRIER(icom,ierr)

    open(9,file=cgrid,access='stream',shared)
    lp=lpos(myid)
 do nn=1,3; lq=(nn-1)*ltomb
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    read(9,pos=nr*(lp+lq+lio(j,k))+1) ss(l:l+lxi,nn)
 end do; end do
 end do
    close(9)
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
     m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)
 
     rr(:,1)=ss(:,2)
     m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)
 
     rr(:,1)=ss(:,3)
     m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     de(:,1)=rr(:,1); de(:,2)=rr(:,2); de(:,3)=rr(:,3)
 
     allocate(xxi(0:ltomb-1),xet(0:ltomb-1),xze(0:ltomb-1))
     allocate(yxi(0:ltomb-1),yet(0:ltomb-1),yze(0:ltomb-1))
     allocate(zxi(0:ltomb-1),zet(0:ltomb-1),zze(0:ltomb-1))

     !xxi=qo(:,1);xet=qo(:,2);xze=qo(:,3)
     !yxi=qa(:,1);yet=qa(:,2);yze=qa(:,3)
     !zxi=de(:,1);zet=de(:,2);zze=de(:,3)
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
    
 !===== COMPUTE JACOBIAN
     yaco(:)=3/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
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
 case(1); rv(:)=yaco(l)*xim(l,:); fctr=1/sqrt(rv(1)**2+rv(2)**2+rv(3)**2); cm1(jk,:,ip)=fctr*rv(:)
 case(2); rv(:)=yaco(l)*etm(l,:); fctr=1/sqrt(rv(1)**2+rv(2)**2+rv(3)**2); cm2(jk,:,ip)=fctr*rv(:)
 case(3); rv(:)=yaco(l)*zem(l,:); fctr=1/sqrt(rv(1)**2+rv(2)**2+rv(3)**2); cm3(jk,:,ip)=fctr*rv(:)
 end select
 end do
 end do
 end do; end do
 end subroutine getMetrics

!====================================================================================
!=====READ RESTART FILE
!====================================================================================
 subroutine readRestart

    open(9,file=crestart,access='stream',shared); lh=0
    read(9,pos=nr*lh+1) n; lh=lh+1
    read(9,pos=nr*lh+1) ndt; lh=lh+1
    read(9,pos=nr*lh+1) dt; lh=lh+1
    read(9,pos=nr*lh+1) dts; lh=lh+1
    read(9,pos=nr*lh+1) dte; lh=lh+1
    read(9,pos=nr*lh+1) timo; lh=lh+1
    ni=n;ndti=ndt;dti=dt;dtsi=dts;dtei=dte;timoi=timo;
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
!=====READ GRID
!====================================================================================
 subroutine readGrid

         open(9,file='data/grid'//cnzone,access='stream',shared); lh=0
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

         open(9,file='data/grid'//cnzone,access='stream',shared); lh=0
         lp=lpos(myid)
      do m=1,3; lq=(m-1)*ltomb
      do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         write(9,pos=nr*(lp+lq+lio(j,k))+1) xyz2(l:l+lxi,m)
      end do; end do
      end do
         close(9)
 end subroutine writeGrid
!====================================================================================
!=====WRITE RESTART FILE
!====================================================================================
 subroutine writeRestart

         open(9,file='i'//crestart,access='stream',shared); lh=0
      if(myid==mo(mb)) then
         write(9,pos=nr*lh+1) ni; lh=lh+1
         write(9,pos=nr*lh+1) ndti; lh=lh+1
         write(9,pos=nr*lh+1) dti; lh=lh+1
         write(9,pos=nr*lh+1) dtsi; lh=lh+1
         write(9,pos=nr*lh+1) dtei; lh=lh+1
         write(9,pos=nr*lh+1) timoi; lh=lh+1
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
     m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     varr=rr(:,1); call joinBlock; fxi=lvarr
     varr=rr(:,2); call joinBlock; fet=lvarr
     varr=rr(:,3); call joinBlock; fze=lvarr
 
     rr(:,1)=rr(:,2)
     m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(1,1)
     varr=rr(:,1); call joinBlock; fetxi=lvarr
 
     rr(:,1)=rr(:,3)
     m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(2,1); call deriv(1,1)
     varr=rr(:,1); call joinBlock; fzexi=lvarr
     varr=rr(:,2); call joinBlock; fzeet=lvarr

     rr(:,1)=rr(:,2)
     m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(1,1)
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
         !elseif (xs(3)<bounds(0,3)) then
         !    qb(l2,n)=outside(n)
         !elseif (xs(3)>bounds(1,3)) then
         !    qb(l2,n)=outside(n)
         else
           if (n==1) then
                 if (i==0) then
                    start(1)=0
                 else
                    l=indx4(i-1,j,k,1)
                    start(1)=(ixis(l,1))
                    if (myid==0) then
                       write(*,*) start(1) 
                    end if
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
! ====READ DATA FOR POST-PROCESSING
!====================================================================================
 subroutine interead(num)
 implicit none
 integer, intent (in) :: num
 integer :: lp,lq,l,k,j
  lp=0
     lq=(num-1)*ltomb
        read(8,pos=nr*(lp+lq)+1) varr(:)
 end subroutine interead


!====================================================================================
! ====READ DATA FOR POST-PROCESSING
!====================================================================================
 subroutine joinBlock
 implicit none
 integer :: i,j,k,l,lp

 if(.not.allocated(lvarr)) allocate(lvarr(0:ltomb-1),lvarr2(0:ltomb-1))
 lvarr=0_nr;lvarr2=0_nr

 do k = 0, lze
    do j = 0, let
       do i = 0, lxi; l=indx3(i,j,k,1)
          lp=i+lio(j,k)+lpos(myid)
          lvarr2(lp)=varr(l)
       end do
    end do
 end do

 CALL MPI_ALLREDUCE(lvarr2,lvarr,ltomb,MPI_REAL8,MPI_SUM,ncom,ierr)
    
 end subroutine joinBlock
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

 f00=hermite(f000,fxi000,f100,fxi100,xi0,xi1,xi);
 f10=hermite(f010,fxi010,f110,fxi110,xi0,xi1,xi);
 fet00=hermite(fet000,fetxi000,fet100,fetxi100,xi0,xi1,xi);
 fet10=hermite(fet010,fetxi010,fet110,fetxi110,xi0,xi1,xi);
 f0=hermite(f00,fet00,f10,fet10,et0,et1,et);
 fzeet00=hermite(fzeet000,fzeetxi000,fzeet100,fzeetxi100,xi0,xi1,xi);
 fzeet10=hermite(fzeet010,fzeetxi010,fzeet110,fzeetxi110,xi0,xi1,xi);

 f01=hermite(f001,fxi001,f101,fxi101,xi0,xi1,xi);
 f11=hermite(f011,fxi011,f111,fxi111,xi0,xi1,xi);
 fet01=hermite(fet001,fetxi001,fet101,fetxi101,xi0,xi1,xi);
 fet11=hermite(fet011,fetxi011,fet111,fetxi111,xi0,xi1,xi);
 f1=hermite(f01,fet01,f11,fet11,et0,et1,et);
 fzeet01=hermite(fzeet001,fzeetxi001,fzeet101,fzeetxi101,xi0,xi1,xi);
 fzeet11=hermite(fzeet011,fzeetxi011,fzeet111,fzeetxi111,xi0,xi1,xi);

 fze00=hermite(fze000,fzexi000,fze100,fzexi100,xi0,xi1,xi);
 fze10=hermite(fze010,fzexi010,fze110,fzexi110,xi0,xi1,xi);
 fze01=hermite(fze001,fzexi001,fze101,fzexi101,xi0,xi1,xi);
 fze11=hermite(fze011,fzexi011,fze111,fzexi111,xi0,xi1,xi);
 fze0=hermite(fze00,fzeet00,fze10,fzeet10,et0,et1,et);
 fze1=hermite(fze01,fzeet01,fze11,fzeet11,et0,et1,et);
 fs=hermite(f0,fze0,f1,fze1,ze0,ze1,ze);

 end function htrilinr

!===== FUNCTION FOR INTERFACE POINTS
 function hermite(k1,k2,k3,k4,x0,x1,x) result(y)

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
 
 end function hermite

end module rptinter
