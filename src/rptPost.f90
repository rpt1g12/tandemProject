!***************************
!***** RPT POST MODULE *****
!***************************

module  rptpost

use mainvar3d
use subroutineso
use subroutines3d
use problemcase
use mpi
use rpt

contains

!====================================================================================
!=====PROBLEM SETUP
!====================================================================================
 subroutine setup
!===== INPUT PARAMETERS

    open(9,file='inputo.dat',shared)
    read(9,*) cinput,mbk
    read(9,*) cinput,nts,nto,iwrec
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
    read(9,*) cinput,forcing,amfor
    read(9,*) cinput,aoa
    read(9,*) cinput,LES,smago1,smago2
    read(9,*) cinput,output,ogrid,osol,oblock
    close(9)

!===== INPUT PARAMETERS POSTPROCESS
    open(9,file='ipost.dat',shared)
    read(9,*) cinput,fparallel,fmblk
    read(9,*) cinput,favg,fwavg,favgu
    read(9,*) cinput,fcoef,fcf,fcp
    read(9,*) cinput,floc
    read(9,*) cinput,fwplus
    read(9,*) cinput,fqcrit,fwss
    read(9,*) cinput,fcurl
    read(9,*) cinput,frms,fwrms
    read(9,*) cinput,fstrip
    close(9)

    cinput=cinput; fltk=pi*fltk; fltkbc=pi*fltkbc
    rhooo=1; poo=1/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1**2+amach2**2+amach3**2)
    srefoo=111.0_k8/tempoo; srefp1dre=(srefoo+1)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=1/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo
    ! rpt-Initialising the record count 
    nwrec=0

    abc(:,0)=(/a01,a02,a03,a04,a05,a06/)
    abc(:,1)=(/a10,a12,a13,a14,a15,a16/)
    abc(:,2)=(/a20,a21,a23,a24,a25,a26/)

    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),npc(0:mbk,3))

    call inputext
    if(fparallel==0) npc(:,:)=1

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
    coutput='out/output'//cnzone//'.plt'
    cgrid='misc/grid'//cnzone//'.dat';
    crestart='rsta/restart'//cnzone//'.dat'

    no(4)=myid/10000; no(3)=mod(myid,10000)/1000;
    no(2)=mod(myid,1000)/100; no(1)=mod(myid,100)/10; no(0)=mod(myid,10)
    cno=achar(no+48); cnnode=cno(4)//cno(3)//cno(2)//cno(1)//cno(0)
    cdata='data/data'//cnnode//'.dat';
    cturb='misc/turb'//cnnode//'.dat'

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

!===== ALLOCATION OF MAIN ARRAYS

    allocate(qo(0:lmx,5),qa(0:lmx,5),qb(0:lmx,5),de(0:lmx,5))
    allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
    allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx))

 if(nviscous==1) then
    allocate(txx(0:lmx),tyy(0:lmx),tzz(0:lmx))
    allocate(txy(0:lmx),tyz(0:lmx),tzx(0:lmx))
    allocate(hxx(0:lmx),hyy(0:lmx),hzz(0:lmx))
 end if

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

    cbca(:,:)=0;                cbca(1,1:2)=albed(1:2,0,0);
    cbca(2,1:3)=albed(0:2,1,0); cbca(3,1:3)=albed(-1:1,2,0)
 if(mbci>=4) then
    cbca(3,4)=albed(2,2,0)
 do i=4,mbci
          cbca(i,i-3:i)=(/beta,alpha,one,alpha/);
          if(i<mbci) cbca(i,i+1)=beta
 end do
 end if
    rbci(:)=0; rbci(1:3)=(/one,albed(-1,1,0),albed(-2,2,0)/)
    call mtrxi(cbca,cbcs,1,mbci); sbci(:)=-matmul(cbcs(:,:),rbci(:))

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


 end subroutine setup

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
 
     allocate(xyz(0:lmx,3))
     do i = 1, 3
        xyz(:,i)=ss(:,i)
     end do
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

    call wallArea
    call walldir

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
!=====COMPUTE Y+,X+ AND Z+
!====================================================================================
 subroutine getwplus(nvar)
 implicit none
 integer, intent(in) :: nvar
 integer :: l,lp1,lm1,lw
 real(k8), dimension(0:lcwall) :: utau
 if (wflag) then
    if(.not.allocated(wplus)) allocate(wplus(0:lcwall,3))
    call gettw(nvar)
    ! COMPUTE Y+ AND UTAU
    lw=-1
    do k = 0, lze
       do i = 0, lxi;l=indx3(i,0,k,1);lp1=indx3(i,1,k,1);lw=lw+1
          utau(lw)=sqrt(tw(lw,1)**2+tw(lw,2)**2+tw(lw,3)**2)/qo(l,1)
          utau(lw)=sqrt(utau(lw))*reoo
          wplus(lw,2)=sqrt((xyz(lp1,1)-xyz(l,1))**2+(xyz(lp1,2)-xyz(l,2))**2)*utau(lw)
       end do
    end do
    ! COMPUTE X+
    lw=-1
    do k = 0, lze
    l=indx3(0,0,k,1);lp1=indx3(1,0,k,1);lw=lw+1
    wplus(lw,1)=sqrt((xyz(lp1,1)-xyz(l,1))**2+(xyz(lp1,2)-xyz(l,2))**2)*utau(lw)
       do i = 1, lxi-1;lp1=indx3(i+1,0,k,1);lm1=indx3(i-1,0,k,1);lw=lw+1
          wplus(lw,1)=sqrt((xyz(lp1,1)-xyz(lm1,1))**2+(xyz(lp1,2)-xyz(lm1,2))**2)*utau(lw)*half
       end do
    l=indx3(lxi,0,k,1);lm1=indx3(lxi-1,0,k,1);lw=lw+1
    wplus(lw,1)=sqrt((xyz(l,1)-xyz(lm1,1))**2+(xyz(l,2)-xyz(lm1,2))**2)*utau(lw)
    end do
    ! COMPUTE Z+
    lw=-1
    do i = 0, lxi;lp1=indx3(i,0,1,1);l=indx3(i,0,0,1);lw=lw+1
       wplus(lw,3)=abs((xyz(lp1,3)-xyz(l,3)))*utau(lw)
    end do
    do k = 1, lze-1
       do i = 0, lxi;lp1=indx3(i,0,k+1,1);lm1=indx3(i,0,k-1,1);lw=lw+1
          wplus(lw,3)=abs((xyz(lp1,3)-xyz(lm1,3)))*utau(lw)*half
       end do
    end do
    do i = 0, lxi;lm1=indx3(i,0,lze-1,1);l=indx3(i,0,lze,1);lw=lw+1
       wplus(lw,3)=abs((xyz(l,3)-xyz(lm1,3)))*utau(lw)
    end do
 end if
 end subroutine getwplus

!====================================================================================
!=====GET THE VALUES OF VARR AT WALL AND SAVE THEM IN WVARR
!====================================================================================
 subroutine getatWall
    implicit none
    integer :: n,l
    if (wflag) then
    if(.not.allocated(wvarr)) allocate(wvarr(0:lcwall))
    do n = 0, lcwall; l=lwall(n)
       wvarr(n)=varr(l)
    end do
    end if
 end subroutine getatWall

!====================================================================================
!=====GET THE VALUES OF WVARR FROM WALL AND SAVE THEM IN VARR
!====================================================================================
 subroutine getfromWall
    implicit none
    integer :: n,l

    varr(:)=0
    if (wflag) then
    do n = 0, lcwall; l=lwall(n)
       varr(l)=wvarr(n)
    end do
    end if
 end subroutine getfromWall

!====================================================================================
!=====AVERAGE IN SPACE OVER WALL SURFACE AND WRITE TO fname.dat
!====================================================================================
 subroutine wavg(fname,dir,wall)
    implicit none
    character(*), intent(in) :: fname
    integer, intent(in) :: dir
    logical, intent(in) :: wall
    integer :: idir,odir,ia,oa,l,ls,le,lp1,lm1,lw,lws,lwe,lfn
    real(k8), dimension(:), allocatable :: delt
    real(k8), dimension(:,:), allocatable :: avg
    character(1) :: str
    character(*), parameter :: cpath='out/',cext='.dat'
    character(len=:), allocatable :: lfname

    lfn=len(fname)
    l=len(cpath)+len(cext)+lfn
    allocate(character(len=l) :: lfname)
    lfname=cpath//fname//cext
    if (wflag) then
       if (.not.wall) call getatWall
       select case(dir)
       case(1); idir=lxi; odir=lze; ia=1; oa=3; str='z'
       case(2); idir=lze; odir=lxi; ia=3; oa=1; str='x'
       end select
       allocate(delt(0:idir))
       allocate(avg(0:odir,2))
       ! OPEN FILE TO OUTPUT RESUT
       open(7,file=lfname)
       write(7,*) str//' '//fname
       do j = 0, odir
          lws=indx2(0,j,dir);lwe=indx2(idir,j,dir)
          ls=lwall(lws);le=lwall(lwe)
          lp1=lwall(indx2(1,j,dir));lm1=lwall(indx2(idir-1,j,dir))
          fctr=half/abs(xyz(le,ia)-xyz(ls,ia))
          delt(0)=fctr*abs(xyz(lp1,ia)-xyz(ls,ia))
          delt(idir)=fctr*abs(xyz(le,ia)-xyz(lm1,ia))
          avg(j,2)=delt(0)*wvarr(lws)+delt(idir)*wvarr(lwe)
          avg(j,1)=xyz(ls,oa)
          do i = 1, idir-1;l=lwall(indx2(i,j,dir));lw=indx2(i,j,dir)
             lm1=lwall(indx2(i-1,j,dir));lp1=lwall(indx2(i+1,j,dir))
             delt(i)=fctr*abs(xyz(lm1,ia)-xyz(lp1,ia))
             avg(j,2)=avg(j,2)+delt(i)*wvarr(lw)
          end do
        write(7,"(f10.5,1x,f10.5)") avg(j,1),avg(j,2)
       end do
       close(7)
    end if

 end subroutine wavg
 
!====================================================================================
!===== SUBROUTINE FOR CALCULATING VORTICITY
!====================================================================================
 subroutine getCurl(nvar)

        implicit none
        integer, intent(in) :: nvar

        !call rdP3dS(nvar,fmblk)
        !p(:)=qo(:,5)
        de(:,1:3)=0

        rr(:,1)=qo(:,2)
        m=1; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1); call deriv(2,1); call deriv(1,1)
        de(:,2)=de(:,2)+rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)
        de(:,3)=de(:,3)-rr(:,1)*xim(:,2)-rr(:,2)*etm(:,2)-rr(:,3)*zem(:,2)

        rr(:,1)=qo(:,3)
        m=2; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1); call deriv(2,1); call deriv(1,1)
        de(:,3)=de(:,3)+rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
        de(:,1)=de(:,1)-rr(:,1)*xim(:,3)-rr(:,2)*etm(:,3)-rr(:,3)*zem(:,3)

        rr(:,1)=qo(:,4)
        m=3; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1); call deriv(2,1); call deriv(1,1)
        de(:,1)=de(:,1)+rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
        de(:,2)=de(:,2)-rr(:,1)*xim(:,1)-rr(:,2)*etm(:,1)-rr(:,3)*zem(:,1)

        qo(:,2)=de(:,1)*yaco(:);
        qo(:,3)=de(:,2)*yaco(:);
        qo(:,4)=de(:,3)*yaco(:)
        !ra0=aoa*pi/180;ra1=cos(ra0);ra2=sin(ra0)
        !ss(:,1)=de(:,1)!*ra1+de(:,2)*ra2
        !ss(:,2)=de(:,2)!*ra1-de(:,1)*ra2
        !ss(:,3)=de(:,3)

 end subroutine getCurl
!====================================================================================
!=====COMPUTE Q-CRITERION
!====================================================================================
 subroutine qcriterion(nvar)
    implicit none
    integer, intent(in) :: nvar

    call rdP3dS(nvar,fmblk)
    p(:)=qo(:,5)
    if (nviscous==1) then
     de(:,2)=qo(:,2)
     de(:,3)=qo(:,3)
     de(:,4)=qo(:,4)
 
     rr(:,1)=de(:,2)
     m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     txx(:)=yaco(:)*(xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3))
     hzz(:)=yaco(:)*(xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3))
     tzx(:)=yaco(:)*(xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3))
 
     rr(:,1)=de(:,3)
     m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     txy(:)=yaco(:)*(xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3))
     tyy(:)=yaco(:)*(xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3))
     hxx(:)=yaco(:)*(xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3))
 
     rr(:,1)=de(:,4)
     m=4; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     hyy(:)=yaco(:)*(xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3))
     tyz(:)=yaco(:)*(xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3))
     tzz(:)=yaco(:)*(xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3))

     qo(:,1)=2*(half*(hzz(:)-txy(:)))**2 + 2*(half*(tzx(:)-tzz(:)))**2 + &
     2*(half*(hxx(:)-tyz(:)))**2 - txx(:)**2 - tyy(:)**2 - tzz(:)**2 -   &
     2*(half*(hzz(:)+txy(:)))**2 - 2*(half*(tzx(:)+hyy(:)))**2 -         &
     2*(half*(hxx(:)+tyz(:)))**2
    end if

 end subroutine qcriterion

!====================================================================================
!=====FIND INDEX FROM X,Y, AND Z COORDINATES
!====================================================================================
 subroutine findll(xpos,ypos,zpos,l,id) 
    implicit none
    real(k8), intent(in) :: xpos,ypos,zpos
    integer(k4), intent (out) :: l,id
    real(k8) :: tmp,tmpall

    id=-1
    varr(:)=sqrt((xyz(:,1)-xpos)**2+(xyz(:,2)-ypos)**2+(xyz(:,3)-zpos)**2)
    l=minloc(varr(:),1)-1
    tmp=varr(l)
    CALL MPI_ALLREDUCE(tmp,tmpall,1,MPI_REAL8,MPI_MIN,icom,ierr)
    if (abs(tmp-tmpall)/tmp<sml) then
       id=myid
    end if

 end subroutine findll
 
!====================================================================================
!=====FILL ARRAY QO(:,:) WITH VALUES FROM OUTPUT = NVAR
!====================================================================================
 subroutine fillqo(nvar)
    implicit none
    integer, intent(in) :: nvar

    call rdP3dS(nvar,fmblk)
    !call p3dread(0,nvar)
    p(:)=qo(:,5)
    
 end subroutine fillqo

!====================================================================================
!=====WRITE LIST OF SOLUTION FILES
!====================================================================================
subroutine flst(mblkin)
   implicit none
   integer, intent(in),optional :: mblkin
   character(10)  :: ctime
   character(3) :: cout
   character(*), parameter :: cs1='ls out/*'
   character(:), allocatable :: cflst,rcout,ofile
   integer :: lmt,stat,i,mblk,foper,wrcom,err
   real(4) :: res

   if(present(mblkin)) then
      mblk=mblkin
   else
      mblk=1
   end if

   selectcase(mblk);
   case(1)
      cout=''
      foper=0
      wrcom=icom
   case(0)
      write(cout,"(a,i2)") 'b',mb
      do i = 0, 1
         l=scan(cout,' ')
         if (l==0) exit
         cout(l:l)='0'
      end do
      foper=mo(mb)
      wrcom=bcom
   case default
      if(myid==0) write(*,*) 'Wrong multiblock option! Aborting...'
      CALL MPI_ABORT(icom,err,ierr)
   end select
   l=len(trim(cout))
   allocate(character(l) :: rcout)
   rcout=trim(cout)

   if (myid==foper) then
      write(*,*) 'Delete previous filelist '//rcout
      open(90,file='out/filelist'//rcout); close(90,status='delete')
      call system('ls out/*'//rcout//'.q -1 > out/filelist'//rcout)
   end if
   CALL MPI_BARRIER(icom,ierr)
   if (myid==foper) then
      open(90,file='out/filelist'//rcout)
      inquire(unit=90,size=l)
      if (l.le.1) then
         close(90,status='delete')
      else
         close(90)
      end if   
   end if
   CALL MPI_BARRIER(icom,ierr)
   
   inquire(file='out/filelist'//rcout,exist=fflag)
   if (fflag) then
      open(90,file='out/filelist'//rcout,shared)
      lmt=-2
      do while (stat.ge.0)
         read(90,*,IOSTAT=stat) 
         lmt=lmt+1
      end do
      close(90)
      ndata=lmt

      l=len('out/sotT')+8+len(trim(rcout))+len('.q')
      allocate(character(l) :: ofiles(0:lmt))
      allocate(character(l) :: ofile)
      allocate(times(0:ndata))

      open(90,file='out/filelist'//rcout,shared)
      do i = 0, lmt
         read(90,"(a)") ofile
         ofiles(i)=ofile
         ctime=(ofile(9:16)//'e0')
         read(ctime,*) times(i)
      end do
      close(90)
      
      !do i = 0, lmt
      !   open(91,file=ofiles(i),access='stream',shared)
      !       lh=1+(mbk+1)*3+3
      !       read(91,pos=4*lh+1) res; times(i)=res
      !   close(91)
      !end do
   else
      l=1
      ndata=0
      allocate(character(l) :: ofiles(0:ndata))
      allocate(character(l) :: ofile)
      allocate(times(0:ndata))
      times=0
      ndata=-1
   end if

end subroutine flst

!====================================================================================
!=====AVERAGE RESULTS IN TIME PLOT3D
!====================================================================================
 subroutine p3daverage
    implicit none
    integer :: foper
    real(k8),dimension(:),allocatable :: delt
    character(3) :: cout

    selectcase(fmblk)
    case(0)
       write(cout,"(a,i2)") 'b',mb
       do i = 0, 1
          l=scan(cout,' ')
          if (l==0) exit
          cout(l:l)='0'
       end do
       foper=mo(mb)
    case(1)
       cout=''
       foper=0
    end select

    if (fflag) then
       if (myid==foper) then
          write(*,"('Total amout of data: ',i3,a)") ndata,cout
       end if
       ! CONSTRUCT THE COEFFICIENTS ARRAY
          ns=0; ne=ndata; allocate(delt(ns:ne))
          fctr=half/(times(ne)-times(ns))
          delt(ns)=fctr*(times(ns+1)-times(ns)); 
          delt(ne)=fctr*(times(ne)-times(ne-1))
       do n=ns+1,ne-1
          delt(n)=fctr*(times(n+1)-times(n-1))
       end do
          qa(:,:)=0
       do n=0,ndata
       if (myid==foper) then
          write(*,"(f5.1,'% Averaged',a)") real(n)*100.0e0/real(ndata),cout
       end if
          call rdP3dS(n,fmblk)
          !call p3dread(gsflag=0,nout=n)
          qa(:,:)=qa(:,:)+delt(n)*qo(:,:)
       end do
    end if
 end subroutine p3daverage

!====================================================================================
!=====RMS OF RESULTS IN TIME PLOT3D
!====================================================================================
 subroutine p3drms
    implicit none
    integer :: foper
    real(k8),dimension(:),allocatable :: delt
    character(3) :: cout

    selectcase(fmblk)
    case(0)
       write(cout,"(a,i2)") 'b',mb
       do i = 0, 1
          l=scan(cout,' ')
          if (l==0) exit
          cout(l:l)='0'
       end do
       foper=mo(mb)
    case(1)
       cout=''
       foper=0
    end select

    if (fflag) then
       if (myid==foper) then
          write(*,"('Total amout of data: ',i3,a)") ndata,cout
       end if
       ! CONSTRUCT THE COEFFICIENTS ARRAY
          ns=0; ne=ndata; allocate(delt(ns:ne))
          fctr=half/(times(ne)-times(ns))
          delt(ns)=fctr*(times(ns+1)-times(ns)); 
          delt(ne)=fctr*(times(ne)-times(ne-1))
       do n=ns+1,ne-1
          delt(n)=fctr*(times(n+1)-times(n-1))
       end do
          call rdP3dS(ndata+1,fmblk)
          !call p3dread(gsflag=0,nout=ndata+1)
          qa(:,:)=qo(:,:)
          qb(:,:)=0
       do n=0,ndata
          if (myid==foper) then
             write(*,"(f5.1,'% Done',a)") real(n)*100.0e0/real(ndata),cout
          end if
          !call p3dread(gsflag=0,nout=n)
          call rdP3dS(n,fmblk)
          de(:,:)=(qo(:,:)-qa(:,:))
          qb(:,:)=qb(:,:)+delt(n)*de(:,:)*de(:,:)
       end do
       qb(:,:)=sqrt(qb(:,:))
    end if
 end subroutine p3drms
!====================================================================================
!=====WRITE AVERAGE RESULTS IN TIME PLOT3D
!====================================================================================
 subroutine p3dwaverage()
 integer :: n

       if (myid==0) then
         open(9,file='out/solA.qa'); close(9,status='delete')
       end if
       CALL MPI_BARRIER(icom,ierr)
       open (unit=9, file='out/solA.qa', access='stream',shared)
       lh=0
       if (myid==0) then
        write(9,pos=4*lh+1) mbk+1; lh=lh+1 ! Number of zones
        do mm = 0, mbk
           write(9,pos=4*lh+1) int4(lximb(mm)+1); lh=lh+1 ! IMax
           write(9,pos=4*lh+1) int4(letmb(mm)+1); lh=lh+1 ! JMax
           write(9,pos=4*lh+1) int4(lzemb(mm)+1); lh=lh+1 ! KMax
        end do
        lhmb(mb)=lh
        do mm = 0, mbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do
       end if
        call MPI_BCAST(lhmb,mbk+1,MPI_INTEGER,0,icom,ierr)
        if (myid==mo(mb)) then
           lh=lhmb(mb)
           write(9,pos=4*lh+1) real(amachoo,kind=4); lh=lh+1 ! Mach Number
           write(9,pos=4*lh+1) real(aoa,kind=4); lh=lh+1  
           write(9,pos=4*lh+1) real(reoo,kind=4); lh=lh+1 ! Reynolds Number
           write(9,pos=4*lh+1) real(times(ndata),kind=4); lh=lh+1 ! Time
        end if
        lp=lpos(myid)+lhmb(mb)+4
        ns=1; ne=5
        do n=ns,ne; lq=(n-ns)*ltomb
           selectcase(n)
           case(1,5); varr(:)=qa(:,n)
           case(2,3,4); varr(:)=qa(:,n)!*qa(:,1)
           end select
        do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
          write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
        end do; end do
        end do
        close(9)
        CALL MPI_BARRIER(icom,ierr)
        if (myid==0) then
           write(*,"('Averaged Solution written!')") 
        end if
 end subroutine p3dwaverage


!====================================================================================
!=====WRITE RMS RESULTS IN TIME PLOT3D
!====================================================================================
 subroutine p3dwrms()
 integer :: n

       if (myid==0) then
         open(9,file='out/solRMS.qa'); close(9,status='delete')
       end if
       CALL MPI_BARRIER(icom,ierr)
       open (unit=9, file='out/solRMS.qa', access='stream',shared)
       lh=0
       if (myid==0) then
        write(9,pos=4*lh+1) mbk+1; lh=lh+1 ! Number of zones
        do mm = 0, mbk
           write(9,pos=4*lh+1) int4(lximb(mm)+1); lh=lh+1 ! IMax
           write(9,pos=4*lh+1) int4(letmb(mm)+1); lh=lh+1 ! JMax
           write(9,pos=4*lh+1) int4(lzemb(mm)+1); lh=lh+1 ! KMax
        end do
        lhmb(mb)=lh
        do mm = 0, mbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do
       end if
        call MPI_BCAST(lhmb,mbk+1,MPI_INTEGER,0,icom,ierr)
        if (myid==mo(mb)) then
           lh=lhmb(mb)
           write(9,pos=4*lh+1) real(amachoo,kind=4); lh=lh+1 ! Mach Number
           write(9,pos=4*lh+1) real(aoa,kind=4); lh=lh+1  
           write(9,pos=4*lh+1) real(reoo,kind=4); lh=lh+1 ! Reynolds Number
           write(9,pos=4*lh+1) real(times(ndata),kind=4); lh=lh+1 ! Time
        end if
        lp=lpos(myid)+lhmb(mb)+4
        ns=1; ne=5
        do n=ns,ne; lq=(n-ns)*ltomb
           selectcase(n)
           case(1,5); varr(:)=qb(:,n)
           case(2,3,4); varr(:)=qb(:,n)!*qa(:,1)
           end select
        do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
          write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
        end do; end do
        end do
        close(9)
        CALL MPI_BARRIER(icom,ierr)
        if (myid==0) then
           write(*,"('RMS Solution written!')") 
        end if
 end subroutine p3dwrms
!====================================================================================
!=====WRITE FUNCTION FILE PLOT3D
!====================================================================================
 subroutine wffile(fname,nout,ndim)
 integer, intent(in) :: nout,ndim
 character(16), intent(in) :: fname
 character(3) :: cout
 integer :: n,l,i
 
   write(cout,"(i3)") nout
   do i = 0, 2
   l=scan(cout,' ')
   if (l==0) exit
   cout(l:l)='0'
   end do

   if (nout==(ndata+1)) then
      cout='AVG'
   elseif (nout==(ndata+2)) then
      cout='RMS'
   end if

   if (myid==0) then
     open(9,file='out/'//trim(fname)//cout//'.f'); close(9,status='delete')
   end if
   CALL MPI_BARRIER(icom,ierr)
   open (unit=9, file='out/'//trim(fname)//cout//'.f', access='stream',shared)
   lh=0
   if (myid==0) then
    write(9,pos=4*lh+1) mbk+1; lh=lh+1 ! Number of zones
    do mm = 0, mbk
       write(9,pos=4*lh+1) int4(lximb(mm)+1); lh=lh+1 ! IMax
       write(9,pos=4*lh+1) int4(letmb(mm)+1); lh=lh+1 ! JMax
       write(9,pos=4*lh+1) int4(lzemb(mm)+1); lh=lh+1 ! KMax
       write(9,pos=4*lh+1) int4(ndim); lh=lh+1 ! #dimensions
    end do
    lhmb(mb)=lh
    do mm = 0, mbk-1
       lhmb(mm+1)=lhmb(mm)+ndim*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
    end do
   end if
   call MPI_BCAST(lhmb,mbk+1,MPI_INTEGER,0,icom,ierr)
   lp=lpos(myid)+lhmb(mb)
   ns=1; ne=ndim
   do n=ns,ne; lq=(n-ns)*ltomb
   varr(:)=qo(:,n)
   do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
     write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
   end do; end do
   end do
   close(9)
   CALL MPI_BARRIER(icom,ierr)
   if (myid==0) then
      write(*,"(a,' funtion written!')") trim(fname)//cout
   end if
 end subroutine wffile

!====================================================================================
!=====  PLOT3D Q FILES WRITE POST
!====================================================================================
  subroutine wrP3dP(nout,mblkin,cname)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nout
     character(*), intent(in),optional :: cname
     character(*),parameter :: fname='sol'
     character(:),allocatable :: lfname
     character(3) :: cout,ncout
     character(8) :: ctime
     character(len=*),parameter :: cext='.qa',cpath='out/'
     integer :: n,l,i,lh,iolen,foper,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

        ! rpt- Set default option to Multiblock
        if(present(mblkin)) then
           mblk=mblkin
        else
           mblk=1
        end if

        selectcase(mblk);
        case(1)
           cout=''
           foper=0
           wrcom=icom
           nbk=mbk
        case(0)
           write(cout,"(a,i2)") 'b',mb
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
           foper=mo(mb)
           wrcom=bcom
           nbk=0
        case default
           if(myid==0) write(*,*) 'Wrong multiblock option! Aborting...'
           CALL MPI_ABORT(icom,err,ierr)
        end select


        ! rpt- Set default option to Multiblock
        if(present(cname)) then
           write(ncout,"(i3)") nout
           do i = 0, 2
              l=scan(ncout,' ')
              if (l==0) exit
              ncout(l:l)='0'
           end do
           ctime=trim(cname)//ncout
        else
           if (nout.le.ndata) then
              write(ctime,"(f8.4)") times(nout)
              do i = 0, 8
                 l=scan(ctime,' ')
                 if (l==0) exit
                 ctime(l:l)='0'
              end do
           else if (nout==ndata+1) then
              ctime='A'
           else if (nout==ndata+2) then
              ctime='RMS'
           end if
        end if

        l=len(cpath)+len(fname)+len(trim(adjustl(ctime)))+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(adjustl(ctime))//trim(cout)//cext
        if(myid==foper) CALL MPI_FILE_DELETE(lfname,info,ierr)

        wrlen=5*(lmx+1)
        if(.not.allocated(q4)) allocate(q4(0:lmx,5))

           q4(:,:)=qo(:,:)

     if (fflag) then
        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.q4flag) then
           gsizes(:)=(/mbijkl(:),5/)
           lsizes(:)=(/mpijkl(:),5/)
           starts(:)=(/mpijks(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                           MPI_ORDER_FORTRAN,MPI_REAL4,q4arr,ierr) 
           CALL MPI_TYPE_COMMIT(q4arr,ierr)
           q4flag=.true.
        end if

        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,q4fh,ierr)

       lh=0
        if (myid==foper) then
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            ibuf=lximb(mm)+1; offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=letmb(mm)+1; offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=lzemb(mm)+1; offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do
        if (myid==mo(mb)) then
           lh=lhmb(mb)
            rbuf=amachoo; offset=lh*iolen ! Mach Number
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=aoa    ; offset=lh*iolen ! AoA
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=reoo   ; offset=lh*iolen ! Reynolds Number
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=timo   ; offset=lh*iolen ! Time
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
        end if
        disp=(lhmb(mb)+4)*iolen
        CALL MPI_FILE_SET_VIEW(q4fh,disp,MPI_REAL4,q4arr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL(q4fh,q4,wrlen,MPI_REAL4,ista,ierr)
        if (myid==foper) then
           write(*,"('Solution written! T= ',a)") trim(ctime)//trim(cout)
        end if
        CALL MPI_FILE_CLOSE(q4fh,ierr)

     end if

  end subroutine wrP3dP
!====================================================================================
!===== CONSTRUCT THE COEFFICIENTS ARRAY FOR INTEGRATION
!====================================================================================
 subroutine integCoef()

    if (myid==0) then
       write(*,"('Total amout of data: ',i3)") ndata
    end if
       ns=0; ne=ndata; if(.not.allocated(delt)) allocate(delt(ns:ne))
       fctr=half/(times(ne)-times(ns))
       delt(ns)=fctr*(times(ns+1)-times(ns)); 
       delt(ne)=fctr*(times(ne)-times(ne-1))
    do n=ns+1,ne-1
       delt(n)=fctr*(times(n+1)-times(n-1))
    end do
    
 end subroutine integCoef
!====================================================================================
!=====WRITE FUNCTION FILE PLOT3D
!====================================================================================
 subroutine getCp(n)
    integer(k4), intent(in) :: n

    call fillqo(n)
    ra0=two/(amachoo**2)
    qo(:,1)=(p(:)-poo)*ra0

 end subroutine getCp
!*****

end module rptpost
