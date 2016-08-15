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

 real(k4), dimension(:,:,:), allocatable :: myarr
contains

!====================================================================================
!=====PROBLEM SETUP
!====================================================================================
 subroutine setup
!===== INPUT PARAMETERS

    open(9,file='inputo.dat')
    read(9,*) cinput,mbk,bkx,bky,bkz
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
    read(9,*) cinput,fintg,rdis,xpos,ypos,atk
    close(9)

    cinput=cinput; fltk=pi*fltk; fltkbc=pi*fltkbc
    rhooo=one; poo=one/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1*amach1+amach2*amach2+amach3*amach3)
    srefoo=111/tempoo; srefp1dre=(srefoo+one)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=one/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo

    abc(:,0)=(/a01,a02,a03,a04,a05,a06/)
    abc(:,1)=(/a10,a12,a13,a14,a15,a16/)
    abc(:,2)=(/a20,a21,a23,a24,a25,a26/)

    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),npc(0:mbk,3))
    allocate(lxibk(0:bkx-1),letbk(0:bky-1),lzebk(0:bkz-1))

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
    m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)
 
     rr(:,1)=ss(:,2)
    m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)
 
     rr(:,1)=ss(:,3)
    m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
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
    yaco(:)=three/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
               +qa(:,1)*xim(:,2)+qa(:,2)*etm(:,2)+qa(:,3)*zem(:,2)&
               +de(:,1)*xim(:,3)+de(:,2)*etm(:,3)+de(:,3)*zem(:,3))

    call wallArea
    call walldir
    if (fintg) then
       call intgUp(rdis,xpos,ypos,atk)
    end if
    !nprob=31
    !call probUp(nprob,(/-0.5_k8,0.15_k8/),(/0.48_k8,0.15_k8/),0.01_k8)

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
     ii=min(4,mbci)
     cbca(3,ii)=albed(2,2,0)
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

        de(:,1:3)=0

        rr(:,1)=qo(:,2)
        m=1; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
        de(:,2)=de(:,2)+rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)
        de(:,3)=de(:,3)-rr(:,1)*xim(:,2)-rr(:,2)*etm(:,2)-rr(:,3)*zem(:,2)

        rr(:,1)=qo(:,3)
        m=2; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
        de(:,3)=de(:,3)+rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
        de(:,1)=de(:,1)-rr(:,1)*xim(:,3)-rr(:,2)*etm(:,3)-rr(:,3)*zem(:,3)

        rr(:,1)=qo(:,4)
        m=3; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
        de(:,1)=de(:,1)+rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
        de(:,2)=de(:,2)-rr(:,1)*xim(:,1)-rr(:,2)*etm(:,1)-rr(:,3)*zem(:,1)

        qo(:,2)=de(:,1)*yaco(:);
        qo(:,3)=de(:,2)*yaco(:);
        qo(:,4)=de(:,3)*yaco(:)

 end subroutine getCurl
!====================================================================================
!===== SUBROUTINE FOR CALCULATING VORTICITY
!====================================================================================
 subroutine getVGrad(nvar)

        implicit none
        integer, intent(in) :: nvar

        if(allocated(fout)) deallocate(fout)
        allocate(fout(0:lmx,9))

        rr(:,1)=qo(:,2)
        m=1; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
        fout(:,1)=rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
        fout(:,2)=rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
        fout(:,3)=rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)

        rr(:,1)=qo(:,3)
        m=2; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
        fout(:,4)=rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
        fout(:,5)=rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
        fout(:,6)=rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)

        rr(:,1)=qo(:,4)
        m=3; call mpigo(ntdrv,nrone,n45no,m);
        call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
        fout(:,7)=rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
        fout(:,8)=rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
        fout(:,9)=rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)

        do m = 1, 9
           fout(:,m)=fout(:,m)*yaco(:)
        end do

 end subroutine getVGrad
!====================================================================================
!=====COMPUTE Q-CRITERION
!====================================================================================
 subroutine qcriterion(nvar)
    implicit none
    integer, intent(in) :: nvar

    call rdP3dS(nvar,fmblk)
    !p(:)=qo(:,5)
    if (nviscous==1) then
     de(:,2)=qo(:,2)
     de(:,3)=qo(:,3)
     de(:,4)=qo(:,4)
 
     rr(:,1)=de(:,2)
     m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     txx(:)=yaco(:)*(xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3))
     hzz(:)=yaco(:)*(xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3))
     tzx(:)=yaco(:)*(xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3))
 
     rr(:,1)=de(:,3)
     m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
     txy(:)=yaco(:)*(xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3))
     tyy(:)=yaco(:)*(xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3))
     hxx(:)=yaco(:)*(xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3))
 
     rr(:,1)=de(:,4)
     m=4; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
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
!=====COMPUTE dui/dxi+Q+VORTICITY
!====================================================================================
 subroutine getAllDs(nvar)
    implicit none
    integer, intent(in) :: nvar
    real(k4),parameter :: r108=-1_k4/108_k4,r27=-1_k4/27_k4,&
                          r4=-1_k4/4_k4,r18=-1_k4/18_k4

    call rdP3dS(nvar,fmblk)
    p(:)=qo(:,5)
    de(:,2)=qo(:,1)
    de(:,1)=1/de(:,2)
    de(:,5)=gam*p(:)*de(:,1)
    de(:,1)=srefp1dre*de(:,5)**1.5_k8/(de(:,5)+srefoo)
    call getVGrad(nvar)

     ! Q = 2nd Invariant of the characteristic eq
     qo(:,1)=2*(half*(fout(:,2)-fout(:,4)))**2 + 2*(half*(fout(:,3)-fout(:,9)))**2 + &
     2*(half*(fout(:,6)-fout(:,8)))**2 - fout(:,1)**2 - fout(:,5)**2 - fout(:,9)**2 -   &
     2*(half*(fout(:,2)+fout(:,4)))**2 - 2*(half*(fout(:,3)+fout(:,7)))**2 -         &
     2*(half*(fout(:,6)+fout(:,8)))**2

     qo(:,2)=fout(:,8)-fout(:,6)
     qo(:,3)=fout(:,3)-fout(:,7)
     qo(:,4)=fout(:,4)-fout(:,2)

     ! P = 1st Invariant of the characteristic eq
     de(:,3)=fout(:,1)+fout(:,5)+fout(:,9)
     ! R = 3rd Invariant of the characteristic eq
     de(:,4)=fout(:,1)*qo(:,2)+fout(:,5)*qo(:,3)+fout(:,9)*qo(:,4)
     ! PQ
     de(:,5)=de(:,3)*qo(:,1)

     ! Compute discriminant DELTA 
     ! (PQ)**2/108
     qo(:,5)=r108*de(:,5)
     ! -(Q**3+P**3R)/27
     qo(:,5)=qo(:,5)-r27*(qo(:,1)*qo(:,1)*qo(:,1)+de(:,3)*de(:,3)*de(:,3)*de(:,4))
     ! -R**2/4
     qo(:,5)=qo(:,5)-r4*de(:,4)*de(:,4)
     ! +(PQR)/6
     qo(:,5)=qo(:,5)+de(:,5)*de(:,4)

    if (wflag) then
    ! READ VARIABLES
    if (nviscous==1) then
 
     fctr=2.0_k8/3
     rr(:,1)=-de(:,1)
     rr(:,2)=1/yaco(:)
     de(:,5)=fctr*(fout(:,1)+fout(:,5)+fout(:,9))*rr(:,2)
 
     txx(:)=rr(:,1)*(2*fout(:,1)-de(:,5))
     tyy(:)=rr(:,1)*(2*fout(:,5)-de(:,5))
     tzz(:)=rr(:,1)*(2*fout(:,9)-de(:,5))
     txy(:)=rr(:,1)*(fout(:,4)+fout(:,2))
     tyz(:)=rr(:,1)*(fout(:,8)+fout(:,6))
     tzx(:)=rr(:,1)*(fout(:,3)+fout(:,7))
 
       if(.not.allocated(tw)) allocate(tw(0:lcwall,3))
       do ll = 0, lcwall; l=lwall(ll)
         tw(ll,1)=de(l,2)*(txx(l)*wnor(ll,1)+txy(l)*wnor(ll,2)+tzx(l)*wnor(ll,3))
         tw(ll,2)=de(l,2)*(txy(l)*wnor(ll,1)+tyy(l)*wnor(ll,2)+tyz(l)*wnor(ll,3))
         tw(ll,3)=de(l,2)*(tzx(l)*wnor(ll,1)+tyz(l)*wnor(ll,2)+tzz(l)*wnor(ll,3))
       end do
    end if
    end if

 end subroutine getAllDs
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
      if (myid==foper) then
         open(90,file='out/filelist'//rcout,shared)
         lmt=-2
         do while (stat.ge.0)
            read(90,*,IOSTAT=stat) 
            lmt=lmt+1
         end do
         close(90)
         ndata=lmt
         l=len('out/sotT')+8+len(trim(rcout))+len('.q')
      end if

      CALL MPI_BCAST(ndata,1,MPI_INTEGER4,foper,wrcom,ierr)
      CALL MPI_BCAST(l,1,MPI_INTEGER4,foper,wrcom,ierr)
      allocate(character(l) :: ofiles(0:ndata))
      allocate(character(l) :: ofile)
      allocate(times(0:ndata))

      if (myid==foper) then
         open(90,file='out/filelist'//rcout)
         do i = 0, ndata
            read(90,"(a)") ofile
            ofiles(i)=ofile
            ctime=(ofile(9:16)//'e0')
            read(ctime,*) times(i)
         end do
         close(90)
      end if
      CALL MPI_BCAST(times,ndata+1,MPI_REAL8,foper,wrcom,ierr)
      CALL MPI_BCAST(ofiles,l*(ndata+1),MPI_CHAR,foper,wrcom,ierr)
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
 subroutine p3dStats
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
       if(.not.allocated(fout)) allocate(fout(0:lmx,6))
       fout(:,:)=0
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
          qa(:,:)=zero
          qb(:,:)=zero
       nn=ndata*0.01_k8
       wts=0
       do n=0,ndata
       if (myid==mo(mb)) then
          if (mod(n,nn)==0) then
             wte=MPI_WTIME(); res=wte-wts
             write(*,"('Block ',i2,x,f5.1,'% Averaged',a,' took ',f5.1,'s')")&
             mb,real(n)*100.0e0/real(ndata),cout,res
             wts=MPI_WTIME()
          end if
       end if
          call rdP3dS(n,fmblk)
          de(:,:)=delt(n)*qo(:,:)
          qa(:,:)=qa(:,:)+de(:,:)
          qb(:,:)=qb(:,:)+de(:,:)*qo(:,:)
          fout(:,4)=fout(:,4)+de(:,2)*qo(:,3)
          fout(:,5)=fout(:,5)+de(:,2)*qo(:,4)
          fout(:,6)=fout(:,6)+de(:,3)*qo(:,4)
       end do
       do i = 2, 4
          fout(:,i-1)=qb(:,i)-qa(:,i)*qa(:,i)
       end do
       fout(:,4)=fout(:,4)-qa(:,2)*qa(:,3)
       fout(:,5)=fout(:,5)-qa(:,2)*qa(:,4)
       fout(:,6)=fout(:,6)-qa(:,3)*qa(:,4)
       !!===== WRITE AVERAGE VALUES 
           qo(:,:)=qa(:,:)
           call wrP3dP(ndata+1,fmblk)
       !!===== WRITE RMS VALUES 
          qo(:,:)=sqrt(qb(:,:))
          call wrP3dP(ndata+2,fmblk)
          call wrP3dF('Rij',0,6,fmblk)
    end if
 end subroutine p3dStats


!====================================================================================
!=====  PLOT3D F FILES Read
!====================================================================================
  subroutine rdP3dF(fname,nout,ndim,mblkin)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nout,ndim
     character(len=*), intent(in) :: fname
     character(len=:),allocatable :: lfname
     character(3) :: cout,cblk
     character(len=*),parameter :: cext='.f',cpath='out/'
     integer :: n,l,i,lh,iolen,foper,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: fh,amode,farr
     integer(k4) :: ibuf
     integer, dimension (4) :: gsizes,lsizes,starts

     ! rpt- Set default option to Multiblock
     if(present(mblkin)) then
        mblk=mblkin
     else
        mblk=1
     end if

     selectcase(mblk);
     case(1)
        cblk=''
        foper=0
        wrcom=icom
        nbk=mbk
     case(0)
        write(cblk,"(i2,a)") mb,'n'
        do i = 0, 1
           l=scan(cblk,' ')
           if (l==0) exit
           cblk(l:l)='0'
        end do
        foper=mo(mb)
        wrcom=bcom
        nbk=0
     case default
        if(myid==0) write(*,*) 'Wrong multiblock option! Aborting...'
        CALL MPI_ABORT(icom,err,ierr)
     end select
     write(cout,"(i3)") nout
     do i = 0, 2
        l=scan(cout,' ')
        if (l==0) exit
        cout(l:l)='0'
     end do
     l=len(cpath)+len(fname)+len(trim(cblk))+len(cout)+len(cext)
     allocate(character(len=l) :: lfname)
     lfname=cpath//trim(fname)//trim(cblk)//cout//cext

     if (allocated(fout)) deallocate(fout)
     if (.not.allocated(fout)) allocate(fout(0:lmx,ndim))
     wrlen=ndim*(lmx+1)
     amode=MPI_MODE_RDONLY

     CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
     gsizes(:)=(/mbijkl(:),ndim/)
     lsizes(:)=(/mpijkl(:),ndim/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL4,farr,ierr) 
     CALL MPI_TYPE_COMMIT(farr,ierr)
     
     CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,fh,ierr)
 
     lh=0
     if (myid==foper) then
      ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
      CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
      do l = 0, nbk
         mm=l+(1-mblk)*mb
         ibuf=lximb(mm)+1; offset=lh*iolen ! IMax
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         ibuf=letmb(mm)+1; offset=lh*iolen ! JMax
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         ibuf=lzemb(mm)+1; offset=lh*iolen ! KMax
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         ibuf=ndim; offset=lh*iolen        ! #Dimensions
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
      end do
     end if
     l=(1-mblk)*mb
     lhmb(l)=1+(nbk+1)*4
     do mm = 0, nbk-1
        lhmb(mm+1)=lhmb(mm)+ndim*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
     end do
     disp=lhmb(mb)*iolen
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,farr,'native',info,ierr)
     CALL MPI_FILE_READ_ALL(fh,fout,wrlen,MPI_REAL4,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(farr,ierr)
     if (myid==foper) then
        write(*,"(a,' funtion read!')") trim(fname)//cout
     end if
  end subroutine rdP3dF

!===================================================================================
!=====  PLOT3D Q FILES READ POST AT POSITION
!===================================================================================
  subroutine rdP3dPat(nout,blk,nxk,pos,mblkin,cname)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nout,blk
     integer, dimension(2), intent(in) :: nxk
     integer, dimension(3,nxk(1),nxk(2)), intent(in) :: pos
     character(*), intent(in),optional :: cname
     character(*),parameter :: fname='solT'
     character(:),allocatable :: lfname
     character(3) :: cout,ncout
     character(8) :: ctime
     character(:), allocatable :: cext
     character(len=*),parameter :: cext0='.q',cext1='.qa',cpath='out/'
     integer :: n,l,i,k,mm,m,lh,iolen,nbk,err,mblk
     integer :: disp
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
           nbk=mbk
        case(0)
           write(cout,"(a,i2)") 'b',blk
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
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
           allocate(character(len=len(cext1)) :: cext)
           cext=cext1
        else
           if (nout.le.ndata) then
              write(ctime,"(f8.4)") times(nout)
              do i = 0, 8
                 l=scan(ctime,' ')
                 if (l==0) exit
                 ctime(l:l)='0'
              end do
              allocate(character(len=len(cext0)) :: cext)
              cext=cext0
           else if (nout==ndata+1) then
              ctime='A'
              allocate(character(len=len(cext1)) :: cext)
              cext=cext1
           else if (nout==ndata+2) then
              ctime='RMS'
              allocate(character(len=len(cext1)) :: cext)
              cext=cext1
           end if
        end if

        l=len(cpath)+len(fname)+len(trim(adjustl(ctime)))+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(adjustl(ctime))//trim(cout)//cext

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)

        l=(1-mblk)*nbk
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do

        open(91,file=trim(lfname),access='stream',form='unformatted') 
        if(.not.allocated(mval)) allocate(mval(nxk(1),nxk(2),5,0:ndata+2))
        lq=(lximb(blk)+1)*(letmb(blk)+1)*(lzemb(blk)+1)
        do m = 0, 4
           disp=(lhmb(blk)+4+lq*m)*iolen
           do k = 1, nxk(2)
              do i = 1, nxk(1)
                 iopos=pos(3,i,k)*(lximb(blk)+1)*(letmb(blk)+1)
                 iopos=iopos+pos(2,i,k)*(lximb(blk)+1)
                 iopos=(iopos+pos(1,i,k))*iolen
                 read(91,pos=disp+iopos+1) mval(i,k,m+1,nout) 
              end do
           end do
        end do

  end subroutine rdP3dPat
!===================================================================================
!=====  PLOT3D Q FILES READ POST
!===================================================================================
  subroutine rdP3dP(nout,mblkin,cname)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nout
     character(*), intent(in),optional :: cname
     character(*),parameter :: fname='solT'
     character(:),allocatable :: lfname
     character(3) :: cout,ncout
     character(8) :: ctime
     character(:), allocatable :: cext
     character(len=*),parameter :: cext0='.q',cext1='.qa',cpath='out/'
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
           allocate(character(len=len(cext1)) :: cext)
           cext=cext1
        else
           if (nout.le.ndata) then
              write(ctime,"(f8.4)") times(nout)
              do i = 0, 8
                 l=scan(ctime,' ')
                 if (l==0) exit
                 ctime(l:l)='0'
              end do
              allocate(character(len=len(cext0)) :: cext)
              cext=cext0
           else if (nout==ndata+1) then
              ctime='A'
              allocate(character(len=len(cext1)) :: cext)
              cext=cext1
           else if (nout==ndata+2) then
              ctime='RMS'
              allocate(character(len=len(cext1)) :: cext)
              cext=cext1
           end if
        end if

        l=len(cpath)+len(fname)+len(trim(adjustl(ctime)))+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(adjustl(ctime))//trim(cout)//cext

        wrlen=5*(lmx+1)
        if(.not.allocated(q4)) allocate(q4(0:lmx,5))

        amode=MPI_MODE_RDONLY

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

        if (fflag) then
           CALL MPI_FILE_OPEN(wrcom,trim(lfname) ,amode ,info ,q4fh,ierr)

           l=(1-mblk)*mb
           lhmb(l)=1+(nbk+1)*3
           do mm = 0, nbk-1
              lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
           end do
           lh=lhmb(mb)
           offset=lh*iolen ! Mach Number
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           amachoo=rbuf;
           offset=lh*iolen ! AoA
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           aoa=rbuf;
           offset=lh*iolen ! Reynolds Number
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           reoo=rbuf;
           offset=lh*iolen ! Time
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           timo=rbuf;

           disp=(lhmb(mb)+4)*iolen
           CALL MPI_FILE_SET_VIEW(q4fh,disp,MPI_REAL4,q4arr,'native',info,ierr)
           CALL MPI_FILE_READ_ALL(q4fh,q4,wrlen,MPI_REAL4,ista,ierr)
           if (myid==foper) then
              if (nout==ndata+1) then
                 write(*,"('AVG Solution read!')") 
              else if(nout==ndata+2) then
                 write(*,"('RMS Solution read!')") 
              else
                 if (present(cname)) then
                    write(*,*) lfname//" Read!"
                 else
                    write(*,"('Solution read! T= ',f8.4)") times(nout)
                 end if
              end if
           end if

           CALL MPI_FILE_CLOSE(q4fh,ierr)
           qo(:,:)=q4(:,:)
        else
           qo(:,:)=0
        end if
  end subroutine rdP3dP

!===================================================================================
!=====  PLOT3D Q FILES WRITE POST
!===================================================================================
  subroutine wrP3dP(nout,mblkin,cname)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nout
     character(*), intent(in),optional :: cname
     character(*),parameter :: fname='solT'
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
!===================================================================================
!=====COMPUTE CP AND STORE IT IN QO(:,1)
!===================================================================================
 subroutine getCp(n)
    integer(k4), intent(in) :: n

    if (fflag) then
       call fillqo(n)
       ra0=two/(amachoo**2)
       qo(:,1)=(p(:)-poo)*ra0
    end if

 end subroutine getCp

!===================================================================================
!=====SET UP INTEGRAL PARAMETERS
!===================================================================================
 subroutine intgUp(rdis,xpos,ypos,k)
 implicit none
    real, intent(in) :: rdis,xpos,ypos
    integer , intent(in) :: k
    integer(k4) :: color,i,j,l,ll
    real(k8) :: g11, g22, g12,rdis2,tmp

    ! Initialise values
    intgflag=.False.
    g11=0;g22=g11;g12=g22
    color=MPI_UNDEFINED

    rdis2=rdis**2
    ll=-1
    do kk = 35, 39
    do j = 45, 60
       do i = 52, 58; l=indx3(i,j,kk,1)
    !do kk = 0, lze
    !do j = 0, let
    !   do i = 0, lxi; l=indx3(i,j,kk,1)
    !      tmp=(xyz(l,1)-xpos)**2+(xyz(l,2)-ypos)**2
          if (mb==7) then
          !if (tmp-rdis2<0) then
             ll=ll+1; de(ll,5)=l+sml
          end if
       end do
    end do
    end do
    lcintg=ll
    if (lcintg.ne.-1) then
       color=1
       intgflag=.True.
       allocate(lintg(0:lcintg),aintg(0:lcintg),vintg(0:lcintg))
       do ll = 0, lcintg; l=de(ll,5); lintg(ll)=l
          g11 = qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
          g22 = qo(l,2)*qo(l,2)+qa(l,2)*qa(l,2)+de(l,2)*de(l,2)
          g12 = qo(l,1)*qo(l,2)+qa(l,1)*qa(l,2)+de(l,1)*de(l,2)
          aintg(ll)=sqrt(g11*g22-g12*g12)
       end do
    end if

    if (myid==0) then
       color=1
    end if
    call MPI_COMM_SPLIT(icom,color,myid,intgcom,ierr)

 end subroutine intgUp
!===================================================================================
!=====PERFORM INTEGRAL
!===================================================================================
subroutine integrate()
implicit none
   integer(k4) :: l,ll
   
   if (intgflag) then
      do ll = 0, lcintg; l=lintg(ll)
         vintg(ll)=aintg(ll)*varr(l)
      end do
      ra0=sum(vintg)/sum(aintg)
   else 
      ra0=0
   end if

   if (fmblk==1) then
      if (intgflag.or.myid==0) then
        CALL MPI_REDUCE(ra0,ra1,1,MPI_REAL8,MPI_SUM,0,intgcom,ierr)
      end if
      if (myid==0) ra0=ra1
   end if

end subroutine integrate

!===================================================================================
!=====SET UP PROBE CYLINDERS
!===================================================================================
subroutine probUp(nprob,sprob,eprob,orprob)
implicit none
   integer(k4),intent(in) :: nprob
   real(k8),intent(in), optional :: orprob
   real(k8),dimension (2),intent(in) :: sprob,eprob
   real(k8),dimension (2)            :: dirprob
   integer :: color,err
   real(k8) :: rprob,r2prob,g11,g22,g12

   allocate(xyprob(2,nprob),nklprob(2,0:lze,nprob),nlprob(2,nprob))
   allocate(mprob(nprob),probcom(nprob),probflag(nprob))

   dirprob=eprob-sprob
   if (nprob>1) then
      ra0=sqrt(dirprob(1)**2+dirprob(2)**2)/(nprob-1)
      rprob=ra0*half
      dirprob(:)=dirprob(:)/(nprob-1)
      if (present(orprob)) then
         rprob=orprob
         if(myid==0) write(*,*) rprob
      end if
   else 
      if (present(orprob)) then
         rprob=orprob
         dirprob(:)=0
      else
         if (myid==0) then
            write(*,"('Must specify probe radius for single probbing!!')") 
            write(*,"('Aborting..')") 
         end if
         CALL MPI_ABORT(icom,err,ierr)
      end if
   end if

   do i = 0, nprob-1
      xyprob(:,i+1)=sprob(:)+i*dirprob(:)
   end do

   if (rprob<0) then
   ! in development
   else
      r2prob=rprob**2;ll=0
      do n = 1, nprob
         nlprob(1,n)=ll
         do k = 0, lze
         nklprob(1,k,n)=ll
         do j = 0, let
            do i = 0, lxi; l=indx3(i,j,k,1)
               tmp=(xyz(l,1)-xyprob(1,n))**2+(xyz(l,2)-xyprob(2,n))**2
               if (tmp-r2prob<0) then
                  ll=ll+1; de(ll,5)=l+sml
               end if
            end do
         end do
         nklprob(2,k,n)=ll
         end do
         nlprob(2,n)=ll
      end do
      lcprob=ll
      if(lcprob>0) then
         allocate(lprob(0:lcprob),aprob(0:lcprob),&
                  vprob(0:lcprob))
         do ll = 0, lcprob; l=de(ll,5); lprob(ll)=l
            g11 = qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
            g22 = qo(l,2)*qo(l,2)+qa(l,2)*qa(l,2)+de(l,2)*de(l,2)
            g12 = qo(l,1)*qo(l,2)+qa(l,1)*qa(l,2)+de(l,1)*de(l,2)
            aprob(ll)=sqrt(g11*g22-g12*g12)
         end do
      end if


      do n = 1, nprob
         ii=nlprob(2,n)-nlprob(1,n)
         if (ii.ne.0) then
            color=n
            probflag(n)=.true.
         else
            color=MPI_UNDEFINED
            probflag(n)=.false.
         end if
         call MPI_COMM_SPLIT(icom,color,myid,probcom(n),ierr)
         if (probflag(n)) then
            CALL MPI_ALLREDUCE(myid,mprob(n),1,MPI_INTEGER4,MPI_MIN,probcom(n),ierr)
            write(*,*)&
            myid,mprob(n),n,sum(aprob(nlprob(1,n):nlprob(2,n)))/(lze+1),nlprob(1,n),nlprob(2,n)
         end if
      end do


   end if
   
end subroutine probUp

subroutine probCirc(ntotal)
implicit none
   integer, intent(in) :: ntotal
   integer :: l,ll,n,m
   real(k8), dimension(:,:,:), allocatable :: probze
   character(3) :: cprob
   
 
   do m = 0, ntotal
      call rdP3dP(m,fmblk)
      varr(:)=qo(:,5)
      do n = 1, nprob
         if (probflag(n)) then
         probze(:,m,n)=0
         if(.not.allocated(probze)) allocate(probze(0:lze,0:ntotal,nprob))
            do k = 0, lze
            ra0=0
               do ll = nklprob(1,k,n), nklprob(2,k,n); l=lprob(ll)
                  ra0=ra0+aprob(ll)*varr(l)
               end do
               ra1=sum(aprob(nklprob(1,k,n):nklprob(2,k,n)))
               ra0=ra0/ra1
               CALL MPI_REDUCE(ra0,probze(k,m,n),1,MPI_REAL8,MPI_SUM,0,probcom(n),ierr)
            end do
         else
            ra0=0
         end if
      end do
   end do

   do n = 1, nprob
      write(cprob,"(i3)") n
      do i = 0, 3
         l=scan(cprob,' ')
         if (l==0) exit
         cprob(l:l)='0'
      end do
      if(myid==mprob(n).and.probflag(n)) then
         write(*,"('Processor ',i2,' writing probe ',i2)") myid,n
         open(unit=50+n, file='out/circ'//cprob//'.dat')
         write(50+n,"('t* circ (',f8.4,',',f8.4,')')") xyprob(1,n),xyprob(2,n)
         do m = 0, ntotal
            write(50+n,"(202(es15.7))") times(m),probze(:,m,n)
         end do
         close(50+n)
         write(*,"('Finished writing probe ',i2)") n
      end if
   end do
   CALL MPI_BARRIER(icom,ierr)

end subroutine probCirc

subroutine getijkMax(blk,nxk,xs,ks)
implicit none
integer, intent(in) :: blk
integer, dimension(2), intent(in) :: nxk
integer(k4) , dimension(nxk(1)), intent (in) :: xs
integer(k4) , dimension(nxk(2)), intent (in) :: ks
integer(k4) :: pos,mymaxid,flg,flg2
integer(k4) , dimension(nxk(1),nxk(2)) :: maxflag,maxcom
integer(k4) , dimension(nxk(2)) :: maxkcom
real(k4), dimension(0:npc(blk,2)-1) :: val
real(k4), dimension(3) :: sdval
integer(k4), dimension(3) :: sdpos
real(k4), dimension(0:let) :: jvarr
real(k4) :: rsp0,rsp1

maxflag=0
color=0
allocate(maxpos(3,nxk(1),nxk(2)),maxxyz(3,nxk(1),nxk(2)))

   if (mb==blk) then
      do k = 1, nxk(2)
         if ((mpijks(3).le.ks(k)).and.(mpijke(3).ge.ks(k))) then
            do i = 1, nxk(1)
               if ((mpijks(1).le.xs(i)).and.(mpijke(1).ge.xs(i))) then
                  maxflag(i,k)=1
               end if
            end do
         end if
      end do
      do k = 1, nxk(2)
        do i = 1, nxk(1)
           color=maxflag(i,k)
           call MPI_COMM_SPLIT(bcom,color,myid,maxcom(i,k),ierr)
        end do
      end do
      do k = 1, nxk(2);kk=ks(k)-mpijks(3)
        do i = 1, nxk(1);ii=xs(i)-mpijks(1)
           if (maxflag(i,k)==1) then
              do j = 0, let;ll=indx3(ii,j,kk,1)
                 jvarr(j)=varr(ll)
              end do
              rsp0=maxval(jvarr(:))
              pos=maxloc(jvarr(:),1)-1
              CALL MPI_ALLGATHER(rsp0,1,MPI_REAL4,val,1,MPI_REAL4,maxcom(i,k),ierr)
              rsp1=maxval(val)
              if (rsp0==rsp1) then
                 maxpos(1,i,k)=xs(i)
                 maxpos(2,i,k)=pos+mpijks(2)
                 maxpos(3,i,k)=ks(k)
              else 
                 maxflag(i,k)=0
              end if
           end if
        end do
      end do
      do k = 1, nxk(2)
         do i = 1, nxk(1)
            call MPI_COMM_FREE(maxcom(i,k),ierr)
         end do
      end do
      do k = 1, nxk(2)
         l=sum(maxflag(:,k))
         if (l>0) then
            color=1
         else
            color=MPI_UNDEFINED
         end if
         call MPI_COMM_SPLIT(bcom,color,myid,maxkcom(k),ierr)
         if (color==1) then
            CALL MPI_COMM_RANK(maxkcom(k),mymaxid,ierr) 
            if (mymaxid==0) then
               do i = 1, nxk(1)
                  if (maxflag(i,k)==1) then
                      ii=maxpos(1,i,k)-mpijks(1)
                      jj=maxpos(2,i,k)-mpijks(2)
                      kk=maxpos(3,i,k)-mpijks(3)
                      l=indx3(ii,jj,kk,1)
                      maxxyz(:,i,k)=xyz(l,:)
                  else
                     CALL MPI_RECV(maxxyz(1,i,k),3,MPI_REAL4,MPI_ANY_SOURCE,i,&
                                   maxkcom(k),ista,ierr)
                     CALL MPI_RECV(maxpos(1,i,k),3,MPI_INTEGER,MPI_ANY_SOURCE,i,&
                                   maxkcom(k),ista,ierr)
                  end if
                  flg=nxk(1)*(k-1)+i
                  flg2=flg+(nxk(1)*nxk(2))
                  CALL MPI_SEND(maxpos(1,i,k),3,MPI_INTEGER,0,flg,icom,ierr)
                  CALL MPI_SEND(maxxyz(1,i,k),3,MPI_REAL4,0,flg2,icom,ierr)
               end do
            else
               do i = 1, nxk(1)
                  if (maxflag(i,k)==1) then
                      ii=maxpos(1,i,k)-mpijks(1)
                      jj=maxpos(2,i,k)-mpijks(2)
                      kk=maxpos(3,i,k)-mpijks(3)
                      l=indx3(ii,jj,kk,1)
                      sdval(:)=xyz(l,:)
                      sdpos(:)=maxpos(:,i,k)
                     CALL MPI_SEND(sdval,3,MPI_REAL4,0,i,maxkcom(k),ierr)
                     CALL MPI_SEND(sdpos,3,MPI_INTEGER,0,i,maxkcom(k),ierr)
                  end if
               end do
            end if
         end if
      end do
   end if

   if (myid==0) then
      do k = 1, nxk(2)
         write(*,"(a,x,i2,7x,a,x,i3)") 'nk =',k,'ze =',ks(k)
         write(*,"(a,7x,a,7x,a,7x,a,7x,a,7x,a)") 'xi','et','ze','x','y','z'
         do i = 1, nxk(1)
            flg=nxk(1)*(k-1)+i
            flg2=flg+(nxk(1)*nxk(2))
            CALL MPI_RECV(maxpos(1,i,k),3,MPI_INTEGER,MPI_ANY_SOURCE,flg,icom,ista,ierr)
            CALL MPI_RECV(maxxyz(1,i,k),3,MPI_REAL4,MPI_ANY_SOURCE,flg2,icom,ista,ierr)
            write(*,"(i3,7x,i3,7x,i3,7x,f7.3,7x,f7.3,7x,f7.3)")&
            maxpos(1,i,k),maxpos(2,i,k),maxpos(3,i,k),maxxyz(1,i,k),maxxyz(2,i,k),maxxyz(3,i,k)
         end do
      end do
   end if

   ll=3*nxk(1)*nxk(2)
   CALL MPI_BCAST(maxpos,ll,MPI_INTEGER,0,icom,ierr)
   CALL MPI_BCAST(maxxyz,ll,MPI_REAL4,0,icom,ierr)

end subroutine getijkMax

subroutine getvalMax(blk,nxk,nout)
implicit none
integer, intent(in) :: blk,nout
integer, dimension(2), intent(in) :: nxk
   
   if(.not.allocated(nose)) allocate(nose(0:npro-1,2))
   lq=nout
   mm=mod((lq+1),npro)
   do m = 0, npro-1
      nose(m,1)=m*((lq+1)/npro)+min(m,mm)
      nose(m,2)=nose(m,1)+(lq+1)/npro+min(1,mm/(m+1))-1
   end do
   if (myid==0) then
   write(*,*) mm
       do i = 0, npro-1
          write(*,"(i3,5x,i5,5x,i5,5x,i5)") i,nose(i,:),nose(i,2)-nose(i,1)+1
       end do
   end if
   lp=(nose(myid,2)-nose(myid,1))/((nose(myid,2)-nose(myid,1))*0.1e0)
   mm=0
   do m = nose(myid,1), nose(myid,2)
      mm=mm+1
      call rdP3dPat(m,blk,nxk,maxpos)
      if (mod(mm,lp)==0) then
         write(*,"(a,x,i4,x,a,x,i4,x,a,f7.3,x,a)") 'Process:',myid,'Snap',m,'@',times(m),'read!'
      end if
   end do
   
   if (myid==0) then
      do m = 1, npro-1
         l=nose(m,1);ll=(nose(m,2)-nose(m,1)+1)*nxk(1)*nxk(2)*5
         CALL MPI_RECV(mval(1,1,1,l),ll,MPI_REAL4,MPI_ANY_SOURCE,m,icom,ista,ierr)
      end do
   else
      l=nose(myid,1);ll=(nose(myid,2)-nose(myid,1)+1)*nxk(1)*nxk(2)*5
      CALL MPI_SEND(mval(1,1,1,l),ll,MPI_REAL4,0,myid,icom,ierr)
   end if
   
   ll=(ndata+1)*nxk(1)*nxk(2)*5
   CALL MPI_BCAST(mval,ll,MPI_REAL4,0,icom,ierr)
   do k = 1, nxk(2)
      if (myid==k-1) then
      write(cinput,"(i3)") k
      cstring='maxpln'//trim(adjustl(cinput))//'.dat'
      open(unit=100, file=trim(adjustl(cstring)))
      write(100,"(a,7x)",advance='no') 't'
      do i = 1, nxk(1)-1
         write(cinput,"(i3)") i
         write(100,"(a,7x,a,7x,a,7x,a,7x,a,7x)",advance='no')&
         'rho'//trim(adjustl(cinput)),'u'//trim(adjustl(cinput)),'v'//trim(adjustl(cinput)),&
         'w'//trim(adjustl(cinput)),'p'//trim(adjustl(cinput))
      end do
         write(cinput,"(i3)") nxk(1)
         write(100,"(a,7x,a,7x,a,7x,a,7x,a,7x)")&
         'rho'//trim(adjustl(cinput)),'u'//trim(adjustl(cinput)),'v'//trim(adjustl(cinput)),&
         'w'//trim(adjustl(cinput)),'p'//trim(adjustl(cinput))
      do m = 0,ndata
         write(100,"(es15.7,7x)",advance='no') times(m)
         do i = 1, nxk(1)-1
            write(100,"(es15.7,7x,es15.7,7x,es15.7,7x,es15.7,7x,es15.7,7x)",advance='no')&
            mval(i,k,:,m)
         end do
         write(100,"(es15.7,7x,es15.7,7x,es15.7,7x,es15.7,7x,es15.7,7x)")&
         mval(i,k,:,m)
      end do
      write(*,"(a,x,a)") trim(adjustl(cstring)),'written!'
      close(100)
      end if
   end do
end subroutine getvalMax
!*****

end module rptpost
