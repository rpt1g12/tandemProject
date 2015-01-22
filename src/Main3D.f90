!*****
!***** 3D PARALLEL SOLVER
!*****

 program main3d

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 use rpt
 implicit none

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

!===== INPUT PARAMETERS

    open(9,file='inputo.dat',shared)
    read(9,*) cinput,mbk
    read(9,*) cinput,nts
    read(9,*) cinput,nscrn,nsgnl
    read(9,*) cinput,ndata,ndatp
    read(9,*) cinput,nkrk
    read(9,*) cinput,nviscous
    read(9,*) cinput,nsmf
    read(9,*) cinput,nfskp
    read(9,*) cinput,nrestart
    read(9,*) cinput,nvarout
    read(9,*) cinput,reoo,tempoo
    read(9,*) cinput,amach1,amach2,amach3
    read(9,*) cinput,wtemp
    read(9,*) cinput,cfl
    read(9,*) cinput,tmax,timf,tsam
    read(9,*) cinput,fltk
    read(9,*) cinput,dto
    close(9)

    cinput=cinput; fltk=pi*fltk
    rhooo=1; poo=1/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1**2+amach2**2+amach3**2)
    srefoo=111.0_nr/tempoo; srefp1dre=(srefoo+1)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=1/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo
    ! rpt-Initialising the record count 
    nwrec=0

    allocate(times(0:ndata))
    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),npc(0:mbk,3))

    call inputext

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

    mo(0)=0
 do mm=1,mbk
    mo(mm)=mo(mm-1)+npc(mm-1,1)*npc(mm-1,2)*npc(mm-1,3)
 end do
 do mm=0,mbk
 if(myid>=mo(mm)) then; mb=mm; end if
 end do
    lxio=lximb(mb); leto=letmb(mb); lzeo=lzemb(mb)

    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10); cno=achar(no+48)
    czone='zone'//cno(2)//cno(1)//cno(0)
    coutput='out/output'//cno(2)//cno(1)//cno(0)//'.plt'
    cgrid='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat'
    crestart='rsta/restart'//cno(2)//cno(1)//cno(0)//'.dat'
    no(4)=myid/10000; no(3)=mod(myid,10000)/1000; no(2)=mod(myid,1000)/100
    no(1)=mod(myid,100)/10; no(0)=mod(myid,10); cno=achar(no+48)
    cdata='misc/data'//cno(4)//cno(3)//cno(2)//cno(1)//cno(0)//'.dat'
    cturb='misc/turb'//cno(4)//cno(3)//cno(2)//cno(1)//cno(0)//'.dat'

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

!===== ALLOCATION OF MAIN ARRAYS

    allocate(qo(0:lmx,5),qa(0:lmx,5),de(0:lmx,5))
    allocate(xim(0:lmx,3),etm(0:lmx,3),zem(0:lmx,3),rr(0:lmx,3),ss(0:lmx,3))
    ! rpt-dA allocation added
    allocate(p(0:lmx),yaco(0:lmx),varr(0:lmx),dA(0:lmx))

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
    allocate(cmm1(0:ii,0:1),cmm2(0:jj,0:1),cmm3(0:kk,0:1))

    allocate(xu(0:lim,3),yu(0:lim,3),xl(0:lim,2),yl(0:lim,2),li(0:lim),sa(0:lim),sb(0:lim))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

    albed(-2:2,0,0)=(/zero,zero,one,alpha01,beta02/)
    albed(-2:2,1,0)=(/zero,alpha10,one,alpha12,beta13/)
    albed(-2:2,2,0)=(/beta20,alpha21,one,alpha23,beta24/)

    albed(-2:2,0,1)=(/zero,zero,one,alpha,beta/)
    albed(-2:2,1,1)=(/zero,alpha,one,alpha,beta/)
    albed(-2:2,2,1)=(/beta,alpha,one,alpha,beta/)

    call fcbcm(fltk,albef(:,:,0),fam(:),fbm(:),fcm(:))
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
 do i=1,mbci
    sbci(i)=sbci(i)*(1-real(i,nr)/(mbci+1))
 end do

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

!===== GRID INPUT & CALCULATION OF GRID METRICS

    allocate(lio(0:let,0:lze))
 do k=0,lze; kp=k*(leto+1)*(lxio+1)
 do j=0,let; jp=j*(lxio+1)
    lio(j,k)=jp+kp
 end do
 end do
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

    call spanLoc(2)

    rr(:,1)=ss(:,1)
    m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)

    rr(:,1)=ss(:,2)
    m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)

    rr(:,1)=ss(:,3)
    m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    de(:,1)=rr(:,1); de(:,2)=rr(:,2); de(:,3)=rr(:,3)

    call clComp(0,1,2)
    call clComp(0,2,2)
    

    xim(:,1)=qa(:,2)*de(:,3)-de(:,2)*qa(:,3)
    xim(:,2)=de(:,2)*qo(:,3)-qo(:,2)*de(:,3)
    xim(:,3)=qo(:,2)*qa(:,3)-qa(:,2)*qo(:,3)
    etm(:,1)=qa(:,3)*de(:,1)-de(:,3)*qa(:,1)
    etm(:,2)=de(:,3)*qo(:,1)-qo(:,3)*de(:,1)
    etm(:,3)=qo(:,3)*qa(:,1)-qa(:,3)*qo(:,1)
    zem(:,1)=qa(:,1)*de(:,2)-de(:,1)*qa(:,2)
    zem(:,2)=de(:,1)*qo(:,2)-qo(:,1)*de(:,2)
    zem(:,3)=qo(:,1)*qa(:,2)-qa(:,1)*qo(:,2)


!    rr(:,3)=qa(:,2)*ss(:,3); rr(:,2)=qa(:,3)*ss(:,3)
!    m=1; call mpigo(ntdrv,nrall,n45go,m); call deriv(3,3); call deriv(2,2); xim(:,m)=rr(:,3)-rr(:,2)
!    rr(:,3)=de(:,2)*ss(:,1); rr(:,2)=de(:,3)*ss(:,1)
!    m=2; call mpigo(ntdrv,nrall,n45go,m); call deriv(3,3); call deriv(2,2); xim(:,m)=rr(:,3)-rr(:,2)
!    rr(:,3)=qo(:,2)*ss(:,2); rr(:,2)=qo(:,3)*ss(:,2)
!    m=3; call mpigo(ntdrv,nrall,n45go,m); call deriv(3,3); call deriv(2,2); xim(:,m)=rr(:,3)-rr(:,2)
!
!    rr(:,1)=qa(:,3)*ss(:,3); rr(:,3)=qa(:,1)*ss(:,3)
!    m=1; call mpigo(ntdrv,nrall,n45go,m); call deriv(1,1); call deriv(3,3); etm(:,m)=rr(:,1)-rr(:,3)
!    rr(:,1)=de(:,3)*ss(:,1); rr(:,3)=de(:,1)*ss(:,1)
!    m=2; call mpigo(ntdrv,nrall,n45go,m); call deriv(1,1); call deriv(3,3); etm(:,m)=rr(:,1)-rr(:,3)
!    rr(:,1)=qo(:,3)*ss(:,2); rr(:,3)=qo(:,1)*ss(:,2)
!    m=3; call mpigo(ntdrv,nrall,n45go,m); call deriv(1,1); call deriv(3,3); etm(:,m)=rr(:,1)-rr(:,3)
!
!    rr(:,2)=qa(:,1)*ss(:,3); rr(:,1)=qa(:,2)*ss(:,3)
!    m=1; call mpigo(ntdrv,nrall,n45go,m); call deriv(2,2); call deriv(1,1); zem(:,m)=rr(:,2)-rr(:,1)
!    rr(:,2)=de(:,1)*ss(:,1); rr(:,1)=de(:,2)*ss(:,1)
!    m=2; call mpigo(ntdrv,nrall,n45go,m); call deriv(2,2); call deriv(1,1); zem(:,m)=rr(:,2)-rr(:,1)
!    rr(:,2)=qo(:,1)*ss(:,2); rr(:,1)=qo(:,2)*ss(:,2)
!    m=3; call mpigo(ntdrv,nrall,n45go,m); call deriv(2,2); call deriv(1,1); zem(:,m)=rr(:,2)-rr(:,1)

    yaco(:)=3/(qo(:,1)*xim(:,1)+qo(:,2)*etm(:,1)+qo(:,3)*zem(:,1)&
              +qa(:,1)*xim(:,2)+qa(:,2)*etm(:,2)+qa(:,3)*zem(:,2)&
              +de(:,1)*xim(:,3)+de(:,2)*etm(:,3)+de(:,3)*zem(:,3))


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

!===== SETTING UP OUTPUT FILE & STORING GRID DATA

    inquire(iolength=lp) varr
    open(0,file=cdata,access='direct',recl=lp)
 do nn=1,3
 ! rpt-Increasigng the record count
 nwrec=nwrec+1
    varr(:)=ss(:,nn); write(0,rec=nwrec) varr(:)
 end do
 ! rpt-If ngrid is 1 then store the grid metrics and increase the record count
 if (ngridv==1) then
 do nn=1,3
 nwrec=nwrec+1
    varr(:)=xim(:,nn); write(0,rec=nwrec) varr(:)
 end do
 do nn=1,3
 nwrec=nwrec+1
    varr(:)=etm(:,nn); write(0,rec=nwrec) varr(:)
 end do
 do nn=1,3
 nwrec=nwrec+1
    varr(:)=zem(:,nn); write(0,rec=nwrec) varr(:)
 end do
 end if

 nrec=nwrec

!===== SETTING UP SPONGE ZONE PARAMETERS

    call spongeup

!===== INITIAL CONDITIONS

 if(nts==0) then
    n=0; ndt=0; dt=0; dts=0; dte=0; timo=0
    call initialo ! Make sure that pressure is calculated.
 else
    CALL MPI_BARRIER(icom,ierr)
    open(3,file=crestart,access='stream',shared)
    lh=0
    read(3,pos=4*lh+1) n;    lh=lh+2
    read(3,pos=4*lh+1) ndt;  lh=lh+2
    read(3,pos=4*lh+1) dt;   lh=lh+2
    read(3,pos=4*lh+1) dts;  lh=lh+2
    read(3,pos=4*lh+1) dte;  lh=lh+2
    read(3,pos=4*lh+1) timo; lh=lh+2

    tmp=(tsam-timo)/timo
    if (abs(tmp)<0.1e0) then
    tsam=timo
    end if

    lp=lpos(myid)+lh
 do m=1,5; lq=(m-1)*ltomb
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    read(3,pos=nr*(lp+lq+lio(j,k))+1) qa(l:l+lxi,m)
 end do; end do
 end do
    close(3)
    p(:)=gamm1*(qa(:,5)-half*(qa(:,2)**2+qa(:,3)**2+qa(:,4)**2)/qa(:,1))
 end if

!============================================
!===== BEGINNING OF TIME MARCHING IN SOLUTION
!============================================

    wts=MPI_WTIME()

 if(myid==0) then
    open(1,file='signal.dat'); close(1,status='delete')
 end if
    call MPI_BARRIER(icom,ierr)
    open(1,file='signal.dat',access='direct',form='formatted',recl=16,shared)

    ndati=-1; nsigi=-1
 do while(timo-tmax<0.and.(dt/=0.or.n<=2))

 call clComp(1,1,2)
 call clComp(1,2,2)
 call clComp(1,1,1)
 call clComp(1,2,1)

 if(myid==0.and.mod(n,nscrn)==0) then
    write(*,"(' n =',i8,'   time =',f12.5,' Cl1 = ',f10.5,', Cl2 = ',f10.5)") &
    n,timo,cl(1,2),cl(2,2)
 end if


    qo(:,:)=qa(:,:)

!-----------------------------------
!----- NKRK-STAGE RUNGE-KUTTA STAGES
!-----------------------------------

 do nk=1,nkrk

!----- MOVING FRAME VELOCITY & ACCELERATION BEFORE TIME ADVANCING

    dtko=min(max(nk-2,0),1)*dt/(nkrk-nk+3); dtk=min(nk-1,1)*dt/(nkrk-nk+2)
    call movef(dtko,dtk)

!----- TEMPORARY STORAGE OF 1/DENSITY, VELOCITY & TEMPERATURE

    de(:,1)=1/qa(:,1)
    de(:,2)=qa(:,2)*de(:,1)
    de(:,3)=qa(:,3)*de(:,1)
    de(:,4)=qa(:,4)*de(:,1)
    de(:,5)=gam*p(:)*de(:,1)
    ss(:,1)=srefp1dre*de(:,5)**1.5_nr/(de(:,5)+srefoo)

!----- DETERMINATION OF TIME STEP SIZE & OUTPUT TIME

 if(nk==1) then
 if(mod(n,10)==1) then; ndt=n; dts=dte
 if(dto<0) then
    rr(:,1)=xim(:,1)*xim(:,1)+xim(:,2)*xim(:,2)+xim(:,3)*xim(:,3)&
           +etm(:,1)*etm(:,1)+etm(:,2)*etm(:,2)+etm(:,3)*etm(:,3)&
           +zem(:,1)*zem(:,1)+zem(:,2)*zem(:,2)+zem(:,3)*zem(:,3)
    rr(:,2)=abs(xim(:,1)*(de(:,2)+umf(1))+xim(:,2)*(de(:,3)+umf(2))+xim(:,3)*(de(:,4)+umf(3)))&
           +abs(etm(:,1)*(de(:,2)+umf(1))+etm(:,2)*(de(:,3)+umf(2))+etm(:,3)*(de(:,4)+umf(3)))&
           +abs(zem(:,1)*(de(:,2)+umf(1))+zem(:,2)*(de(:,3)+umf(2))+zem(:,3)*(de(:,4)+umf(3)))
    ss(:,2)=abs(yaco(:))
    res=maxval((sqrt(de(:,5)*rr(:,1))+rr(:,2))*ss(:,2))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    ra0=cfl/fctr; ra1=ra0
 if(nviscous==1) then
    res=maxval(de(:,1)*ss(:,1)*rr(:,1)*ss(:,2)*ss(:,2))
    call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr)
    ra1=half/fctr
 end if
    dte=min(ra0,ra1)
 else
    dte=dto
 end if
 end if
    dt=dts+(dte-dts)*sin(0.05_nr*pi*(n-ndt))**2

    nout=0; res=tsam+(ndati+1)*(tmax-tsam)/ndata
 if((timo-res)*(timo+dt-res)<=0) then
    nout=1; ndati=ndati+1
 end if
 ! rpt-If cl goes to infty the calculation is crashed and 
 !     needs to save result to check what happened
 if ((cl(1,2).ne.cl(1,2)).or.(cl(2,2).ne.cl(2,2))) then
   nout = 2; ndati = ndati+1;
 end if

 end if

!----- VISCOUS SHEAR STRESSES & HEAT FLUXES

 if(nviscous==1) then
    de(:,1)=ss(:,1)

    rr(:,1)=de(:,2)
    m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    txx(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    hzz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    tzx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,3)
    m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    txy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    tyy(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    hxx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,4)
    m=4; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    hyy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    tyz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    tzz(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,5)
    m=5; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    ss(:,1)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    ss(:,2)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    ss(:,3)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    fctr=2.0_nr/3
    rr(:,1)=de(:,1)*yaco(:)
    rr(:,2)=gamm1prndtli*rr(:,1)
    de(:,5)=fctr*(txx(:)+tyy(:)+tzz(:))

    txx(:)=rr(:,1)*(2*txx(:)-de(:,5))
    tyy(:)=rr(:,1)*(2*tyy(:)-de(:,5))
    tzz(:)=rr(:,1)*(2*tzz(:)-de(:,5))
    txy(:)=rr(:,1)*(txy(:)+hzz(:))
    tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
    tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
    hxx(:)=rr(:,2)*ss(:,1)+de(:,2)*txx(:)+de(:,3)*txy(:)+de(:,4)*tzx(:)
    hyy(:)=rr(:,2)*ss(:,2)+de(:,2)*txy(:)+de(:,3)*tyy(:)+de(:,4)*tyz(:)
    hzz(:)=rr(:,2)*ss(:,3)+de(:,2)*tzx(:)+de(:,3)*tyz(:)+de(:,4)*tzz(:)
 end if

!----- CALCULATION OF FLUX DERIVATIVES

    rr(:,1)=de(:,2)+umf(1)
    rr(:,2)=de(:,3)+umf(2)
    rr(:,3)=de(:,4)+umf(3)
    ss(:,1)=xim(:,1)*rr(:,1)+xim(:,2)*rr(:,2)+xim(:,3)*rr(:,3)
    ss(:,2)=etm(:,1)*rr(:,1)+etm(:,2)*rr(:,2)+etm(:,3)*rr(:,3)
    ss(:,3)=zem(:,1)*rr(:,1)+zem(:,2)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=qa(:,1)*ss(:,1)
    rr(:,2)=qa(:,1)*ss(:,2)
    rr(:,3)=qa(:,1)*ss(:,3)
    m=1; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1); call deriv(2,2); call deriv(3,3)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,2)*ss(:,1)+xim(:,1)*p(:)
    rr(:,2)=qa(:,2)*ss(:,2)+etm(:,1)*p(:)
    rr(:,3)=qa(:,2)*ss(:,3)+zem(:,1)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txx(:)-xim(:,2)*txy(:)-xim(:,3)*tzx(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txx(:)-etm(:,2)*txy(:)-etm(:,3)*tzx(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txx(:)-zem(:,2)*txy(:)-zem(:,3)*tzx(:)
 end if
    m=2; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1); call deriv(2,2); call deriv(3,3)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,3)*ss(:,1)+xim(:,2)*p(:)
    rr(:,2)=qa(:,3)*ss(:,2)+etm(:,2)*p(:)
    rr(:,3)=qa(:,3)*ss(:,3)+zem(:,2)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txy(:)-xim(:,2)*tyy(:)-xim(:,3)*tyz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txy(:)-etm(:,2)*tyy(:)-etm(:,3)*tyz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txy(:)-zem(:,2)*tyy(:)-zem(:,3)*tyz(:)
 end if
    m=3; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1); call deriv(2,2); call deriv(3,3)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,4)*ss(:,1)+xim(:,3)*p(:)
    rr(:,2)=qa(:,4)*ss(:,2)+etm(:,3)*p(:)
    rr(:,3)=qa(:,4)*ss(:,3)+zem(:,3)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*tzx(:)-xim(:,2)*tyz(:)-xim(:,3)*tzz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*tzx(:)-etm(:,2)*tyz(:)-etm(:,3)*tzz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*tzx(:)-zem(:,2)*tyz(:)-zem(:,3)*tzz(:)
 end if
    m=4; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1); call deriv(2,2); call deriv(3,3)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    de(:,5)=qa(:,5)+p(:)
    rr(:,1)=de(:,5)*ss(:,1)-p(:)*(umf(1)*xim(:,1)+umf(2)*xim(:,2)+umf(3)*xim(:,3))
    rr(:,2)=de(:,5)*ss(:,2)-p(:)*(umf(1)*etm(:,1)+umf(2)*etm(:,2)+umf(3)*etm(:,3))
    rr(:,3)=de(:,5)*ss(:,3)-p(:)*(umf(1)*zem(:,1)+umf(2)*zem(:,2)+umf(3)*zem(:,3))
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*hxx(:)-xim(:,2)*hyy(:)-xim(:,3)*hzz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*hxx(:)-etm(:,2)*hyy(:)-etm(:,3)*hzz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*hxx(:)-zem(:,2)*hyy(:)-zem(:,3)*hzz(:)
 end if
    m=5; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1); call deriv(2,2); call deriv(3,3)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

!----- IMPLEMENTATION OF SPONGE CONDITION

    call spongego ! Make sure that "ss(l,1)=0" is specified if sponge is NOT used.

!----- PREPARATION FOR GCBC & GCIC

 do nn=1,3; nz=min(nn-1,1)
 select case(nn)
 case(1); drva=>drva1; cm=>cm1; case(2); drva=>drva2; cm=>cm2; case(3); drva=>drva3; cm=>cm3
 end select
 do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn); iq=1-2*ip
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); call xtq2r(cm(jk,:,ip)); drva(jk,:,ip)=matmul(xt(:,:),yaco(l)*de(l,:))
    rr(l,1)=nz*rr(l,1)+1
 do ii=1,mbci; ll=indx3(i+iq*ii,j,k,nn)
    rr(ll,1)=nz*rr(ll,1)+1
 end do
 end do
 end do
 end if
 end do
 end do

!----- INTERNODE COMMNICATION FOR GCIC

    ir=0; itag=30
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; case(2); drva=>drva2; drvb=>drvb2; case(3); drva=>drva3; drvb=>drvb3
 end select
 do ip=0,1; iq=1-ip; np=nbc(ip,nn)
 if(np==30) then
    ir=ir+1; call MPI_ISEND(drva(:,:,ip),5*nbsize(nn),MPI_REAL8,ncd(ip,nn),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(drvb(:,:,ip),5*nbsize(nn),MPI_REAL8,ncd(ip,nn),itag+ip,icom,ireq(ir),ierr)
 end if
 end do
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end iF

!----- IMPLEMENTATION OF GCBC & GCIC

 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; cm=>cm1
 case(2); drva=>drva2; drvb=>drvb2; cm=>cm2
 case(3); drva=>drva3; drvb=>drvb3; cm=>cm3
 end select
 do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn); iq=1-2*ip
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); cha(:)=drva(jk,:,ip); dha(:)=drvb(jk,:,ip)
 select case(np)
 case(10)
 if(iq*(vn+vs+ao)>0) then; cha(4)=-cha(5)+ss(l,1)*(p(l)-poo); end if
 if(iq*(vn+vs-ao)>0) then; cha(5)=-cha(4)+ss(l,1)*(p(l)-poo); end if
 case(20,25)
    cha(4+ip)=cha(5-ip)+iq*aoi*qa(l,1)*(2*sum(cm(jk,:,ip)*dudtmf(:))+100*(vn+vs))
 case(30)
    cha(:)=half*(cha(:)+dha(:))
 end select
    call xtr2q(cm(jk,:,ip)); res=1/yaco(l); dha(:)=res*matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
    res=1/rr(l,1); de(l,:)=de(l,:)+res*dha(:)
 do ii=1,mbci; ll=indx3(i+iq*ii,j,k,nn)
    res=sbci(ii)/rr(ll,1); de(ll,:)=de(ll,:)+res*dha(:)
 end do
 end do
 end do
 end if
 end do
 end do

!----- UPDATING CONSERVATIVE VARIABLES

    dtko=min(nk-1,1)*dt/(nkrk-nk+2); dtk=dt/(nkrk-nk+1)
    call movef(dtko,dtk)

    rr(:,1)=dtk*yaco(:)
    qa(:,1)=qo(:,1)-rr(:,1)*de(:,1)
    qa(:,2)=qo(:,2)-rr(:,1)*de(:,2)
    qa(:,3)=qo(:,3)-rr(:,1)*de(:,3)
    qa(:,4)=qo(:,4)-rr(:,1)*de(:,4)
    qa(:,5)=qo(:,5)-rr(:,1)*de(:,5)

!----- WALL TEMPERATURE & VELOCITY CONDITION

    ra0=ham*hamm1*wtemp
 do nn=1,3; do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn)
 if((np-20)*(np-25)==0) then; ns=(25-np)/5; ne=1-ns
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    qa(l,2:4)=ns*qa(l,2:4)-ne*umf(:)*qa(l,1)
    qa(l,5)=ra0*qa(l,1)+half*(qa(l,2)*qa(l,2)+qa(l,3)*qa(l,3)+qa(l,4)*qa(l,4))/qa(l,1)
 end do
 end do
 end if
 end do; end do

!----- JUNCTION & INTERFACE AVERAGING

    call junction

    ir=0; itag=30
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; case(2); drva=>drva2; drvb=>drvb2; case(3); drva=>drva3; drvb=>drvb3
 end select
 do ip=0,1; iq=1-ip; np=nbc(ip,nn); i=ip*ijk(1,nn)
 if((np-30)*(np-35)*(np-45)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    drva(jk,:,ip)=qa(l,:); rr(l,1)=1
 end do
 end do
    ir=ir+1; call MPI_ISEND(drva(:,:,ip),5*nbsize(nn),MPI_REAL8,ncd(ip,nn),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(drvb(:,:,ip),5*nbsize(nn),MPI_REAL8,ncd(ip,nn),itag+ip,icom,ireq(ir),ierr)
 end if
 end do
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if
 do nn=1,3
 select case(nn); case(1); drvb=>drvb1; case(2); drvb=>drvb2; case(3); drvb=>drvb3; end select
 do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn)
 if((np-30)*(np-35)*(np-45)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    rr(l,1)=rr(l,1)+1; rr(l,2)=1/rr(l,1); qa(l,:)=rr(l,2)*((rr(l,1)-1)*qa(l,:)+drvb(jk,:,ip))
 end do
 end do
 end if
 end do
 end do

!----- FILTERING

 if(mod(n,nfskp+1)==0.and.nk==nkrk) then
 do m=1,5
    rr(:,1)=qa(:,m)
    call mpigo(ntflt,nrone,n45no,m); call filte(1,1); call filte(2,1); call filte(3,1)
    qa(:,m)=rr(:,1)
 end do
 end if

!----- UPDATING PRESSURE

    p(:)=gamm1*(qa(:,5)-half*(qa(:,2)*qa(:,2)+qa(:,3)*qa(:,3)+qa(:,4)*qa(:,4))/qa(:,1))

!-------------------------------
!----- END OF RUNGE-KUTTA STAGES
!-------------------------------

 end do

!---------------------
!----- ADVANCE IN TIME
!---------------------

    n=n+1
    timo=timo+dt

!----- RECORDING INTERMEDIATE RESULTS
 if(nout==1.or.nout==2) then
      times(ndati)=timo
   if(myid==0) then
      write(*,"('===> saving output ',i3,' at time =',f12.8)") ndati,timo
   end if
   if (nout==2) then
      qa(:,:)=qo(:,:)   
   end if

   !==========SAVING VELOCITY AND DENSITY
   nwrec=nwrec+1
      varr(:)=qa(:,1); write(0,rec=nwrec) varr(:)
   do m = 2, 4
   nwrec=nwrec+1
      varr(:)=((qa(:,m)/qa(:,1))+umf(m-1)); write(0,rec=nwrec) varr(:)
   end do
   !======================================

   select case(nvarout)
   case(0); rr(:,1)=qa(:,1)
   case(1); rr(:,1)=qa(:,2)/qa(:,1)+umf(1)
   case(2); rr(:,1)=qa(:,3)/qa(:,1)+umf(2)
   case(3); rr(:,1)=qa(:,4)/qa(:,1)+umf(3)
   case(4); rr(:,1)=p(:)
   case(5); rr(:,1)=gam*p(:)-1
   case(6)
      ss(:,1)=1/qa(:,1); de(:,1:3)=0

      rr(:,1)=ss(:,1)*qa(:,2)
      m=1; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
      de(:,2)=de(:,2)+rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3)
      de(:,3)=de(:,3)-rr(:,1)*xim(:,2)-rr(:,2)*etm(:,2)-rr(:,3)*zem(:,2)

      rr(:,1)=ss(:,1)*qa(:,3)
      m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
      de(:,3)=de(:,3)+rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1)
      de(:,1)=de(:,1)-rr(:,1)*xim(:,3)-rr(:,2)*etm(:,3)-rr(:,3)*zem(:,3)

      rr(:,1)=ss(:,1)*qa(:,4)
      m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
      de(:,1)=de(:,1)+rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2)
      de(:,2)=de(:,2)-rr(:,1)*xim(:,1)-rr(:,2)*etm(:,1)-rr(:,3)*zem(:,1)

      rr(:,1)=sqrt((de(:,1)*de(:,1)+de(:,2)*de(:,2)+de(:,3)*de(:,3))*yaco(:)*yaco(:))
   case(7)
      nwrec=nwrec+1
      varr(:)=p(:); write(0,rec=nwrec) varr(:)

      ss(:,1)=1/qa(:,1); de(:,1:3)=0

      rr(:,1)=ss(:,1)*qa(:,2)
      m=1; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
      de(:,1)=yaco(:)*(rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1))
      de(:,2)=yaco(:)*(rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2))
      de(:,3)=yaco(:)*(rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3))

      nwrec=nwrec+1
      varr(:)=de(:,1); write(0,rec=nwrec) varr(:)
      nwrec=nwrec+1
      varr(:)=de(:,2); write(0,rec=nwrec) varr(:)
      nwrec=nwrec+1
      varr(:)=de(:,3); write(0,rec=nwrec) varr(:)

      rr(:,1)=ss(:,1)*qa(:,3)
      m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
      de(:,1)=yaco(:)*(rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1))
      de(:,2)=yaco(:)*(rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2))
      de(:,3)=yaco(:)*(rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3))

      nwrec=nwrec+1
      varr(:)=de(:,1); write(0,rec=nwrec) varr(:)
      nwrec=nwrec+1
      varr(:)=de(:,2); write(0,rec=nwrec) varr(:)
      nwrec=nwrec+1
      varr(:)=de(:,3); write(0,rec=nwrec) varr(:)

      rr(:,1)=ss(:,1)*qa(:,4)
      m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
      de(:,1)=yaco(:)*(rr(:,1)*xim(:,1)+rr(:,2)*etm(:,1)+rr(:,3)*zem(:,1))
      de(:,2)=yaco(:)*(rr(:,1)*xim(:,2)+rr(:,2)*etm(:,2)+rr(:,3)*zem(:,2))
      de(:,3)=yaco(:)*(rr(:,1)*xim(:,3)+rr(:,2)*etm(:,3)+rr(:,3)*zem(:,3))

      nwrec=nwrec+1
      varr(:)=de(:,1); write(0,rec=nwrec) varr(:)
      nwrec=nwrec+1
      varr(:)=de(:,2); write(0,rec=nwrec) varr(:)
      rr(:,1)=de(:,3)
   end select

   nwrec=nwrec+1
   varr(:)=rr(:,1); write(0,rec=nwrec) varr(:)
   narec=nwrec-nrec

   if (nout==2) then
      goto 100
   end if

   !===== GENERATING RESTART DATA FILE
   
   if(nrestart==1) then
      if (myid==0) then
         write(*,*) 'Writting restart file..'
      end if
      if(myid==mo(mb)) then
         open(3,file=crestart); close(3,status='delete')
      end if
         call MPI_BARRIER(icom,ierr)
         open(3,file=crestart,access='stream',shared)
      if(myid==mo(mb)) then
      lh=0
         write(3,pos=4*lh+1) n;    lh=lh+2
         write(3,pos=4*lh+1) ndt;  lh=lh+2
         write(3,pos=4*lh+1) dt;   lh=lh+2
         write(3,pos=4*lh+1) dts;  lh=lh+2
         write(3,pos=4*lh+1) dte;  lh=lh+2
         write(3,pos=4*lh+1) timo; lh=lh+2
      end if
         lp=lpos(myid)+12
      do m=1,5; lq=(m-1)*ltomb
      do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         write(3,pos=nr*(lp+lq+lio(j,k))+1) qa(l:l+lxi,m)
      end do; end do
      end do
         close(3)
   end if
 end if

 if(timo-tsam>=0.and.mod(n,nsgnl)==0) then
    nsigi=nsigi+1; call signalgo
 end if

!==========================
!===== END OF TIME MARCHING
!==========================
 end do
! rpt-tag for exiting while loop
100 continue

    close(1)
    ! rpt-Close cl and cd files
    call clComp(2,1,2)
    call clComp(2,2,2)

    wte=MPI_WTIME(); res=wte-wts
    call MPI_ALLREDUCE(res,wtime,1,MPI_REAL8,MPI_SUM,icom,ierr)
 if(myid==0) then
    open(9,file='walltime.dat',position='append')
    write(9,'(2es15.7)') real(npro,nr),wtime/npro
    close(9)
 end if

 call post
 call cpComp(2)

!===== END OF JOB

 if(myid==0) then
    write(*,*) "Finished."
 end if

 if(myid==mo(0)+npc(0,1)*npc(0,2)-1) then
    rr(:,1)=0; rr(:,2)=0
 do k=0,lze; do j=0,let; do i=1,lxi-1; l=indx3(i,j,k,1)
    rv(1)=varr(indx3(i-1,j,k,1)); rv(2)=varr(l); rv(3)=varr(indx3(i+1,j,k,1))
    rr(l,1)=abs(rv(1)-2*rv(2)+rv(3))/sqrt(rv(1)**2+rv(2)**2+rv(3)**2)
 end do; end do; end do
 do k=0,lze; do j=1,let-1; do i=0,lxi; l=indx3(i,j,k,1)
    rv(1)=varr(indx3(i,j-1,k,1)); rv(2)=varr(l); rv(3)=varr(indx3(i,j+1,k,1))
    rr(l,2)=abs(rv(1)-2*rv(2)+rv(3))/sqrt(rv(1)**2+rv(2)**2+rv(3)**2)
 end do; end do; end do
    write(*,*) max(maxval(rr(:,1)),maxval(rr(:,2)))
 end if

    call MPI_FINALIZE(ierr)

 end program main3d

!*****
