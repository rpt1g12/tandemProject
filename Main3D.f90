!*****
!***** 3D PARALLEL SOLVER
!*****

 program main3d

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
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
    read(9,*) cinput,fltk,fltkbco,fltkbcm
    read(9,*) cinput,dto
    close(9)

    cinput=cinput; fltk=pi*fltk; fltkbco=pi*fltkbco; fltkbcm=pi*fltkbcm
    rhooo=1; poo=1/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1**2+amach2**2+amach3**2)
    srefoo=111.0_nr/tempoo; srefp1dre=(srefoo+1)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=1/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo
    ! rpt-Initialising the record count
    nrec=0

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
    coutput='output'//cno(2)//cno(1)//cno(0)//'.plt'
    cgrid='misc/grid'//cno(2)//cno(1)//cno(0)//'.dat'
    crestart='misc/restart'//cno(2)//cno(1)//cno(0)//'.dat'
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
    allocate(cmm1(0:ii,0:1),cmm2(0:jj,0:1),cmm3(0:kk,0:1))

    allocate(xu(0:lim,3),yu(0:lim,3),xl(0:lim,2),yl(0:lim,2),li(0:lim),sa(0:lim),sb(0:lim))

!===== EXTRA COEFFICIENTS FOR DOMAIN BOUNDARIES

    albed(-2:2,0,-1)=(/zero,zero,one,alpha01,beta02/)
    albed(-2:2,1,-1)=(/zero,alpha10,one,alpha12,beta13/)
    albed(-2:2,2,-1)=(/beta20,alpha21,one,alpha23,beta24/)

    albed(:,:,0)=albed(:,:,-1)
    abc(:,0)=(/a01,a02,a03,a04,a05,a06/)
    abc(:,1)=(/a10,a12,a13,a14,a15,a16/)
    abc(:,2)=(/a20,a21,a23,a24,a25,a26/)

    albed(-2:2,0,1)=(/zero,zero,one,alpha,beta/)
    albed(-2:2,1,1)=(/zero,alpha,one,alpha,beta/)
    albed(-2:2,2,1)=(/beta,alpha,one,alpha,beta/)

    call fcbcm(fltk,fltkbcm,albef(:,:,-1),fam(:),fbm(:),fcm(:))
    call fcbco(fltk,fltkbco,albef(:,:,0),fbc(:))
    call fcint(fltk,half,alphf,betf,fa,fb,fc)
    albef(-2:2,0,1)=(/zero,zero,one,alphf,betf/)
    albef(-2:2,1,1)=(/zero,alphf,one,alphf,betf/)
    albef(-2:2,2,1)=(/betf,alphf,one,alphf,betf/)

    pbco(:,:,:)=0; pbci(:,:,:)=0; call sbcco(0,-1)
 do nt=0,1; do j=0,1; ii=lmd+nt*(lmf-lmd)
    pbcot(j,nt)=sum(pbco(0:ii,j,nt))
 end do; end do

!===== PENTADIAGONAL MATRICES FOR DIFFERENCING & FILETERING

 do nn=1,3
 select case(nn)
 case(1); is=0; ie=is+lxi; case(2); is=lxi+1; ie=is+let; case(3); is=lxi+let+2; ie=is+lze
 end select
 do ip=0,1; np=nbc(ip,nn)
 select case(np)
 case(10,20,25,30); ndf(ip,0,nn)=0; ndf(ip,1,nn)=-1
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

    rr(:,1)=ss(:,1)
    m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)

    rr(:,1)=ss(:,2)
    m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)

    rr(:,1)=ss(:,3)
    m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
    de(:,1)=rr(:,1); de(:,2)=rr(:,2); de(:,3)=rr(:,3)

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
 case(1); rv(:)=yaco(l)*xim(l,:)
    cmm1(jk,ip)=sqrt(rv(1)**2+rv(2)**2+rv(3)**2); fctr=1/cmm1(jk,ip); cm1(jk,:,ip)=fctr*rv(:)
 case(2); rv(:)=yaco(l)*etm(l,:)
    cmm2(jk,ip)=sqrt(rv(1)**2+rv(2)**2+rv(3)**2); fctr=1/cmm2(jk,ip); cm2(jk,:,ip)=fctr*rv(:)
 case(3); rv(:)=yaco(l)*zem(l,:)
    cmm3(jk,ip)=sqrt(rv(1)**2+rv(2)**2+rv(3)**2); fctr=1/cmm3(jk,ip); cm3(jk,:,ip)=fctr*rv(:)
 end select
 end do
 end do
 end do; end do

!===== SETTING UP OUTPUT FILE & STORING GRID DATA

    inquire(iolength=lp) varr
    open(0,file=cdata,access='direct',recl=lp)
 do nn=1,3
 ! rpt-Increasigng the record count
 nrec=nrec+1
    varr(:)=ss(:,nn); write(0,rec=nrec) varr(:)
 end do
 ! rpt-If ngrid is 1 then store the grid metrics and increase the record count
 if (ngridv==1) then
 do nn=1,3
 nrec=nrec+1
    varr(:)=xim(:,nn); write(0,rec=nrec) varr(:)
 end do
 do nn=1,3
 nrec=nrec+1
    varr(:)=etm(:,nn); write(0,rec=nrec) varr(:)
 end do
 do nn=1,3
 nrec=nrec+1
    varr(:)=zem(:,nn); write(0,rec=nrec) varr(:)
 end do
 end if

!===== SETTING UP SPONGE ZONE PARAMETERS

    call spongeup

!===== INITIAL CONDITIONS

 if(nts==0) then
    n=0; ndt=0; dt=0; dts=0; dte=0; timo=0
    call initialo ! Make sure that pressure is calculated.
 else
    open(9,file=crestart,access='stream',shared)
    read(9,pos=1) n; read(9,pos=2) ndt; read(9,pos=3) dt
    read(9,pos=4) dts; read(9,pos=5) dte; read(9,pos=6) timo
    lp=lpos(myid)+6
 do m=1,5; lq=(m-1)*ltomb
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    read(9,pos=nr*(lp+lq+lio(j,k))+1) qa(l:l+lxi,m)
 end do; end do
 end do
    close(9)
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

 if(myid==0.and.mod(n,nscrn)==0) then
    write(*,"(' n =',i8,'   time =',f12.5)") n,timo
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
 ! rpt-If it crashes exit the loop and save last results recorded
 if (n>2.and.dt==0) then
    goto 100
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

    call spongego

!----- IMPLEMENTATION OF GCBC & GCIC

    rr(:,1)=0
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; cm=>cm1; cmm=>cmm1
 case(2); drva=>drva2; cm=>cm2; cmm=>cmm2
 case(3); drva=>drva3; cm=>cm3; cmm=>cmm3
 end select
 do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn)
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); call xtq2r(cm(jk,:,ip))
    drva(jk,:,ip)=matmul(xt(:,:),yaco(l)*de(l,:)); rr(l,1)=rr(l,1)+cmm(jk,ip)
 end do
 end do
 end if
 end do
 end do
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
 end if
 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; drvb=>drvb1; cm=>cm1; cmm=>cmm1
 case(2); drva=>drva2; drvb=>drvb2; cm=>cm2; cmm=>cmm2
 case(3); drva=>drva3; drvb=>drvb3; cm=>cm3; cmm=>cmm3
 end select
 do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn); lp=1-2*ip
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); cha(:)=drva(jk,:,ip); dha(:)=drvb(jk,:,ip)
 select case(np)
 case(10)
 if(lp*(vn+vs+ao)>0) then; cha(4)=-cha(5); end if
 if(lp*(vn+vs-ao)>0) then; cha(5)=-cha(4); end if
 case(20,25)
    cha(4+ip)=cha(5-ip)+lp*aoi*qa(l,1)*(2*sum(cm(jk,:,ip)*dudtmf(:))+100*(vn+vs))
 case(30)
    fctr=min(max(half*(1+lp*aoi*(vn+vs)),zero),one)
    cha(1:3)=(1-fctr)*cha(1:3)+fctr*dha(1:3)
 if(lp*(vn+vs+ao)>0) then; cha(4)=dha(4); end if
 if(lp*(vn+vs-ao)>0) then; cha(5)=dha(5); end if
 end select
    call xtr2q(cm(jk,:,ip)); res=cmm(jk,ip)/(yaco(l)*rr(l,1))
    de(l,:)=de(l,:)+res*matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
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

!----- FILTERING

 if(mod(n,nfskp+1)==0.and.nk==nkrk) then
 do m=1,5
    rr(:,1)=qa(:,m)
    call mpigo(ntflt,nrone,n45no,m); call filte(1); call filte(2); call filte(3)
    qa(:,m)=rr(:,1)
 end do
 end if

!----- WALL TEMPERATURE & VELOCITY CONDITION

    ra0=ham*hamm1*wtemp
 do nn=1,3; do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn)
 if((np-20)*(np-25)==0) then; ns=(25-np)/5; ne=1-ns
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    qa(l,2:4)=ns*qa(l,2:4)+ne*umf(:)*qa(l,1)
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

 if(nout==1) then
    times(ndati)=timo
    if(myid==0) then
       write(*,"('**saved results** => n =',i8,'   time =',f12.5)") n,timo
    end if

    !==========SAVING VELOCITY AND DENSITY
    i=nrec
    do m = 2, 4
    i=i+1
       varr(:)=((qa(:,m)/qa(:,1))+umf(m-1)); write(0,rec=i) varr(:)
    end do
    i=i+1
       varr(:)=qa(:,1); write(0,rec=i) varr(:)
    !======================================

 select case(nvarout)
 case(0); rr(:,1)=qa(:,1)
 case(1); rr(:,1)=qa(:,2)/qa(:,1)+umf(1)
 case(2); rr(:,1)=qa(:,3)/qa(:,1)+umf(2)
 case(3); rr(:,1)=qa(:,4)/qa(:,1)+umf(3)
 case(4); rr(:,1)=p(:)
 case(5); rr(:,1)=gam*p(:)-1
 case(6)
    ss(:,1)=one/qa(:,1); de(:,1:3)=zero

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
 end select
    ! rpt-Because of saving metrics, record count changes
    !varr(:)=rr(:,1); write(0,rec=ndati+4) varr(:)
    !call intermed(3)
    i=i+1
    varr(:)=rr(:,1); write(0,rec=ndati+i) varr(:)


    !===== GENERATING RESTART DATA FILE
    
     if(nrestart==1) then
        if(myid==mo(mb)) then
           open(9,file=crestart); close(9,status='delete')
        end if
           call MPI_BARRIER(icom,ierr)
           open(9,file=crestart,access='stream',shared)
        if(myid==mo(mb)) then
           write(9,pos=1) timo
        end if
           lp=lpos(myid)+1
        do m=1,5; lq=(m-1)*ltomb
        do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
           write(9,pos=nr*(lp+lq+lio(j,k))+1) qa(l:l+lxi,m)
        end do; end do
        end do
           close(9)
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
nrec=nrec+4

    close(1)

    wte=MPI_WTIME(); res=wte-wts
    call MPI_ALLREDUCE(res,wtime,1,MPI_REAL8,MPI_SUM,icom,ierr)
 if(myid==0) then
    open(9,file='walltime.dat',position='append')
    write(9,'(2es15.7)') real(npro,nr),wtime/npro
    close(9)
 end if


!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

 if(dt==0) then
    ! rpt-File is not closed so it can still be saved up to the last record
    !close(0,status='delete')
    write(*,*) "Overflow."
 ndata=ndati
 end if
 !else
 if(myid==0) then
    write(*,*) "Writing Output files..."
 end if
    !call finalout
 if(myid==mo(mb)) then
    open(9,file=coutput); close(9,status='delete')
 end if
    call MPI_BARRIER(icom,ierr)
    open(9,file=coutput,access='stream',shared)
    lh=0
 if(myid==mo(mb)) then
    write(9,pos=4*lh+1) '#!TDV112'; lh=lh+2
    write(9,pos=4*lh+1) 1; lh=lh+1 ! Header Section
    write(9,pos=4*lh+1) 0; lh=lh+1 ! File Type
    cinput='title'; call strio(9,lh,cinput) ! File Title
    write(9,pos=4*lh+1) int4(ndata+ndatp+nrec+1); lh=lh+1 ! Number of Variables
    cinput='x'; call strio(9,lh,cinput)
    cinput='y'; call strio(9,lh,cinput)
    cinput='z'; call strio(9,lh,cinput)
    if (ngridv==1) then
       cinput='xix'; call strio(9,lh,cinput)
       cinput='xiy'; call strio(9,lh,cinput)
       cinput='xiz'; call strio(9,lh,cinput)
       cinput='etax'; call strio(9,lh,cinput)
       cinput='etay'; call strio(9,lh,cinput)
       cinput='etaz'; call strio(9,lh,cinput)
       cinput='zetax'; call strio(9,lh,cinput)
       cinput='zetay'; call strio(9,lh,cinput)
       cinput='zetaz'; call strio(9,lh,cinput)
    end if
    cinput='u'; call strio(9,lh,cinput)
    cinput='v'; call strio(9,lh,cinput)
    cinput='w'; call strio(9,lh,cinput)
    cinput='rho'; call strio(9,lh,cinput)
 do n=0,ndata+ndatp
    no(2)=n/100; no(1)=mod(n,100)/10; no(0)=mod(n,10); cno=achar(no+48)
    cinput='var'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
 end do
    write(9,pos=4*lh+1) 299.0; lh=lh+1 ! Zone Marker
    cinput=czone; call strio(9,lh,cinput)
    write(9,pos=4*lh+1) -1; lh=lh+1 ! Parent Zone
    write(9,pos=4*lh+1) -2; lh=lh+1 ! Strand ID
    write(9,pos=4*lh+1) dble(0.0); lh=lh+2 ! Solution Time (Double)
    write(9,pos=4*lh+1) -1; lh=lh+1 ! (Not used. Set to -1.)
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Zone Type
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Specify Var Location
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Raw Local 1-to-1 Face Neighbours Suppliled
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Number of Miscellaneous Face Neighbour Connections
    write(9,pos=4*lh+1) int4(lximb(mb)+1); lh=lh+1 ! IMax
    write(9,pos=4*lh+1) int4(letmb(mb)+1); lh=lh+1 ! JMax
    write(9,pos=4*lh+1) int4(lzemb(mb)+1); lh=lh+1 ! KMax
    write(9,pos=4*lh+1) 0; lh=lh+1 ! No Auxillary Data Pairs
    write(9,pos=4*lh+1) 357.0; lh=lh+1 ! End of Header Marker
    write(9,pos=4*lh+1) 299.0; lh=lh+1 ! Zone Marker
 do n=-nrec,ndata+ndatp
    write(9,pos=4*lh+1) 1; lh=lh+1 ! 1 = Float / 2 = Double
 end do
    write(9,pos=4*lh+1) 0; lh=lh+1 ! No Passive Variables
    write(9,pos=4*lh+1) 0; lh=lh+1 ! No Variable Sharing
    write(9,pos=4*lh+1) -1; lh=lh+1 ! Zero Based Zone Number to Share
 do n=-nrec,ndata+ndatp
    lh=lh+2 ! Minimum Value (Double) of Variables (to be filled)
    lh=lh+2 ! Maximum Value (Double) of Variables (to be filled)
 end do
    lhmb(mb)=lh
 end if
 do mm=0,mbk
    call MPI_BCAST(lhmb(mm),1,MPI_INTEGER,mo(mm),icom,ierr)
 end do
    ns=-nrec; ne=ndata+ndatp; allocate(varmin(ns:ne),varmax(ns:ne))
    lp=lpos(myid)+lhmb(mb)
 do n=ns,ne; lq=(n+nrec)*ltomb
    read(0,rec=n+nrec+1) varr(:)
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
 end do; end do
    varmin(n)=minval(varr(:)); varmax(n)=maxval(varr(:))
 end do
    close(0,status='delete')
 do n=ns,ne
    res=varmin(n); call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MIN,icom,ierr); varmin(n)=fctr
    res=varmax(n); call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr); varmax(n)=fctr
 end do
 if(myid==mo(mb)) then
    l=0; lq=4*(ndata+ndatp+nrec+1)
 do n=ns,ne
    write(9,pos=4*(lp-lq+l)+1) dble(varmin(n)); l=l+2 ! 8-Bytes "Stream"
    write(9,pos=4*(lp-lq+l)+1) dble(varmax(n)); l=l+2 ! 8-Bytes "Stream"
 end do
 end if
    close(9)
 !end if

!===== END OF JOB

 if(myid==0) then
    write(*,*) "Finished."
 end if
    call MPI_FINALIZE(ierr)

 end program main3d

!*****
