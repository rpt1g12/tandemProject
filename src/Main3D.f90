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
 real(k8) :: xpos,ypos,zpos
 integer(k4) :: lsignal,idsignal

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))

	ll=max(npro,12); allocate(ista(MPI_STATUS_SIZE,ll),ireq(ll))

	inquire(iolength=ll) pi; nrec=ll/2

!===== INPUT PARAMETERS

    open(9,file='inputo.dat',shared)
    read(9,*) cinput,mbk,bkx,bky,bkz
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
    read(9,*) cinput,tgustd,tguste
    read(9,*) cinput,aoa,talphas,talphar
    read(9,*) cinput,LES,smago1,smago2
    read(9,*) cinput,output,ogrid,osol,oblock
    close(9)

    cinput=cinput; fltk=pi*fltk; fltkbc=pi*fltkbc
    rhooo=1; poo=1/gam; aoo=sqrt(gam*poo/rhooo); amachoo=sqrt(amach1**2+amach2**2+amach3**2)
    srefoo=111.0_k8/tempoo; srefp1dre=(srefoo+1)/reoo; sqrtrema=sqrt(reoo*amachoo); sqrtremai=1/sqrtrema
    uoo(1)=amach1*aoo; uoo(2)=amach2*aoo; uoo(3)=amach3*aoo
    ! rpt-Initialising the record count 
    nwrec=0
    ! rpt-Do not use postprocessing subroutines
    ispost=.false.
    ! rpt-Position of signal sampling
    !xpos=-1.0_k8;ypos=0.01_k8;zpos=0.005_k8

	abc(:,0)=(/a01,a02,a03,a04,a05,a06/)
	abc(:,1)=(/a10,a12,a13,a14,a15,a16/)
	abc(:,2)=(/a20,a21,a23,a24,a25,a26/)

	ll=3+5*(ndata+1)
    allocate(times(0:ndata),cfilet(-1:ndata),ctecplt(-1:ndata),varm(0:1,0:mpro),varmin(ll),varmax(ll))
    allocate(lximb(0:mbk),letmb(0:mbk),lzemb(0:mbk),lhmb(0:mbk),mo(0:mbk),npc(0:mbk,3))
	allocate(czonet(0:mbk),cthead(0:mbk))

    call inputext

    ! rpt-Forcing parameters
    xfor=-(0.5_k8+3*cos(aoa*pi/180_k8))!cos(delt1)-0.5_k8-1.0_k8+(0.1_k8);
    yfor=-3*sin(aoa*pi/180_k8)!-sin(delt1)+(0.129_k8);
    rfor=5.0e-1
    amfor=amfor*amachoo/100.0e0
    tsfor=151.751e0;tefor=200.000e0
!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

    mo(0)=0
 do mm=1,mbk
    mo(mm)=mo(mm-1)+npc(mm-1,1)*npc(mm-1,2)*npc(mm-1,3)
 end do
 do mm=0,mbk
 if(myid>=mo(mm)) then; mb=mm; end if
 end do
    lxio=lximb(mb); leto=letmb(mb); lzeo=lzemb(mb)

    if (output==2) then
    cfilet(-1)='grid'
 do n=0,ndata
    no(2)=n/100; no(1)=mod(n,100)/10; no(0)=mod(n,10)
    cno=achar(no+48); cfilet(n)='n'//cno(2)//cno(1)//cno(0)
 end do
 do n=-1,ndata
	ctecplt(n)='data/'//cfilet(n)//'.plt'
 end do
 do mm=0,mbk
    no(2)=mm/100; no(1)=mod(mm,100)/10; no(0)=mod(mm,10)
    cno=achar(no+48); czonet(mm)='z'//cno(2)//cno(1)//cno(0)
	cthead(mm)='data/'//czonet(mm)//'.plt'
 end do
    end if
    no(2)=mb/100; no(1)=mod(mb,100)/10; no(0)=mod(mb,10)
    cno=achar(no+48); cnzone=cno(2)//cno(1)//cno(0)
    czone='zone'//cnzone;
    coutput='out/output'//cnzone//'.plt'
    cgrid='misc/grid'//cnzone//'.dat';
    crestart='rsta/restart'//cnzone//'.dat'
    cpostdat='data/postdat'//cnzone//'.dat'

    no(4)=myid/10000; no(3)=mod(myid,10000)/1000;
    no(2)=mod(myid,1000)/100; no(1)=mod(myid,100)/10; no(0)=mod(myid,10)
    cno=achar(no+48); cnnode=cno(4)//cno(3)//cno(2)//cno(1)//cno(0)
    cdata='data/data'//cnnode//'.dat'; cturb='misc/turb'//cnnode//'.dat'

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

    !allocate(qo(0:lmx,5),qa(0:lmx,5),qb(0:lmx,5),de(0:lmx,5))
    !qb(:,:)=0
    ! RPT-take out qb
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

    open(9,file=cgrid,access='stream',form='unformatted',shared)
    lp=lpos(myid)
 do nn=1,3; lq=(nn-1)*ltomb
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    read(9,pos=8*(lp+lq+lio(j,k))+1) ss(l:l+lxi,nn)
 end do; end do
 end do
    close(9)
    call MPI_BARRIER(icom,ierr)
 if(myid==mo(mb)) then
    open(9,file=cgrid); close(9,status='delete')
 end if

    !RPT-FIND POSITION FOR SIGNAL SAMPLING
    !idsignal=-1
    !varr(:)=sqrt((ss(:,1)-xpos)**2+(ss(:,2)-ypos)**2+(ss(:,3)-zpos)**2)
    !lsignal=minloc(varr(:),1)-1
    !ra0=varr(lsignal)
    !CALL MPI_ALLREDUCE(ra0,ra1,1,MPI_REAL8,MPI_MIN,icom,ierr)
    !if (abs(ra0-ra1)/ra0<sml) then
    !   idsignal=myid
    !end if
    !if (myid==idsignal) then
    !open(6,file='data/signal.dat')
    !end if


    rr(:,1)=ss(:,1)
    m=1; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    qo(:,1)=rr(:,1); qo(:,2)=rr(:,2); qo(:,3)=rr(:,3)

    rr(:,1)=ss(:,2)
    m=2; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    qa(:,1)=rr(:,1); qa(:,2)=rr(:,2); qa(:,3)=rr(:,3)

    rr(:,1)=ss(:,3)
    m=3; call mpigo(ntdrv,nrone,n45go,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
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
!    m=1; call mpigo(ntdrv,nrall,n45go,m); call deriv(3,3,m); call deriv(2,2,m); xim(:,m)=rr(:,3)-rr(:,2)
!    rr(:,3)=de(:,2)*ss(:,1); rr(:,2)=de(:,3)*ss(:,1)
!    m=2; call mpigo(ntdrv,nrall,n45go,m); call deriv(3,3,m); call deriv(2,2,m); xim(:,m)=rr(:,3)-rr(:,2)
!    rr(:,3)=qo(:,2)*ss(:,2); rr(:,2)=qo(:,3)*ss(:,2)
!    m=3; call mpigo(ntdrv,nrall,n45go,m); call deriv(3,3,m); call deriv(2,2,m); xim(:,m)=rr(:,3)-rr(:,2)
!
!    rr(:,1)=qa(:,3)*ss(:,3); rr(:,3)=qa(:,1)*ss(:,3)
!    m=1; call mpigo(ntdrv,nrall,n45go,m); call deriv(1,1,m); call deriv(3,3,m); etm(:,m)=rr(:,1)-rr(:,3)
!    rr(:,1)=de(:,3)*ss(:,1); rr(:,3)=de(:,1)*ss(:,1)
!    m=2; call mpigo(ntdrv,nrall,n45go,m); call deriv(1,1,m); call deriv(3,3,m); etm(:,m)=rr(:,1)-rr(:,3)
!    rr(:,1)=qo(:,3)*ss(:,2); rr(:,3)=qo(:,1)*ss(:,2)
!    m=3; call mpigo(ntdrv,nrall,n45go,m); call deriv(1,1,m); call deriv(3,3,m); etm(:,m)=rr(:,1)-rr(:,3)
!
!    rr(:,2)=qa(:,1)*ss(:,3); rr(:,1)=qa(:,2)*ss(:,3)
!    m=1; call mpigo(ntdrv,nrall,n45go,m); call deriv(2,2,m); call deriv(1,1,m); zem(:,m)=rr(:,2)-rr(:,1)
!    rr(:,2)=de(:,1)*ss(:,1); rr(:,1)=de(:,2)*ss(:,1)
!    m=2; call mpigo(ntdrv,nrall,n45go,m); call deriv(2,2,m); call deriv(1,1,m); zem(:,m)=rr(:,2)-rr(:,1)
!    rr(:,2)=qo(:,1)*ss(:,2); rr(:,1)=qo(:,2)*ss(:,2)
!    m=3; call mpigo(ntdrv,nrall,n45go,m); call deriv(2,2,m); call deriv(1,1,m); zem(:,m)=rr(:,2)-rr(:,1)

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
	fctr=pi/(mbci+1)
 do i=1,mbci
	sbci(i)=half*sbci(i)*(1+cos(i*fctr))
 end do
    ll=-1; rr(:,1)=0
 do nn=1,3; do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn); iq=1-2*ip
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    ll=ll+1; res=1/yaco(l); rr(l,1)=rr(l,1)+1; rr(ll,2)=res; rr(ll,3)=l+sml
 do ii=1,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; rr(l,1)=rr(l,1)+1; rr(ll,2)=res*sbci(ii); rr(ll,3)=l+sml
 end do
 end do
 end do
 end if
 end do; end do
    lp=ll; allocate(sbcc(0:lp))
 do ll=0,lp; l=rr(ll,3)
    sbcc(ll)=rr(ll,2)/rr(l,1)
 end do

!===== SETTING UP OUTPUT FILE & STORING GRID DATA

 selectcase(output)
 case(2)
 if(myid==0) then
 do n=-1,ndata
	open(0,file=ctecplt(n)); close(0,status='delete')
 end do
 end if
    open(0,file=cdata,access='direct',form='unformatted',recl=nrec*(lmx+1))
 do nn=1,3
    varr(:)=ss(:,nn); write(0,rec=nn) varr(:); call vminmax(nn)
 end do
 case(0)
       inquire(iolength=lp) varr
       open(0,file=cdata,access='direct',recl=lp)
    if ((1-nto)*nts*nto==0) then
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
          nwrec=nwrec+1
             varr(:)=-1/yaco(:); write(0,rec=nwrec) varr(:)
          do nn=2,3
          nwrec=nwrec+1
             varr(:)=zem(:,nn); write(0,rec=nwrec) varr(:)
          end do
       end if
    ngrec=nwrec
    else
    ngrec=3+9*ngridv
    if (nto==2) then
       nwrec=iwrec
    end if
    end if
 case(1)
   call plot3d(gflag=ogrid,sflag=0,bflag=oblock)
 end select


!===== SETTING UP SPONGE ZONE PARAMETERS

    call spongeup

!===== SETTING UP FORCING PARAMETERS

    if (forcing==1) then
    call forceup
    end if

!===== INITIAL CONDITIONS

 if(nts==0) then
    n=0; ndt=0; dt=0; dts=0; dte=0; timo=0
    call initialo
 else
    open(9,file=crestart,access='stream',form='unformatted',shared); lh=0
    read(9,pos=k8*lh+1) n; lh=lh+1
    read(9,pos=k8*lh+1) ndt; lh=lh+1
    read(9,pos=k8*lh+1) dt; lh=lh+1
    read(9,pos=k8*lh+1) dts; lh=lh+1
    read(9,pos=k8*lh+1) dte; lh=lh+1
    read(9,pos=k8*lh+1) timo; lh=lh+1
    lp=lpos(myid)+lh
    if (tsam<timo) then
       tsam=timo
    end if
 do m=1,5; lq=(m-1)*ltomb
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    read(9,pos=k8*(lp+lq+lio(j,k))+1) qa(l:l+lxi,m)
 end do; end do
 end do
    close(9)
 end if
    !qb(:,:)=0

!============================================
!===== BEGINNING OF TIME MARCHING IN SOLUTION
!============================================

    wts=MPI_WTIME()

 if(myid==0) then
    open(1,file='signal.dat'); close(1,status='delete')
 end if
    call MPI_BARRIER(icom,ierr)
    open(1,file='signal.dat',access='direct',form='formatted',recl=16,shared)

     if (myid==0) then
     write(*,"(3x,'n',8x,'time',9x,'Cl',9x,'Cd',5x)")  
     write(*,"('============================================')")
     end if
    ndati=-1; nsigi=-1; dtsum=0
    if ((nto==2).and.(output==0)) then
       ndati=0
    end if
 do while(timo-tmax<0.and.(dt/=0.or.n<=2))

  if (mod(n,nscrn)==0) then
  if (n.ge.2) then
    call clpost(1,ndati) 
  else
    cl=0
  end if
  end if

  if(myid==0.and.mod(n,nscrn)==0) then
     !Change ra0 to the angle of attack needed!!!
     ra0=aoa*pi/180;ra1=cos(ra0);ra2=sin(ra0)
     write(*,"(i8,f12.5,f12.7,f12.7)") &
     n,timo,cl(1,2)*ra1-cl(1,1)*ra2,cl(1,2)*ra2+cl(1,1)*ra1
  end if

!----- FILTERING

 do m=1,5
    rr(:,1)=qa(:,m)
    call mpigo(ntflt,nrone,n45no,m); call filte(1,1); call filte(2,1); call filte(3,1)
    qa(:,m)=rr(:,1)
 end do

!-------------------------------------
!----- BEGINNING OF RUNGE-KUTTA STAGES
!-------------------------------------

    qo(:,:)=qa(:,:)

 do nk=1,nkrk

!----- MOVING FRAME VELOCITY & ACCELERATION BEFORE TIME ADVANCING

    dtko=min(max(nk-2,0),1)*dt/(nkrk-nk+3); dtk=min(nk-1,1)*dt/(nkrk-nk+2)
    call movef(dtko,dtk)

!----- TEMPORARY STORAGE OF PRIMITIVE VARIABLES & PRESSURE

    de(:,1)=1/qa(:,1)
    de(:,2)=qa(:,2)*de(:,1)
    de(:,3)=qa(:,3)*de(:,1)
    de(:,4)=qa(:,4)*de(:,1)

    p(:)=gamm1*(qa(:,5)-half*(qa(:,2)*de(:,2)+qa(:,3)*de(:,3)+qa(:,4)*de(:,4)))
    de(:,5)=gam*p(:)*de(:,1) ! Temperature
    ss(:,1)=srefp1dre*de(:,5)**1.5_k8/(de(:,5)+srefoo) ! Sutherland's Law

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
    dt=dts+(dte-dts)*sin(0.05_k8*pi*(n-ndt))**2

    nout=0; res=tsam+(ndati+1)*(tmax-tsam)/ndata
 if((timo-res)*(timo+dt-res)<=0) then
    nout=1; ndati=ndati+1
 end if
 end if

!----- VISCOUS SHEAR STRESSES & HEAT FLUXES

 if(nviscous==1) then
    de(:,1)=ss(:,1)

    rr(:,1)=de(:,2)
    m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    txx(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    hzz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    tzx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,3)
    m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    txy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    tyy(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    hxx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,4)
    m=4; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    hyy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    tyz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    tzz(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    rr(:,1)=de(:,5)
    m=5; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1,m); call deriv(2,1,m); call deriv(1,1,m)
    ss(:,1)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
    ss(:,2)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
    ss(:,3)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

    fctr=2.0_k8/3
    rr(:,1)=de(:,1)*yaco(:)
    rr(:,2)=gamm1prndtli*rr(:,1)
    !qb(:,3)=de(:,1)

    selectcase(LES)
    case(1)
    de(:,1)=(txx(:)*txx(:)+tyy(:)*tyy(:)+tzz(:)*tzz(:)+& !rpt- SijSij
             (hzz(:)+txy(:))*(hzz(:)+txy(:))+&
             (hyy(:)+tzx(:))*(hyy(:)+tzx(:))+&
             (hxx(:)+tyz(:))*(hxx(:)+tyz(:)))
    varr(:)=(-1/yaco(:))**1.5 ! rpt- Volume
    rr(:,3)=qa(:,1)*smago1**2*varr(:)*sqrt(2*(de(:,1))) ! rpt-nuSGS
    !qb(:,2)=rr(:,3)
    !qb(:,4)=qb(:,2)/qb(:,3)
    rr(:,1)=rr(:,1)+rr(:,3)*yaco(:)
    rr(:,2)=rr(:,2)+tgamm1prndtli*rr(:,3)*yaco(:)   
    rr(:,3)=fctr*(qa(:,1)*smago2*varr(:)*de(:,1)) ! rpt-2/3*ro*kSGS
    de(:,5)=fctr*(txx(:)+tyy(:)+tzz(:))

    txx(:)=rr(:,1)*(2*txx(:)-de(:,5))-rr(:,3)
    tyy(:)=rr(:,1)*(2*tyy(:)-de(:,5))-rr(:,3)
    tzz(:)=rr(:,1)*(2*tzz(:)-de(:,5))-rr(:,3)
    txy(:)=rr(:,1)*(txy(:)+hzz(:))
    tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
    tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
    case(0)
    de(:,5)=fctr*(txx(:)+tyy(:)+tzz(:))

    txx(:)=rr(:,1)*(2*txx(:)-de(:,5))
    tyy(:)=rr(:,1)*(2*tyy(:)-de(:,5))
    tzz(:)=rr(:,1)*(2*tzz(:)-de(:,5))
    txy(:)=rr(:,1)*(txy(:)+hzz(:))
    tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
    tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
    end select

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
    m=1; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,2)*ss(:,1)+xim(:,1)*p(:)
    rr(:,2)=qa(:,2)*ss(:,2)+etm(:,1)*p(:)
    rr(:,3)=qa(:,2)*ss(:,3)+zem(:,1)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txx(:)-xim(:,2)*txy(:)-xim(:,3)*tzx(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txx(:)-etm(:,2)*txy(:)-etm(:,3)*tzx(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txx(:)-zem(:,2)*txy(:)-zem(:,3)*tzx(:)
 end if
    m=2; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,3)*ss(:,1)+xim(:,2)*p(:)
    rr(:,2)=qa(:,3)*ss(:,2)+etm(:,2)*p(:)
    rr(:,3)=qa(:,3)*ss(:,3)+zem(:,2)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*txy(:)-xim(:,2)*tyy(:)-xim(:,3)*tyz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*txy(:)-etm(:,2)*tyy(:)-etm(:,3)*tyz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*txy(:)-zem(:,2)*tyy(:)-zem(:,3)*tyz(:)
 end if
    m=3; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

    rr(:,1)=qa(:,4)*ss(:,1)+xim(:,3)*p(:)
    rr(:,2)=qa(:,4)*ss(:,2)+etm(:,3)*p(:)
    rr(:,3)=qa(:,4)*ss(:,3)+zem(:,3)*p(:)
 if(nviscous==1) then
    rr(:,1)=rr(:,1)-xim(:,1)*tzx(:)-xim(:,2)*tyz(:)-xim(:,3)*tzz(:)
    rr(:,2)=rr(:,2)-etm(:,1)*tzx(:)-etm(:,2)*tyz(:)-etm(:,3)*tzz(:)
    rr(:,3)=rr(:,3)-zem(:,1)*tzx(:)-zem(:,2)*tyz(:)-zem(:,3)*tzz(:)
 end if
    m=4; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
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
    m=5; call mpigo(ntdrv,nrall,n45no,m); call deriv(1,1,m); call deriv(2,2,m); call deriv(3,3,m)
    de(:,m)=rr(:,1)+rr(:,2)+rr(:,3)

!----- IMPLEMENTATION OF SPONGE CONDITION

    call spongego

!----- IMPLEMENTATION OF FORCING

     if (forcing==1) then
     call forcego
     end if

!----- PREPARATION FOR GCBC & GCIC

 do nn=1,3
 select case(nn)
 case(1); drva=>drva1; cm=>cm1; case(2); drva=>drva2; cm=>cm2; case(3); drva=>drva3; cm=>cm3
 end select
 do ip=0,1; np=nbc(ip,nn); i=ip*ijk(1,nn)
 if((np-10)*(np-20)*(np-25)*(np-30)==0) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(i,j,k,nn)
    call eleme(l,cm(jk,:,ip)); call xtq2r(cm(jk,:,ip)); drva(jk,:,ip)=matmul(xt(:,:),yaco(l)*de(l,:))
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
 end if

!----- IMPLEMENTATION OF GCBC & GCIC

    ll=-1
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
    if(iq*(vn+vs)>0) then; cha(1:3)=0; end if
    if(iq*(vn+vs+ao)>0) then; cha(4)=0; end if
    if(iq*(vn+vs-ao)>0) then; cha(5)=0; end if
 case(20,25)
    cha(4+ip)=cha(5-ip)+iq*aoi*qa(l,1)*(2*sum(cm(jk,:,ip)*dudtmf(:))+100*(vn+vs))
 case(30)
    if(iq*(vn+vs)>0) then; cha(1:3)=dha(1:3); end if
    if(iq*(vn+vs+ao)>0) then; cha(4)=dha(4); end if
    if(iq*(vn+vs-ao)>0) then; cha(5)=dha(5); end if
 end select
    call xtr2q(cm(jk,:,ip)); dha(:)=matmul(xt(:,:),(cha(:)-drva(jk,:,ip)))
    ll=ll+1; de(l,:)=de(l,:)+sbcc(ll)*dha(:)
 do ii=1,mbci; l=indx3(i+iq*ii,j,k,nn)
    ll=ll+1; de(l,:)=de(l,:)+sbcc(ll)*dha(:)
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
 !if(timo-tsam+(tmax-tsam)/ndata>0) then
 !   dtsum=dtsum+dt; qb(:,:)=qb(:,:)+half*dt*(qo(:,:)+qa(:,:))
 !if(nout==1) then
 !   times(ndati)=timo-half*dtsum
 !if(n==1) then
 !   qb(:,:)=qo(:,:)
 !else
 !   fctr=1/dtsum; qb(:,:)=fctr*qb(:,:)
 !end if
 !   rr(:,1)=1/qb(:,1)
 !do m=1,5
 !select case(m)
 !case(1); varr(:)=qb(:,m); case(2:4); varr(:)=rr(:,1)*qb(:,m)+umf(m-1)
 !case(5); varr(:)=gamm1*(qb(:,m)-half*rr(:,1)*(qb(:,2)*qb(:,2)+qb(:,3)*qb(:,3)+qb(:,4)*qb(:,4)))
 !end select
 !   write(0,rec=5*ndati+m+3) varr(:)
 !end do
 !   dtsum=0; qb(:,:)=0
 !end if
 if(nout==1) then
      times(ndati)=timo
   if(myid==0) then
      write(*,"('===> saving output ',i3,' at time =',f12.8)") ndati,timo
   end if
   selectcase(output)
   case(2)
    rr(:,1)=1/qb(:,1)
 do m=1,5
 select case(m)
 case(1); varr(:)=qb(:,m); case(2:4); varr(:)=rr(:,1)*qb(:,m)+umf(m-1)
 case(5); varr(:)=gam*gamm1*(qb(:,m)-half*rr(:,1)*(qb(:,2)*qb(:,2)+qb(:,3)*qb(:,3)+qb(:,4)*qb(:,4)))
 end select
    nn=3+5*ndati+m; write(0,rec=nn) varr(:); call vminmax(nn)
 end do
   case(1)
      call plot3d(gflag=0,sflag=osol,bflag=oblock)
   case(0)
      !==========SAVING INSTANTANEUS DENSITY, VELOCITY AND PRESSURE
      nwrec=nwrec+1
         varr(:)=qa(:,1); write(0,rec=nwrec) varr(:)
      do m = 2, 4
      nwrec=nwrec+1
        varr(:)=((qa(:,m)/qa(:,1))+umf(m-1)); write(0,rec=nwrec) varr(:)
      end do
      nwrec=nwrec+1
      varr(:)=p(:); write(0,rec=nwrec) varr(:)
      if(myid==0) then
         write(*,"('===>nwrec= ',i8)") nwrec
      end if
   end select


   !===== GENERATING RESTART DATA FILE
   
    if(nrestart==1) then
       call wRestart()
    end if
 end if

 if(timo-tsam>=0.and.mod(n,nsgnl)==0) then
    nsigi=nsigi+1; call signalgo
 end if

 !if (myid==idsignal) then
 !if (timo.le.25.5_k8) then
 !   ra0=qa(lsignal,2)/qa(lsignal,1)+umf(1)
 !   ra1=qa(lsignal,3)/qa(lsignal,1)+umf(2)
 !   ra2=qa(lsignal,4)/qa(lsignal,1)+umf(3)
 !   write(6,"(es15.7,1x,es15.7,1x,es15.7,1x,es15.7)") ra0,ra1,ra2,timo
 !end if
 !end if

!==========================
!===== END OF TIME MARCHING
!==========================

 end do
    close(1)
    if (output==2) then
	close(0)
    end if
    !if (myid==idsignal) then
    !   close(6)
    !end if

    wte=MPI_WTIME(); res=wte-wts
    call MPI_ALLREDUCE(res,wtime,1,MPI_REAL8,MPI_SUM,icom,ierr)
 if(myid==0) then
    if (nto==2) then
       open(9,file='data/timeouts.dat',position='append')
    else
       open(9,file='data/timeouts.dat')
    end if
    write(9,'(es15.7)') times(:)
    close(9)
    open(9,file='walltime.dat',position='append')
    write(9,'(2es15.7)') dble(npro),wtime/npro
    close(9)
 end if


!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

 if(dt==0) then
    if (myid==0) then
       write(*,*) "Overflow."
    end if
    ndata=ndati
 else
    if(myid==0) then
       write(*,'("Simulation time was ",f6.2," hours")') wtime/(3600_k8*npro)
    end if
   selectcase(output)
   case(0)
    if (nto==2) then
       ndata=ndati+(iwrec-ngrec)/5-1
       if (myid==0) then
          write(*,*) ndata
       end if
    end if
      if (output==0) then
      call post(average=.false.)
      end if
   case(2)
    deallocate(qo,qa,qb,de,xim,etm,zem,rr,ss,p,yaco)
 if(nviscous==1) then
    deallocate(txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz)
 end if
	nlmx=(3+5*(ndata+1))*(lmx+1)-1; ll=5*(lmx+1)-1; allocate(vart(0:nlmx),vmean(0:ll))
	open(9,file=cdata,access='direct',form='unformatted',recl=nrec*(nlmx+1))
	read(9,rec=1) vart(:)
	close(9,status='delete')

!----- CALCULATING UNSTEADY FLUCTUATIONS

	fctr=half/(times(ndata)-times(0)); vmean(:)=0
 do n=0,ndata; lis=(3+5*n)*(lmx+1); lie=lis+ll
	if(n==0) then; ra0=fctr*(times(n+1)-times(n)); end if
    if(n==ndata) then; ra0=fctr*(times(n)-times(n-1)); end if
	if(n>0.and.n<ndata) then; ra0=fctr*(times(n+1)-times(n-1)); end if
	vmean(:)=vmean(:)+ra0*vart(lis:lie)
 end do
 do n=0,ndata; lis=(3+5*n)*(lmx+1); lie=lis+ll
	vart(lis:lie)=vart(lis:lie)-vmean(:)
 do m=1,5; nn=3+5*n+m; l=lis+(m-1)*(lmx+1)
	varr(:)=vart(l:l+lmx); call vminmax(nn)
 end do
 end do

!----- COLLECTING DATA FROM SUBDOMAINS & BUILDING TECPLOT OUTPUT FILES

	lje=-1
 do n=-1,ndata
	mq=3+2*min(n+1,1); llmb=mq*ltomb-1; allocate(vara(0:llmb),varb(0:llmb))
	ljs=lje+1; lje=ljs+mq*(lmx+1)-1
 if(myid==mo(mb)) then !===========================================================================
	mps=mo(mb); mpe=mps+npc(mb,1)*npc(mb,2)*npc(mb,3)-1
	lis=0; lie=mq*(lmx+1)-1; vara(lis:lie)=vart(ljs:lje)
 do mp=mps+1,mpe
	lis=lie+1; lie=lis+mq*(lxim(mp)+1)*(letm(mp)+1)*(lzem(mp)+1)-1
    itag=1; call MPI_RECV(vara(lis:lie),lie-lis+1,MPI_REAL4,mp,itag,icom,ista,ierr)
 end do
	lis=0
 do mp=mps,mpe; do m=1,mq; do k=0,lzem(mp); do j=0,letm(mp)
	ljs=lpos(mp)+(m-1)*ltomb+k*(leto+1)*(lxio+1)+j*(lxio+1)
	varb(ljs:ljs+lxim(mp))=vara(lis:lis+lxim(mp)); lis=lis+lxim(mp)+1
 end do; end do; end do; end do
	open(9,file=cthead(mb),access='stream',form='unformatted')
	call techead(9,n,mb,lh)
	deallocate(vara); allocate(vara(0:lh+llmb)); read(9,pos=1) vara(0:lh-1)
	close(9,status='delete')
	lhmb(mb)=lh+llmb+1; vara(lh:lh+llmb)=varb(:)
 if(mb==0) then !----------------------------------------------------------------------------------
 do mm=1,mbk
	itag=2; call MPI_RECV(lhmb(mm),1,MPI_INTEGER8,mo(mm),itag,icom,ista,ierr)
 end do
	llmo=sum(lhmb(:))-1; deallocate(varb); allocate(varb(0:llmo))
	lis=0; lie=lhmb(mb)-1; varb(lis:lie)=vara(:)
 do mm=1,mbk
	lis=lie+1; lie=lis+lhmb(mm)-1
	itag=3; call MPI_RECV(varb(lis:lie),lie-lis+1,MPI_REAL4,mo(mm),itag,icom,ista,ierr)
 end do
 	open(0,file=ctecplt(n),access='direct',form='unformatted',recl=nrec*(llmo+1))
	write(0,rec=1) varb(:)
	close(0)
 else !--------------------------------------------------------------------------------------------
	itag=2; call MPI_SEND(lhmb(mb),1,MPI_INTEGER8,mo(0),itag,icom,ierr)
	itag=3; call MPI_SEND(vara(:),lhmb(mb),MPI_REAL4,mo(0),itag,icom,ierr)
 end if !------------------------------------------------------------------------------------------
 else !============================================================================================
	itag=1; call MPI_SEND(vart(ljs:lje),lje-ljs+1,MPI_REAL4,mo(mb),itag,icom,ierr)
 end if !==========================================================================================
	deallocate(vara,varb)
 end do

!-----

   end select
 end if

if (myid==0) then
   open(9,file='data/post.dat')
   write(9,*) 'ngridv ',ngridv
   write(9,*) 'ndata  ',ndata
   write(9,*) 'ngrec   ',ngrec
   write(9,*) 'nwrec  ',nwrec
   write(9,*) 'lhmb   ',lhmb(mb)
   close(9)
end if

!===== END OF JOB

 if(myid==0) then
    write(*,*) "Finished."
 end if

    call MPI_FINALIZE(ierr)

 end program main3d

!*****
