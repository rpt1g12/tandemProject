!*****
!***** 3D PARALLEL SOLVER
!*****

 program mainpost

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 use rpt
 implicit none
 integer :: lsta
 logical :: tecplot

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

call setup
tecplot=.false.

!===== GRID INPUT & CALCULATION OF GRID METRICS

!===== INPUT RECORDING PARAMETERS
    open(9,file='data/post.dat',shared)
    read(9,*) cinput,ngridv
    read(9,*) cinput,nvarout
    read(9,*) cinput,ndata
    read(9,*) cinput,ndatp
    read(9,*) cinput,nrec
    read(9,*) cinput,nwrec
    read(9,*) cinput
    do n = 0, ndata
    read(9,*) times(n)
    end do
    close(9)

 select case(nvarout)
 case(0,1,2,3,4,5,6)
   totVar=5
 case(7)
   totVar=14
 end select
 if (tecplot) then
 lsta=17+ngridv*45+(ndata+1)*25+ndatp*10+24+nwrec+3+4*nwrec
 end if
!===== INQUIRE RECORD LENGTH
 open(3,file=cpostdat,access='stream',shared)
 open(9,file=coutput,access='stream',shared)
 nread=0

!===== READ X,Y,Z COORDINATES
 do nn=1,3
    nread=nread+1
    if (tecplot) then
    call tpostread(nread,lsta)
    else 
    call postread(nread)
    end if
    ss(:,nn)=varr(:)
 end do
!===== READ METRICS
 if (ngridv==1) then
 do nn=1,3
    nread=nread+1
    if (tecplot) then
    call tpostread(nread,lsta)
    else 
    call postread(nread)
    end if
    xim(:,nn)=varr(:)
 end do
 do nn=1,3
    nread=nread+1
    if (tecplot) then
    call tpostread(nread,lsta)
    else 
    call postread(nread)
    end if
    etm(:,nn)=varr(:)
 end do
 do nn=1,3
    nread=nread+1
    if (tecplot) then
    call tpostread(nread,lsta)
    else 
    call postread(nread)
    end if
    zem(:,nn)=varr(:)
 end do
 end if

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

    do i = 1, 3
       xyz(:,i)=ss(:,i)
    end do


!===== COMPUTE METRICS IF NOT AVAILABLE YET
    if(ngridv.ne.1) then
    xim(:,1)=qa(:,2)*de(:,3)-de(:,2)*qa(:,3)
    xim(:,2)=de(:,2)*qo(:,3)-qo(:,2)*de(:,3)
    xim(:,3)=qo(:,2)*qa(:,3)-qa(:,2)*qo(:,3)
    etm(:,1)=qa(:,3)*de(:,1)-de(:,3)*qa(:,1)
    etm(:,2)=de(:,3)*qo(:,1)-qo(:,3)*de(:,1)
    etm(:,3)=qo(:,3)*qa(:,1)-qa(:,3)*qo(:,1)
    zem(:,1)=qa(:,1)*de(:,2)-de(:,1)*qa(:,2)
    zem(:,2)=de(:,1)*qo(:,2)-qo(:,1)*de(:,2)
    zem(:,3)=qo(:,1)*qa(:,2)-qa(:,1)*qo(:,2)
    end if

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

!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
 if (ndatp.ne.1) then
 ltz=-1
 call average
 ndatp=1
 end if

!===== Save average values into array
 nread=nwrec-totVar
 do nn=1,5
    nread=nread+1
    if (tecplot) then
    call tpostread(nread,lsta)
    else 
    call postread(nread)
    end if
    qa(:,nn)=varr(:)
 end do

call clpost(1,2,1,.true.)

if (myid==11) then
varr(:)=qa(:,5)
call spcavg(2,3,0,1,2)
end if

if (tecplot) then
   call tdatawrite(lsta)
else
   call datawrite
end if
close(3)
close(9)

!===== WRITE TECPLOT FILE
 call post

!===== END OF JOB
if (myid==11) then
 write(*,*) sum(area)
end if

 if(myid==0) then
    write(*,*) "Finished."
 end if


    call MPI_FINALIZE(ierr)


 contains


subroutine setup
   
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
    cpostdat='data/postdat'//cno(2)//cno(1)//cno(0)//'.dat'
    no(4)=myid/10000; no(3)=mod(myid,10000)/1000; no(2)=mod(myid,1000)/100
    no(1)=mod(myid,100)/10; no(0)=mod(myid,10); cno=achar(no+48)
    cdata='data/data'//cno(4)//cno(3)//cno(2)//cno(1)//cno(0)//'.dat'
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
    allocate(xyz(0:lmx,3))
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

    allocate(lio(0:let,0:lze))
 do k=0,lze; kp=k*(leto+1)*(lxio+1)
 do j=0,let; jp=j*(lxio+1)
    lio(j,k)=jp+kp
 end do
 end do

end subroutine setup

!=====AVERAGE RESULTS
 subroutine average
 implicit none
 integer :: totVar
 real(nr),dimension(:),allocatable :: delt

 select case(nvarout)
 case(0,1,2,3,4,5,6)
   totVar=5
 case(7)
   totVar=14
 end select
 
    ns=0; ne=ndata; allocate(delt(ns:ne))
    fctr=half/(times(ne)-times(ns))
    delt(ns)=fctr*(times(ns+1)-times(ns)); delt(ne)=fctr*(times(ne)-times(ne-1))
 do n=ns+1,ne-1
    delt(n)=fctr*(times(n+1)-times(n-1))
 end do
 do m = 1, totVar
    rr(:,1)=0
 do n=0,ndata
    call postread((n*totVar+nrec+m))
    rr(:,1)=rr(:,1)+delt(n)*varr(:)
 end do
    varr(:)=rr(:,1)
    nwrec=nwrec+1
    call postwrite(nwrec)
 end do

 end subroutine average

 ! ---SPACE AVERAGING
 subroutine spcavg(plane,dir,pos,nvar,nout)
    implicit none
    integer, intent(in) :: plane,dir,pos,nvar,nout
    integer :: odir,idir,loop,axis
    real(nr), dimension(:,:), allocatable :: avg,delt
    integer, dimension(:,:), allocatable :: lavg

    select case(plane)
    case(1)
     select case(dir)
     case(2); odir=ijk(3,1); idir=ijk(2,1); loop=1; axis=3
     case(3); odir=ijk(2,1); idir=ijk(3,1); loop=2; axis=2
     end select                                  
    case(2)                                     
     select case(dir)                          
     case(1); odir=ijk(2,2); idir=ijk(3,2); loop=2; axis=3
     case(3); odir=ijk(3,2); idir=ijk(2,2); loop=1; axis=1 
     end select                               
    case(3)                                  
     select case(dir)                       
     case(1); odir=ijk(3,3); idir=ijk(2,3); loop=1; axis=2 
     case(2); odir=ijk(2,3); idir=ijk(3,3); loop=2; axis=1 
     end select
    end select

    allocate(lavg(0:idir,0:odir))
    allocate(avg(0:odir,2),delt(0:idir,0:odir))
    open(7,file='out/avg.dat')
    write(7,"('VARIABLES= dir, var')")
    write(7,"('ZONE T= avg')")

    !! READ VARIABLE TO AVERAGE
    !nread=nrec+(totVar*nvar)+nout
    !write(*,*) nread
    !call postread(nread)
    
    select case(loop)
    case(1);
    do k = 0, odir
       do j = 0, idir 
          lavg(j,k)=indx3(pos,j,k,plane)
       end do
    end do
    case(2);
    do k = 0, odir
       do j = 0, idir
          lavg(j,k)=indx3(pos,j,k,plane)
       end do
    end do
    end select

    do k = 0, odir
    ns=lavg(0,k); ne=lavg(idir,k)
    fctr=half/abs((xyz(ne,dir)-xyz(ns,dir)))
    delt(0,k)=fctr*abs((xyz(lavg(1,k),dir)-xyz(ns,dir)))
    delt(idir,k)=fctr*abs((xyz(ne,dir)-xyz(lavg(idir-1,k),dir)))
    avg(k,2)=delt(0,k)*varr(ns)+delt(idir,k)*varr(ne)
    avg(k,1)=xyz(ns,axis)
       do j = 1, idir-1; l=lavg(j,k)
         delt(j,k)=fctr*abs((xyz(lavg(j+1,k),dir)-xyz(lavg(j-1,k),dir)))
         avg(k,2)=avg(k,2)+delt(j,k)*varr(l)
       end do
    write(7,'(f10.5,"   ",f10.5)') avg(k,1),avg(k,2)
    end do
    close(7)
 end subroutine spcavg
!*****
end program mainpost
