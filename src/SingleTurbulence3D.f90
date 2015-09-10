!*****
!***** 3D SINGLE AEROFOIL-TURBULENCE INTERACTION
!*****

 module problemcase

 use mpi
 use subroutineso
 use gridgen
 implicit none

 integer(k4) :: nbody,nthick,ngridv,nito,nits,litr,lsz,ltz,ntz
 integer(k4),dimension(:),allocatable :: iit,idsgnl,lsgnl
 integer(k4),dimension(4) :: mjct
 real(k8),dimension(:,:,:),allocatable :: vito
 real(k8),dimension(:,:),allocatable :: vit
 real(k8),dimension(:),allocatable :: tt
 real(k8),dimension(0:1,3) :: szr,szp
 real(k8),dimension(3,0:1) :: tam
 real(k8),dimension(0:1) :: tl0,tl1,tlw
 real(k8) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,szco
 ! rpt-tandem variables
 real(k8) :: gap,c1,c2,delt1,delt2
 ! rpt-grid space modifiers for grid generation
 real(k8) :: ximod,etamod
 real(k8) :: tgusto,eps,ck1,ck2,ck3,amp1,amp2,amp3,vk1,vk2,slit,gaus,cfit,tla,tlb,cutlb
 real(k8) :: denxit

 contains

!===== SUBINPUT PARAMETERS

 subroutine inputext

 integer(k4), dimension(:) :: ilximb(0:bkx-1),iletmb(0:bky-1),ilzemb(0:bkz-1)
 integer(k4), dimension(:) :: npcx(0:bkx-1),npcy(0:bky-1),npcz(0:bkz-1)

    open(9,file='inputp.dat',shared)
    ! rpt-two extra blocks added
    read(9,*) cinput,ilximb(0:bkx-1)
    read(9,*) cinput,iletmb(0:bky-1)
    read(9,*) cinput,ilzemb(0:bkz-1)
    read(9,*) cinput,nbody,nthick
    read(9,*) cinput,ngridv
    read(9,*) cinput,npcx(0:bkx-1)
    read(9,*) cinput,npcy(0:bky-1)
    read(9,*) cinput,npcz(0:bkz-1)
    read(9,*) cinput,nito,nits,litr
    read(9,*) cinput,smgrid,domlen
    read(9,*) cinput,ximod,etamod
    read(9,*) cinput,span,wlew,wlea
    read(9,*) cinput,gap,c1,c2
    read(9,*) cinput,delt1,delt2
    read(9,*) cinput,szth1,szth2,szxt,szco
    read(9,*) cinput,tgusto
    read(9,*) cinput,ck1,ck2,ck3
    read(9,*) cinput,amp2,amp3
    read(9,*) cinput,vk1,vk2
    read(9,*) cinput,slit,gaus,cfit
    read(9,*) cinput,tl0(0),tl1(0),tlw(0)
    read(9,*) cinput,tl0(1),tl1(1),tlw(1)
    read(9,*) cinput,tam(:,0)
    read(9,*) cinput,tam(:,1)
    read(9,*) cinput,cutlb
    close(9)

    ck1=2*pi*ck1; ck2=2*pi*ck2; ck3=2*pi*ck3
    amp1=-(amp2*ck2+amp3*ck3)/ck1
    eps=maxval((/abs(amp2),sum(tam)/6,sml/))
    tla=minval((/tl0(:),tl1(:)/)); tlb=maxval((/tl0(:),tl1(:)/))
    ! rpt-from degrees to radians
    delt1=delt1*pi/180.0e0; delt2=delt2*pi/180.0e0

    tsam=max(tmax-2*pi/max(ck1*uoo(1)+ck2*uoo(2)+ck3*uoo(3),sml),tsam)

    do k = 0, bkz-1
       do j = 0, bky-1
          do i = 0, bkx-1; l=k*(bkx*bky)+j*bkx+i
             lximb(l)=ilximb(i)
             letmb(l)=iletmb(j)
             lzemb(l)=ilzemb(k)
             
             npc(l,1)=npcx(i)
             npc(l,2)=npcy(j)
             npc(l,3)=npcz(k)
          end do
       end do
    end do

    ! Dependent on the grid topology!!
    lxi0=lximb(0);lxi1=lximb(1);lxi2=lximb(2)
    let0=letmb(0);let1=letmb(bkx)
    lze0=lzemb(0)

    allocate(mxc(nits),ran(nits,3),sit(nits,3),ait(nits,3),xit(nits),yit(nits),zit(nits))
    allocate(iit(0:lze0),idsgnl(0:lze0),lsgnl(0:lze0))

    open(9,file='randnum.dat',shared)
 do m=1,3
    read(9,*) ran(:,m); read(9,*) ait(:,m)
 end do
!    read(9,*) xit(:); read(9,*) yit(:); read(9,*) zit(:)
    read(9,*) yit(:); read(9,*) zit(:)
    close(9)
 do nn=1,nits
	if(ran(nn,1)-gaus<=0) then; mxc(nn)=0; else; mxc(nn)=1; end if
 end do
    res=(2*(2-cutlb)*tlb*span*slit/nits)**(one/3)
 do m=1,3; do nn=1,nits
    ran(nn,m)=tl0(mxc(nn))+(tl1(mxc(nn))-tl0(mxc(nn)))*ran(nn,m)**tlw(mxc(nn))
    sit(nn,m)=1/ran(nn,m)**2
    ait(nn,m)=res*tam(m,mxc(nn))*(2*ait(nn,m)-1)*sit(nn,m)**1.5_k8
 end do; end do
    denxit=slit/(2*sum(ran(:,1)))
    xit(1)=-domlen+half*(szth1+slit)-denxit*ran(1,1)
 do nn=2,nits
    xit(nn)=xit(nn-1)-denxit*(ran(nn-1,1)+ran(nn,1))
 end do
 do nn=1,nits
    res=3*span-2*max(ran(nn,1),ran(nn,2),ran(nn,3))
!    xit(nn)=-domlen+half*(szth1+slit)-slit*xit(nn)
    yit(nn)=(1-cutlb)*ran(nn,2)*(2*yit(nn)-1)
    zit(nn)=min(span,res)*(zit(nn)-half)
 end do

 end subroutine inputext

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

 subroutine domdcomp

 integer(k4) :: ip1,jp1,ip2,jp2

    !ip=30*nthick+35*(1-nthick); jp=35*(1-nbody)+nbody*(20+5*nviscous)

    select case(nthick)
    case(0); ip1=35; ip2=35; ip=35; jp=35
    case(1); ip1=30; ip2=30; ip=30; jp=30
    case(2); ip1=30; ip2=30; ip=30; jp=30
    case(3); ip1=30; ip2=30; ip=30; jp=30
    end select

    select case(nbody)
    case(0); jp1=35; jp2=35
    case(1); jp1=20+5*nviscous; jp2=20+5*nviscous
    case(2); jp1=20+5*nviscous; jp2=30
    case(3); jp1=30; jp2=20+5*nviscous
    end select

 select case(mb)
 case(0);  nbcs(1)=10;  nbce(1)=ip1; nbcs(2)=10;  nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(1);  nbcs(1)=ip1; nbce(1)=ip1; nbcs(2)=10;  nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(2);  nbcs(1)=ip1; nbce(1)=10;  nbcs(2)=10;  nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(3);  nbcs(1)=10;  nbce(1)=ip1; nbcs(2)=ip1; nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(4);  nbcs(1)=ip1; nbce(1)=ip1; nbcs(2)=ip1; nbce(2)=jp1; nbcs(3)=45; nbce(3)=45
 case(5);  nbcs(1)=ip1; nbce(1)=10;  nbcs(2)=ip1; nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(6);  nbcs(1)=10;  nbce(1)=ip1; nbcs(2)=ip1; nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(7);  nbcs(1)=ip1; nbce(1)=ip1; nbcs(2)=jp1; nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(8);  nbcs(1)=ip1; nbce(1)=10;  nbcs(2)=ip1; nbce(2)=ip1; nbcs(3)=45; nbce(3)=45
 case(9);  nbcs(1)=10;  nbce(1)=ip1; nbcs(2)=ip1; nbce(2)=10;  nbcs(3)=45; nbce(3)=45
 case(10); nbcs(1)=ip1; nbce(1)=ip1; nbcs(2)=ip1; nbce(2)=10;  nbcs(3)=45; nbce(3)=45
 case(11); nbcs(1)=ip1; nbce(1)=10;  nbcs(2)=ip1; nbce(2)=10;  nbcs(3)=45; nbce(3)=45
 end select

 ! rpt-setting the neighbouring blocks
 select case(mb)
 case(0,3,6,9); ms(1)=mb+2; me(1)=mb+1;
 case(1,4,7,10); ms(1)=mb-1; me(1)=mb+1;
 case(2,5,8,11); ms(1)=mb-1; me(1)=mb-2
 end select
 select case(mb)
 case(0,1,2); ms(2)=mb+9; me(2)=mb+3;
 case(3,4,5,6,7,8); ms(2)=mb-3; me(2)=mb+3;
 case(9,10,11); ms(2)=mb-3; me(2)=mb-9
 end select
    ms(3)=mb; me(3)=mb

 end subroutine domdcomp

!===== GRID GENERATION

 subroutine makegrid

    ! rpt-new grid generation file, see GridTandemAerofoil.90
    ! rpt-there are new arguments
    call gridaerofoil(ngridv,nthick,litr,smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,cutlb,c1,delt1,ximod,etamod)

 end subroutine makegrid

!===== SETTING UP SPONGE ZONE PARAMETERS

 subroutine spongeup

 do nn=1,3
 ! rpt-set nsz=1 if contains boundary region
 if(nbcs(nn)==10) then; nsz(0,nn)=1; else; nsz(0,nn)=0; end if
 if(nbce(nn)==10) then; nsz(1,nn)=1; else; nsz(1,nn)=0; end if
 select case(nn)
 case(1); szr(0,nn)=1/szth1; szp(0,nn)=-domlen+szth1; szr(1,nn)=1/(szth1+szxt)
 ! rpt-the sponge zone lenght changes due to just one element
 szp(1,nn)=domlen-szth1
 case(2); szr(0,nn)=1/szth2; szp(0,nn)=-domlen+szth2; szr(1,nn)=1/szth2; szp(1,nn)=domlen-szth2
 case(3); szr(0,nn)=0; szp(0,nn)=0; szr(1,nn)=0; szp(1,nn)=0
 end select
 end do

    ll=-1; ra2=pi/(2*domlen); ra3=pi/szth1
 do l=0,lmx
    rr(l,:)=nsz(0,:)*szr(0,:)*max(szp(0,:)-ss(l,:),zero)+nsz(1,:)*szr(1,:)*max(ss(l,:)-szp(1,:),zero)
    ! rpt-this is sigma(x,y,z)
    de(l,1)=szco*(1-0.125*(1+cos(pi*rr(l,1)))*(1+cos(pi*rr(l,2)))*(1+cos(pi*rr(l,3))))
    ! rpt-this is lambda(x)
    de(l,2)=half*(1+cos(ra2*(min(ss(l,1),domlen)+domlen)))
    de(l,3)=sin(ra3*min(ss(l,1)+domlen,szth1))**2
 if(de(l,1)-sml>0) then
    ll=ll+1; de(ll,5)=l+sml ! rpt-this gives the l's containing sponge points
 end if
 end do
    lsz=ll ! rpt-total number of points in sponge zone
 if(lsz/=-1) then
    allocate(lcsz(0:lsz),asz(0:lsz),bsz(0:lsz),csz(0:lsz))
 do ll=0,lsz; l=de(ll,5); lcsz(ll)=l
    ! rpt-asz=sigma(x,y,z) and bsc=sigma(x,y,z)*lambda(x)
    asz(ll)=de(l,1)*de(l,2)/yaco(l); bsz(ll)=hamm1*de(l,1)*(1-de(l,2))/yaco(l)
    csz(ll)=de(l,3)
 end do
 end if

    !===GUST===
    ll=-1; ra0=tlb*(2.5_k8-cutlb)
 do lh=0,lsz; l=lcsz(lh)
 ! rpt-mark points where the gust is going to happen
 if(ss(l,1)-szp(0,1)<0.and.abs(ss(l,2))-ra0<0) then
    ll=ll+1; de(ll,5)=l+sml
 end if
 end do
    ltz=ll; ntz=litr*slit/tla
 if(ltz/=-1) then
    allocate(lctz(0:ltz),tt(0:ntz),vit(0:ltz,3),vito(0:ltz,0:ntz,3)); lp=2*nrec*(ltz+1)*(ntz+1)
 do ll=0,ltz; l=de(ll,5); lctz(ll)=l
    vit(ll,:)=ss(l,:)
 end do
    fctr=slit/(ntz*uoo(1)); tt(:)=fctr*(/0:ntz/); vito(:,:,:)=0
 if(nito==0) then
 do nn=1,nits
    ve(:)=(-9+mxc(nn))*sit(nn,:)*sit(nn,:); ra1=-36-60*mxc(nn); ra2=2*mxc(nn)/3.0_k8
 do i=0,1; do k=-1,1
    dm(:)=(/xit(nn)-i*slit,yit(nn),zit(nn)-k*span/)
 do ii=0,ntz
    rv(:)=dm(:)+tt(ii)*uoo(:)
    de(0:ltz,1)=vit(:,1)-rv(1); de(0:ltz,2)=vit(:,2)-rv(2); de(0:ltz,3)=vit(:,3)-rv(3)
    de(0:ltz,4)=de(0:ltz,1)*de(0:ltz,1)+de(0:ltz,2)*de(0:ltz,2)+de(0:ltz,3)*de(0:ltz,3)
    de(0:ltz,5)=de(0:ltz,4)*de(0:ltz,4)
    rr(0:ltz,1)=ve(1)*de(0:ltz,5); rr(0:ltz,2)=ve(2)*de(0:ltz,5); rr(0:ltz,3)=ve(3)*de(0:ltz,5)
    de(0:ltz,5)=ra1*de(0:ltz,4)
    rr(0:ltz,1)=ait(nn,1)*de(0:ltz,5)*exp(rr(0:ltz,1))*(1+ra2*rr(0:ltz,1))
    rr(0:ltz,2)=ait(nn,2)*de(0:ltz,5)*exp(rr(0:ltz,2))*(1+ra2*rr(0:ltz,2))
    rr(0:ltz,3)=ait(nn,3)*de(0:ltz,5)*exp(rr(0:ltz,3))*(1+ra2*rr(0:ltz,3))
    vito(:,ii,1)=vito(:,ii,1)+rr(0:ltz,2)*de(0:ltz,3)-rr(0:ltz,3)*de(0:ltz,2)
    vito(:,ii,2)=vito(:,ii,2)+rr(0:ltz,3)*de(0:ltz,1)-rr(0:ltz,1)*de(0:ltz,3)
    vito(:,ii,3)=vito(:,ii,3)+rr(0:ltz,1)*de(0:ltz,2)-rr(0:ltz,2)*de(0:ltz,1)
 end do
 end do; end do
 if(myid==mo(9)) then; write(*,"('Vortex',i3,' Done')") nn; end if
 end do
    open(9,file=cturb); close(9,status='delete')
    open(9,file=cturb,access='direct',form='unformatted',recl=lp)
    write(9,rec=1) vito(:,:,1); write(9,rec=2) vito(:,:,2); write(9,rec=3) vito(:,:,3)
    close(9)
    vito(:,:,:)=cfit*vito(:,:,:)
 else
    open(9,file=cturb,access='direct',form='unformatted',recl=lp)
    read(9,rec=1) vito(:,:,1); read(9,rec=2) vito(:,:,2); read(9,rec=3) vito(:,:,3)
    close(9)
    vito(:,:,:)=cfit*vito(:,:,:)
 end if
 end if

 if(myid==mo(9)) then
    fctr=one/lze0
 do l=0,lze0
    ra1=-domlen; ra2=0; ra3=(-half+l*fctr)*span
    iit(l)=minloc((vit(:,1)-ra1)**2+(vit(:,2)-ra2)**2+(vit(:,3)-ra3)**2,1)-1
 end do
    open(9,file='inflowsignal.dat',status='replace',access='direct',form='formatted',recl=16)
 do ii=0,ntz; lp=(1+3*lze0)*ii
    write(9,'(es15.7)',rec=lp+1) tt(ii)
 do l=0,lze0-2; i=iit(l)
    write(9,'(es15.7)',rec=lp+3*l+2) vito(i,ii,1)
    write(9,'(es15.7)',rec=lp+3*l+3) vito(i,ii,2)
    write(9,'(es15.7)',rec=lp+3*l+4) vito(i,ii,3)
 end do
    l=lze0-1; i=iit(l)
    write(9,'(es15.7)',rec=lp+3*l+2) vito(i,ii,1)
    write(9,'(es15.7)',rec=lp+3*l+3) vito(i,ii,2)
    write(9,'(es15.7,a)',rec=lp+3*l+4) vito(i,ii,3),achar(10)
 end do
    close(9)
 end if

 end subroutine spongeup

!===== INITIAL CONDITIONS

 subroutine initialo

    itag=1; fctr=one/lze0
 do l=0,lze0
 if(l==0) then
    ra1=0; ra2=domlen-szth2; ra3=0
 else
    ra1=-half; ra2=0; ra3=(-half+l*fctr)*span
 end if
    rr(:,1)=(ss(:,1)-ra1)**2+(ss(:,2)-ra2)**2+(ss(:,3)-ra3)**2; vmpi(myid)=minval(rr(:,1))
 if(myid==0) then
 do mp=1,mpro
    ir=mp; call MPI_IRECV(vmpi(mp),1,MPI_REAL8,mp,itag,icom,ireq(ir),ierr)
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if
 else
	call MPI_SEND(vmpi(myid),1,MPI_REAL8,0,itag,icom,ierr)
 end if
    call MPI_BCAST(vmpi(:),npro,MPI_REAL8,0,icom,ierr)
    idsgnl(l)=minloc(vmpi,1)-1; lsgnl(l)=minloc(rr(:,1),1)-1
 end do

    ra1=1/vk1; ra2=-half*ra1**2
 do m=0,1
    mm=m*domlen-half*domlen

    de(:,2)=ra1*ss(:,2)
    de(:,3)=-ra1*(ss(:,1)+mm)
    de(:,4)=vk2*exp(ra2*((ss(:,1)+mm)**2+ss(:,2)**2))
    de(:,5)=m*de(:,5)+de(:,4)

    qa(:,2)=m*qa(:,2)+de(:,2)*de(:,4)
    qa(:,3)=m*qa(:,3)+de(:,3)*de(:,4)
 end do
    p(:)=poo*(1-half*gamm1*de(:,5)**2)**(gam*hamm1)

    qa(:,1)=rhooo*(p(:)/poo)**ham
    qa(:,2)=qa(:,2)*qa(:,1)
    qa(:,3)=qa(:,3)*qa(:,1)
    qa(:,4)=0
    qa(:,5)=hamm1*p(:)+half*(qa(:,2)**2+qa(:,3)**2+qa(:,4)**2)/qa(:,1)

 end subroutine initialo

!===== SPONGE IMPLEMENTATION

 subroutine spongego

 if(ltz/=-1) then ! rpt-ltz=# of points involved in inflow gust
    vit(:,:)=0
 if(timo-tgusto+dtk>0) then
    ra0=timo-tgusto+dtk; ra1=slit/uoo(1); ra2=ra0/ra1; ra3=ra0-ra1*int(ra2)
    is=0; ie=ntz; ii=minloc(abs(tt(:)-ra3),1)-1
 do jj=-2,2
    ilag(jj)=min(max(ii+jj,is),ie); tlag(jj)=tt(ilag(jj))
 end do
 if(ii-is==0) then; ilag(-2:-1)=(/ie-2,ie-1/); tlag(-2:-1)=tt(ilag(-2:-1))-ra1; end if
 if(ii-is==1) then; ilag(-2)=ie-1; tlag(-2)=tt(ilag(-2))-ra1; end if
 if(ie-ii==1) then; ilag(2)=is+1; tlag(2)=tt(ilag(2))+ra1; end if
 if(ie-ii==0) then; ilag(1:2)=(/is+1,is+2/); tlag(1:2)=tt(ilag(1:2))+ra1; end if
    alag(:)=ra3-tlag(:); fctr=sin(pi*min(0.1_k8*ra0,half))**2
 do jj=-2,2
    blag(:)=tlag(jj)-tlag(:); ao=fctr; bo=1
 do ii=-2,2
 if(ii/=jj) then
    ao=ao*alag(ii); bo=bo*blag(ii)
 end if
 end do
    ii=ilag(jj); res=ao/bo
    vit(:,:)=vit(:,:)+res*vito(:,ii,:)
 end do
 end if
 end if

 do ll=0,lsz; l=lcsz(ll)
    rr(l,:)=0
 end do
 if (timo+dtk>tgustd) then
    ra0=timo-tgustd+dtk; ra1=min(ra0*half*pi/tguste,half*pi)
    fctr=cos(ra1)**2
 else
    fctr=1
 end if
 do ll=0,ltz; l=lctz(ll)
    rr(l,:)=csz(ll)*vit(ll,:)*fctr
 end do
    fctr=half*gamm1
 do ll=0,lsz; l=lcsz(ll)
    res=(1-fctr*(rr(l,1)**2+rr(l,2)**2+rr(l,3)**2))**hamm1
    de(l,1)=de(l,1)+asz(ll)*(qa(l,1)-res)
    de(l,2:4)=de(l,2:4)+asz(ll)*(qa(l,2:4)-qa(l,1)*rr(l,:))
    de(l,5)=de(l,5)+bsz(ll)*(p(l)-poo)
 end do

 end subroutine spongego

!===== JUNCTION AVERAGING

 subroutine junction

 integer(k4) :: njct

    is=mbk; ie=0; kk=5*(lze+1)-1
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))
    ns=mo(is)+(kp+1)*npc(is,2)*npc(is,1)-1; ne=mo(ie)+kp*npc(ie,2)*npc(ie,1)

 do nn=0,5; 
    select case(nn);
    case(0);  mp=ns; njct=nn  ;itag=njct
    case(1);  mp=ne; njct=nn  ;itag=njct
    case(2);  mp=ns; njct=nn+1;itag=njct
    case(3);  mp=ne; njct=nn+1;itag=njct
    case(4);  mp=ns; njct=nn+2;itag=njct
    case(5);  mp=ne; njct=nn+2;itag=njct
    end select
    
    do np=1,4
      select case(np)
      case(1); ip=1; jp=1; mm=njct; case(2); ip=0; jp=1;   mm=njct+1
      case(3); ip=1; jp=0; mm=njct+3; case(4); ip=0; jp=0; mm=njct+4
      end select
      mjct(np)=mo(mm)+kp*npc(mm,2)*npc(mm,1)+jp*(npc(mm,2)-1)*npc(mm,1)+ip*(npc(mm,1)-1)
      if(myid==mjct(np)) then
           i=ip*lxi; j=jp*let
        do k=0,lze; l=indx3(i,j,k,1)
        do m=1,5; jk=(m-1)*(lze+1)+k
           drva1(jk,1,0)=qa(l,m)
        end do
        end do
           call MPI_SEND(drva1(0:kk,1,0),kk+1,MPI_REAL8,mp,itag,icom,ierr)
           call MPI_RECV(drva1(0:kk,1,1),kk+1,MPI_REAL8,mp,itag,icom,ista,ierr)
        do k=0,lze; l=indx3(i,j,k,1)
        do m=1,5; jk=(m-1)*(lze+1)+k
           qa(l,m)=drva1(jk,1,1)
        end do
        end do
      end if
    end do

    if(myid==mp) then
       call MPI_RECV(drva1(0:kk,1,0),kk+1,MPI_REAL8,mjct(1),itag,icom,ista,ierr)
       call MPI_RECV(drva1(0:kk,1,1),kk+1,MPI_REAL8,mjct(2),itag,icom,ista,ierr)
       call MPI_RECV(drva1(0:kk,2,0),kk+1,MPI_REAL8,mjct(3),itag,icom,ista,ierr)
       call MPI_RECV(drva1(0:kk,2,1),kk+1,MPI_REAL8,mjct(4),itag,icom,ista,ierr)
       drva1(0:kk,3,0)=0.25_k8*(drva1(0:kk,1,0)+drva1(0:kk,1,1)+drva1(0:kk,2,0)+drva1(0:kk,2,1))
    do np=1,4
       call MPI_SEND(drva1(0:kk,3,0),kk+1,MPI_REAL8,mjct(np),itag,icom,ierr)
    end do
    end if
 end do

 end subroutine junction

!===== SIGNAL RECORDING

 subroutine signalgo

    lp=(2+3*lze0)*nsigi

    m=0; l=lsgnl(m)
 if(myid==idsgnl(m)) then
    write(1,'(es15.7)',rec=lp+1) timo
    write(1,'(es15.7)',rec=lp+2) gam*p(l)-1
 end if
 do m=1,lze0-1; l=lsgnl(m)
 if(myid==idsgnl(m)) then; ve(:)=qa(l,2:4)/qa(l,1)
    write(1,'(es15.7)',rec=lp+3*m) ve(1)
    write(1,'(es15.7)',rec=lp+3*m+1) ve(2)
    write(1,'(es15.7)',rec=lp+3*m+2) ve(3)
 end if
 end do
    m=lze0; l=lsgnl(m)
 if(myid==idsgnl(m)) then; ve(:)=qa(l,2:4)/qa(l,1)
    write(1,'(es15.7)',rec=lp+3*m) ve(1)
    write(1,'(es15.7)',rec=lp+3*m+1) ve(2)
    write(1,'(es15.7,a)',rec=lp+3*m+2) ve(3),achar(10)
 end if

 end subroutine signalgo

!===== FINAL RESULTS

 subroutine finalout

 integer(k4) :: totVar
 real(k8),dimension(:),allocatable :: delt


   totVar=5
 
 if(ltz/=-1) then; deallocate(lctz,tt,vit,vito); end if

    ns=0; ne=ndata; allocate(delt(ns:ne))
    fctr=half/(times(ne)-times(ns))
    delt(ns)=fctr*(times(ns+1)-times(ns)); delt(ne)=fctr*(times(ne)-times(ne-1))
 do n=ns+1,ne-1
    delt(n)=fctr*(times(n+1)-times(n-1))
 end do
 do m = 1, totVar
    rr(:,1)=0
 do n=0,ndata
    read(0,rec=(n*totVar)+ngrec+m) varr(:)
    rr(:,1)=rr(:,1)+delt(n)*varr(:)
 end do
 !   ss(:,1)=0
 !do n=0,ndata
 !   varr(:)=varr(:)-rr(:,1)
 !   write(0,rec=(n*5)+ngrec+m) varr(:)
 !   ss(:,1)=ss(:,1)+delt(n)*varr(:)**2
 !end do
    !varr=sqrt(ss(:,1))
    varr(:)=rr(:,1)
    nwrec=nwrec+1
    write(0,rec=nwrec) varr(:)
 end do

 end subroutine finalout

!=====

 end module problemcase

!*****
