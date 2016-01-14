!*****
!***** 3D SINGLE AEROFOIL-TURBULENCE INTERACTION
!*****

 module problemcase

 use mpi
 use subroutineso
 use gridgen
 implicit none

 integer(kind=ni) :: nbody,nthick,ngridv,nito,nits,litr,lsz,ltz,ntz
 integer(kind=ni),dimension(:),allocatable :: iit,idsgnl,lsgnl
 integer(kind=ni),dimension(4) :: mjct
 real(kind=nr),dimension(:,:,:),allocatable :: vito
 real(kind=nr),dimension(:,:),allocatable :: vit
 real(kind=nr),dimension(:),allocatable :: tt
 real(kind=nr),dimension(0:1,3) :: szr,szp
 real(kind=nr),dimension(3,0:1) :: tam
 real(kind=nr),dimension(0:1) :: tl0,tl1,tlw
 real(kind=nr) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,szco
 real(kind=nr) :: tgusto,slit,gaus,cfit,tla,tlb,cutlb
 real(kind=nr) :: denxit
 ! rpt-tandem variables
 real(k8) :: gap,c1,c2,delt1,delt2
 ! rpt-grid space modifiers for grid generation
 real(k8) :: ximod,etamod
 ! rpt- Declare lze0 which is = lzebk(0)
 integer(kind=ni) :: lze0

 contains

!===== SUBINPUT PARAMETERS

 subroutine inputext

 integer(k4), dimension(:) :: npcx(0:bkx-1),npcy(0:bky-1),npcz(0:bkz-1)

    open(9,file='inputp.dat')
    read(9,*) cinput,lxibk(0:bkx-1) ! # points in xi per columns
    read(9,*) cinput,letbk(0:bky-1) ! # points in eta per row
    read(9,*) cinput,lzebk(0:bkz-1) ! # points in zeta per plane
    read(9,*) cinput,nbody,nthick
    read(9,*) cinput,ngridv
    read(9,*) cinput,npcx(0:bkx-1)  ! # processors in xi per columns
    read(9,*) cinput,npcy(0:bky-1)  ! # processors in eta per row
    read(9,*) cinput,npcz(0:bkz-1)  ! # processors in zeta per plane
    read(9,*) cinput,nito,nits,litr
    read(9,*) cinput,smgrid,domlen
    read(9,*) cinput,ximod,etamod   ! xi and eta modifiers for LE and TE spacing
    read(9,*) cinput,span,wlew,wlea
    read(9,*) cinput,gap,c1,c2
    read(9,*) cinput,delt1,delt2
    read(9,*) cinput,szth1,szth2,szxt,szco
    read(9,*) cinput,tgusto
    read(9,*) cinput,slit,gaus,cfit
    read(9,*) cinput,tl0(0),tl1(0),tlw(0)
    read(9,*) cinput,tl0(1),tl1(1),tlw(1)
    read(9,*) cinput,tam(:,0)
    read(9,*) cinput,tam(:,1)
    read(9,*) cinput,cutlb
    close(9)

    tla=minval((/tl0(:),tl1(:)/)); tlb=maxval((/tl0(:),tl1(:)/))
    ! rpt-from degrees to radians
    delt1=delt1*pi/180.0e0; delt2=delt2*pi/180.0e0

    ! Use input data to fill old arrays
    do k = 0, bkz-1
       do j = 0, bky-1
          do i = 0, bkx-1; l=k*(bkx*bky)+j*bkx+i
             lximb(l)=lxibk(i)
             letmb(l)=letbk(j)
             lzemb(l)=lzebk(k)
             
             npc(l,1)=npcx(i)
             npc(l,2)=npcy(j)
             npc(l,3)=npcz(k)
          end do
       end do
    end do
    lze0 = lzebk(0)

    lze0=lzemb(0)

    allocate(mxc(nits),ran(nits,3),sit(nits,3),ait(nits,3),xit(nits),yit(nits),zit(nits))
    allocate(iit(0:lze0),idsgnl(0:lze0),lsgnl(0:lze0))

    open(9,file='randnum.dat')
 do m=1,3
    read(9,*) ran(:,m); read(9,*) ait(:,m)
 end do
!    read(9,*) xit(:); read(9,*) yit(:); read(9,*) zit(:)
    read(9,*) yit(:); read(9,*) zit(:)
    close(9)
 do nn=1,nits
	if(ran(nn,1)<=gaus) then; mxc(nn)=0; else; mxc(nn)=1; end if
 end do
    res=(two*(two-cutlb)*tlb*span*slit/nits)**one_three
 do m=1,3; do nn=1,nits
    ran(nn,m)=tl0(mxc(nn))+(tl1(mxc(nn))-tl0(mxc(nn)))*ran(nn,m)**tlw(mxc(nn))
    sit(nn,m)=one/ran(nn,m)**two
    ait(nn,m)=res*tam(m,mxc(nn))*(two*ait(nn,m)-one)*sit(nn,m)**1.5_nr
 end do; end do
    denxit=slit/(two*sum(ran(:,1)))
    xit(1)=-domlen+half*(szth1+slit)-denxit*ran(1,1)
 do nn=2,nits
    xit(nn)=xit(nn-1)-denxit*(ran(nn-1,1)+ran(nn,1))
 end do
 do nn=1,nits
    res=three*span-two*max(ran(nn,1),ran(nn,2),ran(nn,3))
!    xit(nn)=-domlen+half*(szth1+slit)-slit*xit(nn)
    yit(nn)=(one-cutlb)*ran(nn,2)*(two*yit(nn)-one)
    zit(nn)=min(span,res)*(zit(nn)-half)
 end do

 end subroutine inputext

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

 subroutine domdcomp

 integer(k4) :: ll,bcw,bcinout,bcperiod,bcinter
 bcw=20+5*nviscous
 bcinout=10
 bcperiod=45
 bcinter=30

 ! rpt-setting the neighbouring blocks and BCS
 
 ll=mod(mb,bkx)
 if (ll==0) then ! rpt-Left boundary blocks
    ms(1)=mb+(bkx-1); me(1)=mb+1
    nbcs(1)=bcinout; nbce(1)=bcinter
 elseif (ll==(bkx-1)) then ! rpt-Right boundary blocks
    ms(1)=mb-1; me(1)=mb-(bkx-1)
    nbcs(1)=bcinter; nbce(1)=bcinout
 else ! rpt-Middle blocks
    ms(1)=mb-1; me(1)=mb+1
    nbcs(1)=bcinter; nbce(1)=bcinter
 end if
 if (mb<bkx) then ! rpt-Bottom boundary blocks
    ms(2)=mb+bkx*(bky-1); me(2)=mb+bkx
    nbcs(2)=bcinout; nbce(2)=bcinter
 elseif(mb>(mbk-bkx)) then ! rpt-Top boundary blocks
    ms(2)=mb-bkx; me(2)=mb-(bkx*(bky-1))
    nbcs(2)=bcinter; nbce(2)=bcinout
 else ! rpt-Middle blocks
    ms(2)=mb-bkx; me(2)=mb+bkx
    nbcs(2)=bcinter; nbce(2)=bcinter
 end if

    ! rpt-Periodic in span
    nbcs(3)=bcperiod; nbce(3)=bcperiod
    ms(3)=mb; me(3)=mb

    ! rpt-Wall Boundary condition
    if(mb==1) nbce(2)=bcw
    if(mb==4) nbcs(2)=bcw
    ! rpt-Interface condition
    if(mb==2) nbce(2)=35
    if(mb==5) nbcs(2)=35

 end subroutine domdcomp

!===== GRID GENERATION

 subroutine makegrid

    ! rpt-new grid generation file, see Grid.90
    call gridaerofoil(ngridv,nthick,smgrid,&
         domlen,span,wlew,wlea,szth1,szth2,szxt,&
         c1,delt1,ximod,etamod)

 end subroutine makegrid

!===== INITIAL CONDITIONS

 subroutine initialo

    itag=1; fctr=one/lze0
 do l=0,lze0
 if(l==0) then
    ra1=zero; ra2=domlen-szth2; ra3=zero
 else
    ra1=-half; ra2=zero; ra3=(-half+l*fctr)*span
 end if
    rr(:,1)=(ss(:,1)-ra1)**two+(ss(:,2)-ra2)**two+(ss(:,3)-ra3)**two; vmpi(myid)=minval(rr(:,1))
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

    qa(:,1)=rhooo; qa(:,2:4)=zero; qa(:,5)=hamm1*poo

 end subroutine initialo

!===== SETTING UP SPONGE ZONE PARAMETERS

 subroutine spongeup

 do nn=1,3
 ! rpt-set nsz=1 if contains boundary region
 if(nbcs(nn)==10) then; nsz(0,nn)=1; else; nsz(0,nn)=0; end if
 if(nbce(nn)==10) then; nsz(1,nn)=1; else; nsz(1,nn)=0; end if
 select case(nn)
 case(1); szr(0,nn)=1/szth1;        szp(0,nn)=-0.7e0*domlen+szth1;  
          szr(1,nn)=1/(szth1+szxt); szp(1,nn)=domlen-szth1
 case(2); szr(0,nn)=1/szth2;        szp(0,nn)=-0.7e0*domlen+szth2; 
          szr(1,nn)=1/szth2;        szp(1,nn)=domlen+szxt-szth2
 case(3); szr(0,nn)=0; szp(0,nn)=0; szr(1,nn)=0; szp(1,nn)=0
 end select
 end do

    ll=-1; ra1=0.25_nr; ra2=pi/(two*domlen)
 do l=0,lmx
    rr(l,:)=nsz(0,:)*szr(0,:)*max(szp(0,:)-ss(l,:),zero)&
           +nsz(1,:)*szr(1,:)*max(ss(l,:)-szp(1,:),zero)
    ! rpt-this is sigma(x,y,z)
    de(l,1)=szco*(one-0.125_nr*(one+cos(pi*rr(l,1)))*(one+cos(pi*rr(l,2)))*(one+cos(pi*rr(l,3))))
    ! rpt-this is lambda(x)
    de(l,2)=ra1+half*(one-ra1)*(one+cos(ra2*(min(ss(l,1),domlen)+domlen)))
 
 if(de(l,1)>zero) then
    ll=ll+1; de(ll,5)=l+sml ! rpt-this gives the l's containing sponge points
 end if
 end do
    ! rpt- Output sponge 
    if ((ngridv==1).and.(output==1)) then
       if(.not.allocated(fout)) allocate(fout(0:lmx,2))
       fout(:,1:2)=de(:,1:2)
    end if
    lsz=ll ! rpt-total number of points in sponge zone
 if(lsz/=-1) then
    allocate(lcsz(0:lsz),asz(0:lsz),bsz(0:lsz))
 do ll=0,lsz; l=de(ll,5); lcsz(ll)=l; res=one/yaco(l)
    ! rpt-asz=sigma(x,y,z) and bsc=sigma(x,y,z)*lambda(x)
    asz(ll)=res*de(l,1)*de(l,2); bsz(ll)=res*hamm1*de(l,1)
    end do
 end if

    !===GUST===
    ll=-1; ra0=tlb*(2.5_nr-cutlb)
 do lh=0,lsz; l=lcsz(lh)
 ! rpt-mark points where the gust is going to happen
 if(ss(l,1)<szp(0,1).and.abs(ss(l,2))<ra0) then
    ll=ll+1; de(ll,5)=l+sml
 end if
 end do
    ltz=ll; ntz=litr*slit/tla; ra0=pi/szth1
    !ltz=-1
 if(ltz/=-1) then
    allocate(lctz(0:ltz),atz(0:ltz),tt(0:ntz),vit(0:ltz,3),vito(0:ltz,0:ntz,3)); lp=nrecd*(ltz+1)*(ntz+1)
    do ll=0,ltz; l=de(ll,5); lctz(ll)=l
       vit(ll,:)=ss(l,:); atz(ll)=sin(min(ra0*(ss(l,1)+domlen),halfpi))**two
    end do
       fctr=slit/(ntz*uoo(1)); do ii=0,ntz; tt(ii)=ii*fctr; end do; vito(:,:,:)=zero
    if(nito==0) then
       do nn=1,nits
          ve(:)=(-9+mxc(nn))*sit(nn,:)*sit(nn,:); ra1=-36-60*mxc(nn); ra2=two_three*mxc(nn)
       do i=0,1; do k=-1,1
          dm(:)=(/xit(nn)-i*slit,yit(nn),zit(nn)-k*span/)
       do ii=0,ntz
          rv(:)=dm(:)+tt(ii)*uoo(:)
          de(0:ltz,1)=vit(:,1)-rv(1); de(0:ltz,2)=vit(:,2)-rv(2); de(0:ltz,3)=vit(:,3)-rv(3)
          de(0:ltz,4)=de(0:ltz,1)*de(0:ltz,1)+de(0:ltz,2)*de(0:ltz,2)+de(0:ltz,3)*de(0:ltz,3)
          de(0:ltz,5)=de(0:ltz,4)*de(0:ltz,4)
          rr(0:ltz,1)=ve(1)*de(0:ltz,5);
          rr(0:ltz,2)=ve(2)*de(0:ltz,5);
          rr(0:ltz,3)=ve(3)*de(0:ltz,5)
          de(0:ltz,5)=ra1*de(0:ltz,4)
          rr(0:ltz,1)=ait(nn,1)*de(0:ltz,5)*exp(rr(0:ltz,1))*(1+ra2*rr(0:ltz,1))
          rr(0:ltz,2)=ait(nn,2)*de(0:ltz,5)*exp(rr(0:ltz,2))*(1+ra2*rr(0:ltz,2))
          rr(0:ltz,3)=ait(nn,3)*de(0:ltz,5)*exp(rr(0:ltz,3))*(1+ra2*rr(0:ltz,3))
          vito(:,ii,1)=vito(:,ii,1)+rr(0:ltz,2)*de(0:ltz,3)-rr(0:ltz,3)*de(0:ltz,2)
          vito(:,ii,2)=vito(:,ii,2)+rr(0:ltz,3)*de(0:ltz,1)-rr(0:ltz,1)*de(0:ltz,3)
          vito(:,ii,3)=vito(:,ii,3)+rr(0:ltz,1)*de(0:ltz,2)-rr(0:ltz,2)*de(0:ltz,1)
       end do
       end do; end do
       if(myid==mo(mbk+1-bkx)) then; write(*,"('Vortex',i3,' Done')") nn; end if
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

 !if((myid==mo(mbk+1-bkx)).and.(ltz/=-1)) then
 if(myid==mo(3)) then
    fctr=one/lze0
 do l=0,lze0
    ra1=-domlen; ra2=zero; ra3=(-half+l*fctr)*span
    iit(l)=minloc((vit(:,1)-ra1)**two+(vit(:,2)-ra2)**two+(vit(:,3)-ra3)**two,1)-1
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

!===== SPONGE IMPLEMENTATION

 subroutine spongego

 if(ltz/=-1) then ! rpt-ltz=# of points involved in inflow gust
    vit(:,:)=zero
 if(timo>tgusto) then
    ra0=timo-tgusto+dtk; ra1=slit/uoo(1); ra2=ra0/ra1; ra3=ra0-ra1*int(ra2,kind=ni)
    is=0; ie=ntz; ii=minloc(abs(tt(:)-ra3),1)-1
 do jj=-2,2
    ilag(jj)=min(max(ii+jj,is),ie); tlag(jj)=tt(ilag(jj))
 end do
 if(ii-is==0) then; ilag(-2:-1)=(/ie-2,ie-1/); tlag(-2:-1)=tt(ilag(-2:-1))-ra1; end if
 if(ii-is==1) then; ilag(-2)=ie-1; tlag(-2)=tt(ilag(-2))-ra1; end if
 if(ie-ii==1) then; ilag(2)=is+1; tlag(2)=tt(ilag(2))+ra1; end if
 if(ie-ii==0) then; ilag(1:2)=(/is+1,is+2/); tlag(1:2)=tt(ilag(1:2))+ra1; end if
    alag(:)=ra3-tlag(:)
 do jj=-2,2
    blag(:)=tlag(jj)-tlag(:); ao=one; bo=one
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
    rr(l,1)=one; ss(l,:)=zero
 end do
    fctr=half*gamm1
 do ll=0,ltz; l=lctz(ll)
    rr(l,1)=(one-fctr*(vit(ll,1)*vit(ll,1)+vit(ll,2)*vit(ll,2)+vit(ll,3)*vit(ll,3)))**hamm1
    ss(l,:)=atz(ll)*rr(l,1)*vit(ll,:)
 end do
 do ll=0,lsz; l=lcsz(ll)
    de(l,1)=de(l,1)+asz(ll)*(qa(l,1)-rr(l,1))
    de(l,2:4)=de(l,2:4)+asz(ll)*(qa(l,2:4)-ss(l,:))
    de(l,5)=de(l,5)+bsz(ll)*(p(l)-ham*rr(l,1)**gam)
 end do

 end subroutine spongego

!===== JUNCTION AVERAGING

 subroutine junction

 integer(k4) :: njct

    is=mbk; ie=0; kk=5*(lze+1)-1
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))
    ns=mo(is)+(kp+1)*npc(is,2)*npc(is,1)-1; ne=mo(ie)+kp*npc(ie,2)*npc(ie,1)

 do nn=0,(bkx-1)*(bky-1)-1; 
    if (mod(nn,2)==0) then
       mp=ns
    else 
       mp=ne
    end if
    njct=nn+nn/(bkx-1); itag=njct
    do np=1,4
      select case(np)
      case(1); ip=1; jp=1; mm=njct; case(2); ip=0; jp=1;   mm=njct+1
      case(3); ip=1; jp=0; mm=njct+bkx; case(4); ip=0; jp=0; mm=njct+bkx+1
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
       !drva1(0:kk,3,0)=0.25_k8*(drva1(0:kk,1,0)+drva1(0:kk,1,1)+drva1(0:kk,2,0)+drva1(0:kk,2,1))
       !?????????
       ! rpt- this won't work for other block structures
       do jk=0,lze ! rpt - Density
          drva1(jk,3,0)=0.25_nr*(drva1(jk,1,0)+drva1(jk,1,1)+drva1(jk,2,0)+drva1(jk,2,1))
       end do
       do m=1,3; ll=m*(lze+1) ! rpt- Velocities
       do jk=ll,ll+lze
!          drva1(jk,3,0)=-umf(m)*drva1(jk-ll,3,0)
          drva1(jk,3,0)=half*(drva1(jk,1,1-nn)+drva1(jk,2,1-nn))
       end do
       end do
       do jk=4*(lze+1),kk ! rpt- Energy
          drva1(jk,3,0)=0.25_nr*(drva1(jk,1,0)+drva1(jk,1,1)+drva1(jk,2,0)+drva1(jk,2,1))
       end do
       !?????????
       do np=1,4
          call MPI_SEND(drva1(0:kk,3,0),kk+1,MPI_REAL8,mjct(np),itag,icom,ierr)
       end do
    end if
 end do

 end subroutine junction

!===== WALL SLOT BOUNDARY IDENTIFICATION

 subroutine extrabcc(flag)

 real(kind=nr),intent(inout) :: flag
 real(kind=nr) :: slitl,slitw

    flag=zero
 if(ss(l,1)<ra0) then
    slitl=wlea; slitw=0.4_nr*slitl
    ra0=-half+wlea+slitl; ra1=-half*span+0.25_nr*wlew-half*slitw; ra2=ra1+slitw; ie=(span+sml)/wlew-1
 do ii=0,ie; ra3=ii*wlew
    res=(ss(l,3)-(ra1+ra3))*(ss(l,3)-(ra2+ra3))
 if(res<zero.and.nextrabc==1) then
    flag=one
 end if
 end do
 end if

 end subroutine extrabcc

!===== WALL SLIT BOUNDARY IMPLEMENTATION

 subroutine extrabcs

    if(ra0*(vn+vs)>zero) then; cha(1:3)=dha(1:3); end if
    if(ra0*(vn+vs+ao)>zero) then; cha(4)=dha(4); end if
    if(ra0*(vn+vs-ao)>zero) then; cha(5)=dha(5); end if

 end subroutine extrabcs

!===== SIGNAL RECORDING

 subroutine signalgo

    lp=(2+3*lze0)*nsigi

    m=0; l=lsgnl(m)
 if(myid==idsgnl(m)) then
    write(1,'(es15.7)',rec=lp+1) timo
    write(1,'(es15.7)',rec=lp+2) gam*p(l)-one
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

!=====

 end module problemcase

!*****
