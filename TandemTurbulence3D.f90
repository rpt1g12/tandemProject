!*****
!***** 3D TANDEM AEROFOIL-TURBULENCE INTERACTION
!*****

 module problemcase

 use mpi
 use subroutineso
 use gridgen
 implicit none

 integer :: ngridv,nits,litr,lsz,ltz,ntz
 integer,dimension(4) :: mjct
 real(nr),dimension(:,:,:),allocatable :: vito
 real(nr),dimension(:,:),allocatable :: vit
 real(nr),dimension(:),allocatable :: tt
 real(nr),dimension(0:1,3) :: szr,szp
 real(nr),dimension(3,0:1) :: tam
 real(nr),dimension(0:1) :: tl0,tl1,tlw
 real(nr) :: smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,szco
 real(nr) :: gap,c1,c2,delt1,delt2
 real(nr) :: ximod,etamod
 real(nr) :: tgusto,eps,ck1,ck2,ck3,amp1,amp2,amp3,vk1,vk2,slit,tla,tlb
 real(nr) :: denxit

 contains

!===== SUBINPUT PARAMETERS

 subroutine inputext

    open(1,file='inputp.dat',shared)
    read(1,*) cinput,lxi0,lxi1,lxi2,lxi3,lxi4
    read(1,*) cinput,let0
    read(1,*) cinput,lze0
    read(1,*) cinput,ngridv
    read(1,*) cinput,npc(:,1)
    read(1,*) cinput,npc(:,2)
    read(1,*) cinput,npc(:,3)
    read(1,*) cinput,nits,litr
    read(1,*) cinput,smgrid,domlen
    read(1,*) cinput,ximod,etamod
    read(1,*) cinput,span,wlew,wlea
    read(1,*) cinput,gap,c1,c2
    read(1,*) cinput,delt1,delt2
    read(1,*) cinput,szth1,szth2,szxt,szco
    read(1,*) cinput,tgusto
    read(1,*) cinput,ck1,ck2,ck3
    read(1,*) cinput,amp2,amp3
    read(1,*) cinput,vk1,vk2
    read(1,*) cinput,slit
    read(1,*) cinput,tl0(0),tl1(0),tlw(0)
    read(1,*) cinput,tl0(1),tl1(1),tlw(1)
    read(1,*) cinput,tam(:,0)
    read(1,*) cinput,tam(:,1)
    close(1)

    ck1=2*pi*ck1; ck2=2*pi*ck2; ck3=2*pi*ck3
    amp1=-(amp2*ck2+amp3*ck3)/ck1
    eps=maxval((/abs(amp2),sum(tam)/6,sml/))
    tla=minval((/tl0(:),tl1(:)/)); tlb=maxval((/tl0(:),tl1(:)/))
    delt1=delt1*pi/180.0e0; delt2=delt2*pi/180.0e0

    tsam=max(tmax-2*pi/max(ck1*uoo(1)+ck2*uoo(2)+ck3*uoo(3),sml),tsam)

    lximb(:)=(/lxi0,lxi1,lxi2,lxi3,lxi4,lxi0,lxi1,lxi2,lxi3,lxi4/)
    letmb(:)=(/let0,let0,let0,let0,let0,let0,let0,let0,let0,let0/)
    lzemb(:)=(/lze0,lze0,lze0,lze0,lze0,lze0,lze0,lze0,lze0,lze0/)

    allocate(mxc(nits),ran(nits),xit(nits),yit(nits),zit(nits),sit(nits),ait(nits),bit(nits),cit(nits))
    call random_number(ran); call random_number(ait); call random_number(bit); call random_number(cit)
    call random_number(yit); call random_number(zit)
!    open(1,file='randnum.dat',shared)
!    read(1,*) ran(:); read(1,*) ait(:); read(1,*) bit(:); read(1,*) cit(:)
!    read(1,*) yit(:); read(1,*) zit(:)
!    close(1)
    mxc(:)=mod(int(1.0e8*ran(:)),2)
    ran(:)=tl0(mxc(:))+(tl1(mxc(:))-tl0(mxc(:)))*ran(:)**tlw(mxc(:))
    sit(:)=1/ran(:)**2
    denxit=slit/(2*sum(ran(:))); res=sqrt(2*tlb*span/nits)
 do nn=1,nits
    fctr=res*sit(nn)**1.5_nr
    ait(nn)=fctr*tam(1,mxc(nn))*(1-2*mod(int(1.0e8*ait(nn)),2))
    bit(nn)=fctr*tam(2,mxc(nn))*(1-2*mod(int(1.0e8*bit(nn)),2))
    cit(nn)=fctr*tam(3,mxc(nn))*(1-2*mod(int(1.0e8*cit(nn)),2))
 end do
    zit(:)=half*span*zit(:)*(1-2*mod(int(1.0e8*zit(:)),2))
    yit(:)=half*ran(:)*yit(:)*(1-2*mod(int(1.0e8*yit(:)),2))
    xit(1)=-domlen+half*(szth1+slit)-denxit*ran(1)
 do nn=2,nits
    xit(nn)=xit(nn-1)-denxit*(ran(nn-1)+ran(nn))
 end do
    call MPI_BCAST(mxc,nits,MPI_INTEGER,0,icom,ierr)
    call MPI_BCAST(sit,nits,MPI_REAL8,0,icom,ierr)
    call MPI_BCAST(ait,nits,MPI_REAL8,0,icom,ierr)
    call MPI_BCAST(bit,nits,MPI_REAL8,0,icom,ierr)
    call MPI_BCAST(cit,nits,MPI_REAL8,0,icom,ierr)
    call MPI_BCAST(zit,nits,MPI_REAL8,0,icom,ierr)
    call MPI_BCAST(yit,nits,MPI_REAL8,0,icom,ierr)
    call MPI_BCAST(xit,nits,MPI_REAL8,0,icom,ierr)

 end subroutine inputext
 

!===== DOMAIN DECOMPOSITION & BOUNDARY INFORMATION

 subroutine domdcomp

    ip=30; jp=20+5*nviscous
 select case(mb)
 case(0); nbcs(1)=10; nbce(1)=ip; nbcs(2)=10; nbce(2)=ip; nbcs(3)=45; nbce(3)=45
 case(1); nbcs(1)=ip; nbce(1)=ip; nbcs(2)=10; nbce(2)=jp; nbcs(3)=45; nbce(3)=45
 case(2); nbcs(1)=ip; nbce(1)=ip; nbcs(2)=10; nbce(2)=ip; nbcs(3)=45; nbce(3)=45
 case(3); nbcs(1)=ip; nbce(1)=ip; nbcs(2)=10; nbce(2)=jp; nbcs(3)=45; nbce(3)=45
 case(4); nbcs(1)=ip; nbce(1)=10; nbcs(2)=10; nbce(2)=35; nbcs(3)=45; nbce(3)=45
 case(5); nbcs(1)=10; nbce(1)=ip; nbcs(2)=ip; nbce(2)=10; nbcs(3)=45; nbce(3)=45
 case(6); nbcs(1)=ip; nbce(1)=ip; nbcs(2)=jp; nbce(2)=10; nbcs(3)=45; nbce(3)=45
 case(7); nbcs(1)=ip; nbce(1)=ip; nbcs(2)=ip; nbce(2)=10; nbcs(3)=45; nbce(3)=45
 case(8); nbcs(1)=ip; nbce(1)=ip; nbcs(2)=jp; nbce(2)=10; nbcs(3)=45; nbce(3)=45
 case(9); nbcs(1)=ip; nbce(1)=10; nbcs(2)=35; nbce(2)=10; nbcs(3)=45; nbce(3)=45
 end select

 select case(mb)
 case(0,5); ms(1)=mb+4; me(1)=mb+1; case(1,2,3,6,7,8); ms(1)=mb-1; me(1)=mb+1; case(4,9); ms(1)=mb-1; me(1)=mb-4
 end select
 select case(mb)
 case(0,1,2,3,4); ms(2)=mb+5; me(2)=mb+5; case(5,6,7,8,9); ms(2)=mb-5; me(2)=mb-5
 end select
    ms(3)=mb; me(3)=mb

 end subroutine domdcomp

!===== GRID GENERATION

 subroutine makegrid

    call gridaerofoil(ngridv,litr,smgrid,domlen,span,wlew,wlea,szth1,szth2,szxt,tla,tlb,gap,c1,c2,delt1,delt2,ximod,etamod)

 end subroutine makegrid

!===== SETTING UP SPONGE ZONE PARAMETERS

 subroutine spongeup

 do nn=1,3
 if(nbcs(nn)==10) then; nsz(0,nn)=1; else; nsz(0,nn)=0; end if
 if(nbce(nn)==10) then; nsz(1,nn)=1; else; nsz(1,nn)=0; end if
 select case(nn)
 case(1); szr(0,nn)=1/szth1; szp(0,nn)=-domlen+szth1; szr(1,nn)=1/(szth1+szxt)
 szp(1,nn)=half*(c1+c2)+gap+domlen-szth1
 case(2); szr(0,nn)=1/szth2; szp(0,nn)=-domlen+szth2; szr(1,nn)=1/szth2; szp(1,nn)=domlen-szth2
 case(3); szr(0,nn)=0; szp(0,nn)=0; szr(1,nn)=0; szp(1,nn)=0
 end select
 end do

    ll=-1
 do l=0,lmx
    rr(l,:)=nsz(0,:)*szr(0,:)*max(szp(0,:)-ss(l,:),0.0_nr)+nsz(1,:)*szr(1,:)*max(ss(l,:)-szp(1,:),0.0_nr)
    de(l,1)=half*szco*(1+cos(pi*(1-rr(l,1))*(1-rr(l,2))*(1-rr(l,3))))
    de(l,2)=0.25_nr*(1-tanh(ss(l,1)))+half
 if(de(l,1)-sml>0) then
    ll=ll+1; de(ll,5)=l+sml
 end if
 end do
    lsz=ll
 if(lsz/=-1) then
    allocate(lcsz(0:lsz),asz(0:lsz),bsz(0:lsz),csz(0:lsz),dsz(0:lsz))
 do ll=0,lsz; l=de(ll,5); lcsz(ll)=l
    asz(ll)=de(l,1)*de(l,2)/yaco(l); bsz(ll)=asz(ll)*de(l,2)
    csz(ll)=ss(l,1)+domlen; dsz(ll)=ck1*csz(ll)+ck2*ss(l,2)+ck3*ss(l,3)
 end do
 end if

    ll=-1; ra0=2*tlb
 do lh=0,lsz; l=lcsz(lh)
 if(ss(l,1)-szp(0,1)<0.and.abs(ss(l,2))-ra0<0) then
    ll=ll+1; de(ll,5)=l+sml
 end if
 end do
    ltz=ll; ntz=litr*slit/tla
 if(ltz/=-1) then
    allocate(lctz(0:ltz),tt(0:ntz),vit(0:ltz,3),vito(0:ltz,0:ntz,3))
 do ll=0,ltz; l=de(ll,5); lctz(ll)=l
    vit(ll,:)=ss(l,:)
 end do
    fctr=slit/(ntz*uoo(1)); tt(:)=fctr*(/0:ntz/); vito(:,:,:)=0
 do nn=1,nits
    ra1=(-9+mxc(nn))*sit(nn)*sit(nn); ra2=-36-60*mxc(nn); ra3=2*mxc(nn)/3.0_nr
 do i=0,1; do k=-1,1
    dm(:)=(/xit(nn)-i*slit,yit(nn),zit(nn)-k*span/)
 do ii=0,ntz
    rv(:)=dm(:)+tt(ii)*uoo(:)
    de(0:ltz,1)=vit(:,1)-rv(1); de(0:ltz,2)=vit(:,2)-rv(2); de(0:ltz,3)=vit(:,3)-rv(3)
    de(0:ltz,4)=de(0:ltz,1)*de(0:ltz,1)+de(0:ltz,2)*de(0:ltz,2)+de(0:ltz,3)*de(0:ltz,3)
    de(0:ltz,5)=ra1*de(0:ltz,4)*de(0:ltz,4)
    rr(0:ltz,1)=ra2*de(0:ltz,4)*exp(de(0:ltz,5))*(1+ra3*de(0:ltz,5))
    vito(:,ii,1)=vito(:,ii,1)+rr(0:ltz,1)*(cit(nn)*de(0:ltz,2)-bit(nn)*de(0:ltz,3))
    vito(:,ii,2)=vito(:,ii,2)+rr(0:ltz,1)*(ait(nn)*de(0:ltz,3)-cit(nn)*de(0:ltz,1))
    vito(:,ii,3)=vito(:,ii,3)+rr(0:ltz,1)*(bit(nn)*de(0:ltz,1)-ait(nn)*de(0:ltz,2))
 end do
 end do; end do
 if(myid==mo(5)) then; write(*,"('Vortex',i3,' Done')") nn; end if
 end do
 end if

 if(myid==mo(5)) then
    open(1,file='signal.dat')
 do ii=0,ntz
    write(1,"(4e16.8)") tt(ii),vito(0,ii,1),vito(0,ii,2),vito(0,ii,3)
 end do
    close(1)
 end if

 end subroutine spongeup

!===== INITIAL CONDITIONS

 subroutine initialo

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

 if(ltz/=-1) then
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
    alag(:)=ra3-tlag(:); fctr=sin(pi*min(ra2,half))**2
 do jj=-2,2
    blag(:)=tlag(jj)-tlag(:); ao=fctr; bo=1
 do ii=-2,2; if(ii/=jj) then
    ao=ao*alag(ii); bo=bo*blag(ii)
 end if; end do
    res=ao/bo
    vit(:,:)=vit(:,:)+res*vito(:,ilag(jj),:)
 end do
 end if
 end if

 do ll=0,lsz; l=lcsz(ll)
    rr(l,:)=0
 end do
 do ll=0,ltz; l=lctz(ll)
    rr(l,:)=qa(l,1)*vit(ll,:)
 end do
 do ll=0,lsz; l=lcsz(ll)
    de(l,1)=de(l,1)+asz(ll)*(qa(l,1)-rhooo)
    de(l,2:4)=de(l,2:4)+bsz(ll)*(qa(l,2:4)-rr(l,:))
    de(l,5)=de(l,5)+asz(ll)*(p(l)-poo)
 end do

 end subroutine spongego

!===== JUNCTION AVERAGING

 subroutine junction

    is=mbk; ie=0; kk=5*(lze+1)-1
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))
    ns=mo(is)+(kp+1)*npc(is,2)*npc(is,1)-1; ne=mo(ie)+kp*npc(ie,2)*npc(ie,1)
 do nn=0,3; itag=nn
 select case(nn); case(0); mp=ns; case(1); mp=ne; end select
 
 do np=1,4
 select case(np)
 case(1); ip=1; jp=1; mm=nn; case(2); ip=0; jp=1; mm=nn+1
 case(3); ip=1; jp=0; mm=nn+5; case(4); ip=0; jp=0; mm=nn+6
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
! if(nn==0) then
! do jk=0,lze
!    drva1(jk,3,0)=0.25_nr*(drva1(jk,1,0)+drva1(jk,1,1)+drva1(jk,2,0)+drva1(jk,2,1))
! end do
! do m=1,3; ll=m*(lze+1)
! do jk=ll,ll+lze
!    drva1(jk,3,0)=-umf(m)*drva(jk-ll,3,0)
! end do
! end do
! do jk=4*(lze+1),kk
!    drva1(jk,3,0)=0.25_nr*(drva1(jk,1,0)+drva1(jk,1,1)+drva1(jk,2,0)+drva1(jk,2,1))
! end do
! else
    drva1(0:kk,3,0)=0.25_nr*(drva1(0:kk,1,0)+drva1(0:kk,1,1)+drva1(0:kk,2,0)+drva1(0:kk,2,1))
! end if
 do np=1,4
    call MPI_SEND(drva1(0:kk,3,0),kk+1,MPI_REAL8,mjct(np),itag,icom,ierr)
 end do
 end if
 end do

 end subroutine junction

!===== INTERMEDIATE RESULTS

 subroutine intermed(no)
 integer :: no

select case(no)
case(1);rr(:,1)=qa(:,1)
case(2);rr(:,1)=sqrt(qa(:,2)*qa(:,2)+qa(:,3)*qa(:,3)+qa(:,4)*qa(:,4))/qa(:,1)
case(3);rr(:,1)=p(:)
case(4);rr(:,1)=gam*p(:)-1
case(5);rr(:,1)=gamm1*(qo(:,5)-half*(qo(:,2)*qo(:,2)+qo(:,3)*qo(:,3)+qo(:,4)*qo(:,4))/qo(:,1))
end select


 end subroutine intermed

!===== FINAL RESULTS

 subroutine finalout

 real(nr),dimension(:),allocatable :: delt
 
    ns=0; ne=ndata; allocate(delt(ns:ne))
    fctr=half/(times(ne)-times(ns))
    delt(ns)=fctr*(times(ns+1)-times(ns)); delt(ne)=fctr*(times(ne)-times(ne-1))
 do n=ns+1,ne-1
    delt(n)=fctr*(times(n+1)-times(n-1))
 end do
    rr(:,1)=0
 do n=0,ndata
    read(0,rec=n+4) varr(:)
    rr(:,1)=rr(:,1)+delt(n)*varr(:)
 end do
    ss(:,1)=0; fctr=1
 do n=0,ndata
    read(0,rec=n+4) varr(:)
    varr(:)=fctr*(varr(:)-rr(:,1))
    write(0,rec=n+4) varr(:)
    ss(:,1)=ss(:,1)+delt(n)*varr(:)*varr(:)
 end do
 if(ndatp==1) then
    varr=sqrt(ss(:,1))
    write(0,rec=ndata+5) varr(:)
 end if

 end subroutine finalout

!=====

 end module problemcase

!*****
