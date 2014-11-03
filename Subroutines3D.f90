!*****
!***** SUBROUTINES FOR 3D SOLVER
!*****

 module subroutines3d

 use mpi
 use mainvar3d
 use subroutineso
 implicit none

 contains

!===== SUBROUTINE FOR MPI IMPLEMENTATION

 subroutine mpigo(nt,nrt,n45,itag)

 integer,intent(in) :: nt,nrt,n45,itag

 select case(nt); case(0); mp=lmd; case(1); mp=lmf; end select

    ir=0
 do nn=1,3; nz=(1-nrt)*(nn-1)+1
 select case(nn)
 case(1); send=>send1; recv=>recv1; case(2); send=>send2; recv=>recv2; case(3); send=>send3; recv=>recv3
 end select
 do ip=0,1; iq=1-ip; is=ip*ijk(1,nn); ie=1-2*ip
 select case(nbc(ip,nn)); case(35); np=0; ii=1; case(40); np=0; ii=0; case(45); np=n45; ii=1; end select
 if(ndf(ip,nt,nn)==1) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(is,j,k,nn)
    res=np*rr(l,nz)
 do i=0,mp; l=indx3(is+ie*(i+ii),j,k,nn)
    sap(i)=rr(l,nz)
 end do
    send(jk,0,ip)=sum(pbco(0:mp,0,nt)*sap(0:mp))-res*pbcot(0,nt)
    send(jk,1,ip)=sum(pbco(0:mp,1,nt)*sap(0:mp))-res*pbcot(1,nt)
    send(jk,2,ip)=sap(0)-res
 end do
 end do
    ir=ir+1; call MPI_ISEND(send(:,:,ip),3*nbsize(nn),MPI_REAL8,ncd(ip,nn),itag+iq,icom,ireq(ir),ierr)
    ir=ir+1; call MPI_IRECV(recv(:,:,ip),3*nbsize(nn),MPI_REAL8,ncd(ip,nn),itag+ip,icom,ireq(ir),ierr)
 end if
 end do
 end do
 if(ir/=0) then
    call MPI_WAITALL(ir,ireq,ista,ierr)
 end if

 if(n45==n45go) then
 do nn=1,3; nz=(1-nrt)*(nn-1)+1
 select case(nn); case(1); recv=>recv1; case(2); recv=>recv2; case(3); recv=>recv3; end select
 do ip=0,1; is=ip*ijk(1,nn)
 if(nbc(ip,nn)==45) then
 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j; l=indx3(is,j,k,nn)
    recv(jk,0,ip)=recv(jk,0,ip)+rr(l,nz)*pbcot(0,nt)
    recv(jk,1,ip)=recv(jk,1,ip)+rr(l,nz)*pbcot(1,nt)
    recv(jk,2,ip)=recv(jk,2,ip)+rr(l,nz)
 end do
 end do
 end if
 end do
 end do
 end if

 end subroutine mpigo

!===== SUBROUTINE FOR COMPACT FINITE DIFFERENTIATING

 subroutine deriv(nn,nz)

 integer,intent(in) :: nn,nz

    nt=0; ns=ndf(0,0,nn); ne=ndf(1,0,nn)

 select case(nn)
 case(1); is=0; ie=is+lxi; recv=>recv1
 case(2); is=lxi+1; ie=is+let; recv=>recv2
 case(3); is=lxi+let+2; ie=is+lze; recv=>recv3
 end select

 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j
 do i=is,ie; l=indx3(i-is,j,k,nn)
    li(i)=l; sa(i)=rr(l,nz)
 end do
 select case(ns)
 case(0)
    sb(is)=dot_product(abc(:,0),(/sa(is+1),sa(is+2),sa(is+3),sa(is+4),sa(is+5),sa(is+6)/)-sa(is))
    sb(is+1)=dot_product(abc(:,1),(/sa(is),sa(is+2),sa(is+3),sa(is+4),sa(is+5),sa(is+6)/)-sa(is+1))
    sb(is+2)=dot_product(abc(:,2),(/sa(is),sa(is+1),sa(is+3),sa(is+4),sa(is+5),sa(is+6)/)-sa(is+2))
 case(1)
    sb(is)=sum(pbci(0:lmd,0,nt)*sa(is:is+lmd))+recv(jk,0,0)
    sb(is+1)=sum(pbci(0:lmd,1,nt)*sa(is:is+lmd))+recv(jk,1,0)
    sb(is+2)=aa*(sa(is+3)-sa(is+1))+ab*(sa(is+4)-sa(is))+ac*(sa(is+5)-recv(jk,2,0))
 end select
 do i=is+3,ie-3
    sb(i)=aa*(sa(i+1)-sa(i-1))+ab*(sa(i+2)-sa(i-2))+ac*(sa(i+3)-sa(i-3))
 end do
 select case(ne)
 case(0)
    sb(ie)=dot_product(abc(:,0),sa(ie)-(/sa(ie-1),sa(ie-2),sa(ie-3),sa(ie-4),sa(ie-5),sa(ie-6)/))
    sb(ie-1)=dot_product(abc(:,1),sa(ie-1)-(/sa(ie),sa(ie-2),sa(ie-3),sa(ie-4),sa(ie-5),sa(ie-6)/))
    sb(ie-2)=dot_product(abc(:,2),sa(ie-2)-(/sa(ie),sa(ie-1),sa(ie-3),sa(ie-4),sa(ie-5),sa(ie-6)/))
 case(1)
    sb(ie)=-sum(pbci(0:lmd,0,nt)*sa(ie:ie-lmd:-1))-recv(jk,0,1)
    sb(ie-1)=-sum(pbci(0:lmd,1,nt)*sa(ie:ie-lmd:-1))-recv(jk,1,1)
    sb(ie-2)=aa*(sa(ie-1)-sa(ie-3))+ab*(sa(ie)-sa(ie-4))+ac*(recv(jk,2,1)-sa(ie-5))
 end select
    sa(is)=sb(is)
    sa(is+1)=sb(is+1)-xl(is+1,2)*sa(is)
 do i=is+2,ie
    sa(i)=sb(i)-xl(i,1)*sa(i-2)-xl(i,2)*sa(i-1)
 end do
    sb(ie)=xu(ie,1)*sa(ie)
    sb(ie-1)=xu(ie-1,1)*sa(ie-1)-xu(ie-1,2)*sb(ie)
 do i=ie-2,is,-1
    sb(i)=xu(i,1)*sa(i)-xu(i,2)*sb(i+1)-xu(i,3)*sb(i+2)
 end do
 do i=is,ie
    l=li(i); rr(l,nn)=sb(i)
 end do
 end do
 end do

 end subroutine deriv

!===== SUBROUTINE FOR COMPACT FILTERING

 subroutine filte(nn,nz)

 integer,intent(in) :: nn,nz

    nt=1; ns=ndf(0,1,nn); ne=ndf(1,1,nn)

 select case(nn)
 case(1); is=0; ie=is+lxi; recv=>recv1
 case(2); is=lxi+1; ie=is+let; recv=>recv2
 case(3); is=lxi+let+2; ie=is+lze; recv=>recv3
 end select

 do k=0,ijk(3,nn); kp=k*(ijk(2,nn)+1)
 do j=0,ijk(2,nn); jk=kp+j
 do i=is,ie; l=indx3(i-is,j,k,nn)
    li(i)=l; sa(i)=rr(l,nz)
 end do
 select case(ns)
 case(0); ra0=2*sa(is); ra1=2*sa(is+1); ra2=2*sa(is+2)
    res=sum(fex(:)*(sa(is+mfbi*(/1,2,3/))-sa(is))); rv(:)=sa(is)-res*(/1,2,3/)
    sb(is)=fam(0)*(rv(1)+sa(is+1)-ra0)+fbm(0)*(rv(2)+sa(is+2)-ra0)+fcm(0)*(rv(3)+sa(is+3)-ra0)
    sb(is+1)=fam(1)*(sa(is)+sa(is+2)-ra1)+fbm(1)*(rv(1)+sa(is+3)-ra1)+fcm(1)*(rv(2)+sa(is+4)-ra1)
    sb(is+2)=fam(2)*(sa(is+1)+sa(is+3)-ra2)+fbm(2)*(sa(is)+sa(is+4)-ra2)+fcm(2)*(rv(1)+sa(is+5)-ra2)
 case(1); ra2=2*sa(is+2)
    sb(is)=sum(pbci(0:lmf,0,nt)*sa(is:is+lmf))+recv(jk,0,0)
    sb(is+1)=sum(pbci(0:lmf,1,nt)*sa(is:is+lmf))+recv(jk,1,0)
    sb(is+2)=fa*(sa(is+1)+sa(is+3)-ra2)+fb*(sa(is)+sa(is+4)-ra2)+fc*(recv(jk,2,0)+sa(is+5)-ra2)
 end select
 do i=is+3,ie-3
    res=2*sa(i); sb(i)=fa*(sa(i-1)+sa(i+1)-res)+fb*(sa(i-2)+sa(i+2)-res)+fc*(sa(i-3)+sa(i+3)-res)
 end do
 select case(ne)
 case(0); ra0=2*sa(ie); ra1=2*sa(ie-1); ra2=2*sa(ie-2)
    res=sum(fex(:)*(sa(ie)-sa(ie-mfbi*(/1,2,3/)))); rv(:)=sa(ie)+res*(/1,2,3/)
    sb(ie)=fam(0)*(sa(ie-1)+rv(1)-ra0)+fbm(0)*(sa(ie-2)+rv(2)-ra0)+fcm(0)*(sa(ie-3)+rv(3)-ra0)
    sb(ie-1)=fam(1)*(sa(ie-2)+sa(ie)-ra1)+fbm(1)*(sa(ie-3)+rv(1)-ra1)+fcm(1)*(sa(ie-4)+rv(2)-ra1)
    sb(ie-2)=fam(2)*(sa(ie-3)+sa(ie-1)-ra2)+fbm(2)*(sa(ie-4)+sa(ie)-ra2)+fcm(2)*(sa(ie-5)+rv(1)-ra2)
 case(1); ra2=2*sa(ie-2)
    sb(ie)=sum(pbci(0:lmf,0,nt)*sa(ie:ie-lmf:-1))+recv(jk,0,1)
    sb(ie-1)=sum(pbci(0:lmf,1,nt)*sa(ie:ie-lmf:-1))+recv(jk,1,1)
    sb(ie-2)=fa*(sa(ie-3)+sa(ie-1)-ra2)+fb*(sa(ie-4)+sa(ie)-ra2)+fc*(sa(ie-5)+recv(jk,2,1)-ra2)
 end select
    sa(is)=sb(is)
    sa(is+1)=sb(is+1)-yl(is+1,2)*sa(is)
 do i=is+2,ie
    sa(i)=sb(i)-yl(i,1)*sa(i-2)-yl(i,2)*sa(i-1)
 end do
    sb(ie)=yu(ie,1)*sa(ie)
    sb(ie-1)=yu(ie-1,1)*sa(ie-1)-yu(ie-1,2)*sb(ie)
 do i=ie-2,is,-1
    sb(i)=yu(i,1)*sa(i)-yu(i,2)*sb(i+1)-yu(i,3)*sb(i+2)
 end do
 do i=is,ie
    l=li(i); rr(l,nz)=rr(l,nz)+sb(i)
 end do
 end do
 end do

 end subroutine filte

!===== SUBROUTINE FOR ELEMENTARY VARIABLES IN GCBC/GCIC

 subroutine eleme(l,cm)

 integer,intent(in) :: l
 real(nr),dimension(3),intent(in) :: cm

    rhoi=1/qa(l,1); ao=sqrt(gam*rhoi*p(l)); aoi=1/ao
    ve(:)=rhoi*qa(l,2:4); hv2=half*(ve(1)*ve(1)+ve(2)*ve(2)+ve(3)*ve(3))
    vn=cm(1)*ve(1)+cm(2)*ve(2)+cm(3)*ve(3); vs=cm(1)*umf(1)+cm(2)*umf(2)+cm(3)*umf(3)

 end subroutine eleme

!===== SUBROUTINE FOR TRANSFORMATION FROM Q TO R IN GCBC/GCIC

 subroutine xtq2r(cm)

 real(nr),dimension(3),intent(in) :: cm

    ho=gamm1*aoi*aoi; bo=1-ho*hv2; co=aoi*vn; dm(:)=aoi*cm(:); rv(:)=ho*ve(:)

    xt(1,1)=bo*cm(1)+dm(2)*ve(3)-dm(3)*ve(2)
    xt(1,2)=cm(1)*rv(1)
    xt(1,3)=cm(1)*rv(2)+dm(3)
    xt(1,4)=cm(1)*rv(3)-dm(2)
    xt(1,5)=-ho*cm(1)

    xt(2,1)=bo*cm(2)+dm(3)*ve(1)-dm(1)*ve(3)
    xt(2,2)=cm(2)*rv(1)-dm(3)
    xt(2,3)=cm(2)*rv(2)
    xt(2,4)=cm(2)*rv(3)+dm(1)
    xt(2,5)=-ho*cm(2)

    xt(3,1)=bo*cm(3)+dm(1)*ve(2)-dm(2)*ve(1)
    xt(3,2)=cm(3)*rv(1)+dm(2)
    xt(3,3)=cm(3)*rv(2)-dm(1)
    xt(3,4)=cm(3)*rv(3)
    xt(3,5)=-ho*cm(3)

    xt(4,1)=1-bo-co
    xt(4,2)=dm(1)-rv(1)
    xt(4,3)=dm(2)-rv(2)
    xt(4,4)=dm(3)-rv(3)
    xt(4,5)=ho

    xt(5,1)=1-bo+co
    xt(5,2)=-dm(1)-rv(1)
    xt(5,3)=-dm(2)-rv(2)
    xt(5,4)=-dm(3)-rv(3)
    xt(5,5)=ho

 end subroutine xtq2r

!===== SUBROUTINE FOR INVERSE TRANSFORMATION FROM R TO Q IN GCBC/GCIC

 subroutine xtr2q(cm)

 real(nr),dimension(3),intent(in) :: cm

    bo=hv2+hamm1*ao*ao; co=ao*vn; dm(:)=ao*cm(:)

    xt(1,1)=cm(1)
    xt(1,2)=cm(2)
    xt(1,3)=cm(3)
    xt(1,4)=half
    xt(1,5)=half

    xt(2,1)=cm(1)*ve(1)
    xt(2,2)=cm(2)*ve(1)-dm(3)
    xt(2,3)=cm(3)*ve(1)+dm(2)
    xt(2,4)=half*(ve(1)+dm(1))
    xt(2,5)=xt(2,4)-dm(1)

    xt(3,1)=cm(1)*ve(2)+dm(3)
    xt(3,2)=cm(2)*ve(2)
    xt(3,3)=cm(3)*ve(2)-dm(1)
    xt(3,4)=half*(ve(2)+dm(2))
    xt(3,5)=xt(3,4)-dm(2)

    xt(4,1)=cm(1)*ve(3)-dm(2)
    xt(4,2)=cm(2)*ve(3)+dm(1)
    xt(4,3)=cm(3)*ve(3)
    xt(4,4)=half*(ve(3)+dm(3))
    xt(4,5)=xt(4,4)-dm(3)

    xt(5,1)=hv2*cm(1)+dm(3)*ve(2)-dm(2)*ve(3)
    xt(5,2)=hv2*cm(2)+dm(1)*ve(3)-dm(3)*ve(1)
    xt(5,3)=hv2*cm(3)+dm(2)*ve(1)-dm(1)*ve(2)
    xt(5,4)=half*(bo+co)
    xt(5,5)=xt(5,4)-co

 end subroutine xtr2q

!=====

 end module subroutines3d

!*****
