MODULE interpVar

!USE Ray_Casting_Algo
!USE Points_Module
!USE Polygons
INTEGER,PARAMETER             :: IOLenUnit = 4
INTEGER,PARAMETER             :: DBLELenUnit = 2  
PRIVATE
PUBLIC :: findMNearestPoints2D,createInterpPoints,interpIDW,writeSurfPress

CONTAINS

SUBROUTINE writeSurfPress(varOut,nStep,finalize)

!--Assuming plate is on the grid block of 1 and 4
!--
USE MainVar3D,      ONLY: ijk,myid,mo,mb,npc,nsampSurf,lxim,lzem,lximb,lzemb,nr,timo,pSurfFileUnit,&
                          writeVel,writeOm,lmx,writeRho
USE Subroutineso,   ONLY: indx3
IMPLICIT NONE 

INTEGER,INTENT(IN)  :: nStep
LOGICAL,INTENT(IN)  :: finalize
REAL(nr),INTENT(IN) :: varOut(0:lmx,5)
INTEGER             :: bsize,i,k,j,mp,mmp,l,ll,maxInd
INTEGER,PARAMETER   :: topBlock = 4, bottBlock = 1
INTEGER,SAVE        :: cpos,hpos,lpos,ip,jp,kp,ltos
REAL(nr),ALLOCATABLE,SAVE:: buff(:,:)

bsize = (ijk(2,2)+1)*(ijk(3,2)+1)

cpos = 0

maxInd = writevel+writeRho+1

IF(.NOT. finalize)THEN
  IF(nStep == 0)THEN
    ip=mod(myid-mo(mb),npc(mb,1))
    jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))
    ltos = maxInd*(lximb(mb)+1)*(lzemb(mb)+1)
    ALLOCATE(buff(bsize,maxInd))
    IF(mb == bottBlock .AND. jp == npc(mb,2)-1 )THEN
      OPEN(pSurfFileUnit,file='pbSurfHist.dat',access='stream',form='unformatted',shared)
      IF(ip == 0 .AND. kp == 0)WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)nsampSurf,writeVel,writeOm,writeRho
      cpos = 4+3*(kp*npc(mb,1)+ip)
      WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)bsize,ijk(3,2),ijk(2,2)
      hpos = 4+3*(npc(mb,1)*npc(mb,3))
      lpos = 0
      DO k = 0,kp-1
        DO i = 0,npc(mb,1)-1
           mp=mo(mb)+k*npc(mb,1)*npc(mb,2)+jp*npc(mb,1)+i
           lpos = lpos + maxInd*(lxim(mp)+1)*(lzem(mp)+1)
        ENDDO
      ENDDO
      DO i = 0,ip-1
        mp=mo(mb)+kp*npc(mb,1)*npc(mb,2)+jp*npc(mb,1)+i
        lpos = lpos + maxInd*(lxim(mp)+1)*(lzem(mp)+1)
      ENDDO
    ELSEIF(mb == topBlock .AND. jp == 0 )THEN
      OPEN(pSurfFileUnit,file='ptSurfHist.dat',access='stream',form='unformatted',shared)
      IF(ip == 0 .AND. kp == 0)WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)nsampSurf,writeVel,writeOm,writeRho
      cpos = 4+3*(kp*npc(mb,1)+ip)
      WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)bsize,ijk(3,2),ijk(2,2)
      hpos = 4+3*(npc(mb,1)*npc(mb,3))  
      lpos = 0
      DO k = 0,kp-1
        DO i = 0,npc(mb,1)-1
           mp=mo(mb)+k*npc(mb,1)*npc(mb,2)+jp*npc(mb,1)+i
           lpos = lpos + maxInd*(lxim(mp)+1)*(lzem(mp)+1)
        ENDDO
      ENDDO
      DO i = 0,ip-1
        mp=mo(mb)+kp*npc(mb,1)*npc(mb,2)+jp*npc(mb,1)+i
        lpos = lpos + maxInd*(lxim(mp)+1)*(lzem(mp)+1)
      ENDDO
    ENDIF 
         
  ELSE
    IF(mb == bottBlock .AND. jp == npc(mb,2)-1 )THEN
      l = 0; j = ijk(1,2)
      DO k = 0,ijk(2,2)
        DO i = 0,ijk(3,2)
          l = l+1
          ll= indx3(j,k,i,2)
          buff(l,1:maxInd)=varOut(ll,1:maxInd)
        ENDDO
      ENDDO
      cpos = hpos+DBLELenUnit*((nStep-1)*(ltos+1)+lpos)
      IF(ip == 0 .AND. kp == 0)WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)timo
      cpos = cpos+DBLELenUnit
      WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)buff
    ELSEIF(mb == topBlock .AND. jp == 0 )THEN
      l = 0; j = 0
      DO k = 0,ijk(2,2)
        DO i = 0,ijk(3,2)
          l = l+1
          ll= indx3(j,k,i,2)
          buff(l,1:maxInd) = varOut(ll,1:maxInd)
        ENDDO
      ENDDO
      cpos = hpos+DBLELenUnit*((nStep-1)*(ltos+1)+lpos)
      IF(ip == 0 .AND. kp == 0)WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)timo
      cpos = cpos+DBLELenUnit
      WRITE(pSurfFileUnit,pos = IOLenUnit*cpos+1)buff
    ENDIF
  ENDIF
ELSE
  CLOSE(pSurfFileUnit)
  DEALLOCATE(buff)
ENDIF

END SUBROUTINE writeSurfPress

SUBROUTINE createInterpPoints(x,y,z,xp,yp,zp,interpInd,order,np,Rc,zl)
!--Create a circular region on x-y plane with z=z(zConst)
   
USE MPI
USE MainVar3D,     ONLY : nr,two,pi,zero,half,npro,ierr,myid,ijk,six,three,sml,icom,lmx,lze,narc,snp,npro,npc,mb,mo,lzemb
USE Subroutineso,  ONLY : indx3
USE problemcase,   ONLY : domlen,wlew,wlea,span,smgrid

IMPLICIT NONE 
REAL(KIND = nr)               :: x(0:lmx),y(0:lmx),z(0:lmx)
REAL(KIND = nr),ALLOCATABLE   :: xp(:,:),yp(:,:),zp(:,:)    !--Can be unallocated initially, location of the points
INTEGER                       :: np(narc)             !--number of intepolation points on the circ  
                                                      !--On return this will be the number of point on the current processor
                                                      !--zConst is the index of a z-coordinate of the x-y plane 
INTEGER                       :: interpInd(narc,npro),order(narc,npro) !--Accumulative number of interpolation points on each processor
REAL(KIND = nr)               :: Rc(narc),zl                !--Radius of the circle

REAL(KIND = nr)               :: nnp(narc),maxx,minx,maxy,miny,XX(4),YY(4),p0x,p0y,p0z,maxz,minz,sdist(npro),cdist
INTEGER                       :: mm,ll,i,l(narc),nInt,ns(narc,npro),j,k,zConst,inds(npro),kp,ls,le !,vt(4*(let+lxi))
!TYPE(point)                   :: pts(2*(let+lxi)),cp
!TYPE(polygon)                 :: poly

IF(ALLOCATED(xp))DEALLOCATE(xp)
IF(ALLOCATED(yp))DEALLOCATE(yp)
IF(ALLOCATED(zp))DEALLOCATE(zp)

ALLOCATE(xp(snp,narc),yp(snp,narc),zp(snp,narc))

maxx = MAXVAL(x)
minx = MINVAL(x)
maxy = MAXVAL(y)
miny = MINVAL(y)
maxz = MAXVAL(z)
minz = MINVAL(z)

zConst = -1
ns(:,myid+1) = np+1
nnp = REAL(np,nr)

kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

DO k=0,lzemb(mb)
  IF(zl+sml>span*(real(lzemb(mb)-k,nr)/lzemb(mb)-half))then
    mm=k
    EXIT
  ENDIF
ENDDO

ls=0;le=0
DO k=1,npc(mb,3)
  IF(k==1)THEN
    le=le+lzemb(mb)-((lzemb(mb)+1)/npc(mb,3))*(npc(mb,3)-1)
  ELSE
    le=le+(lzemb(mb)+1)/npc(mb,3)
  ENDIF
  IF(mm>=ls.AND.mm<=le)EXIT
  ls=le+1
ENDDO
IF(k-1==kp)zConst=mm-ls
!
!  DO i=0,lze-1
!    ll = indx3(1,1,i,1)
!    IF(zl - sml > z(ll))THEN
!      ll = indx3(1,1,i,1)
!      IF(zl - sml < z(ll))THEN
!        zConst = i
!        EXIT
!      ENDIF
!    ENDIF
!  ENDDO   
!   
IF(zconst>=0) p0z = z(indx3(1,1,zConst,1))
IF(zconst==-1)p0z = z(indx3(1,1,0,1))

  DO k=1,narc
    IF(k<narc)THEN
      l(k)=0
      DO i=0,np(k)-1
        p0x = Rc(k)*COS(pi*i/(nnp(k)-1)) 
        p0y = Rc(k)*SIN(pi*i/(nnp(k)-1))
  !      IF(ABS(p0x) < 1d-15)p0x = 0
  !      IF(ABS(p0y) < 1d-15)p0y = 0
         sdist = domlen
         IF(zconst>=0)THEN
           DO j = 0,lmx
             cdist = SQRT( (x(j) - p0x)**2 + (y(j) - p0y)**2 + (z(j) - p0z)**2 ) !--Calculate Distance
             IF(cdist < sdist(myid+1))THEN !--Find the position of the current point in the list
               sdist(myid+1) = cdist
             ENDIF 
           ENDDO
         ENDIF
         CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,sdist,1,MPI_REAL8,icom,ierr)
         IF(MINLOC(sdist,DIM=1) == myid+1)THEN
           IF(l(k) == 0)ns(k,myid+1) = i+1
           l(k) = l(k) + 1
           xp(l(k),k) = p0x
           yp(l(k),k) = p0y
           zp(l(k),k) = p0z
         ENDIF
      ENDDO   
    ELSE
      l(k)=0
      DO i=0,np(k)-1
        p0x = -half-(wlea+2*smgrid) 
        p0y = two*EPSILON(half)
        p0z = span*(real(np(k)-i,nr)/np(k)-half)
  !      IF(ABS(p0x) < 1d-15)p0x = 0
  !      IF(ABS(p0y) < 1d-15)p0y = 0
         sdist = domlen
         DO j = 0,lmx
           cdist = SQRT( (x(j) - p0x)**2 + (y(j) - p0y)**2 + (z(j) - p0z)**2 ) !--Calculate Distance
           IF(cdist < sdist(myid+1))THEN !--Find the position of the current point in the list
             sdist(myid+1) = cdist
           ENDIF 
         ENDDO
         CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,sdist,1,MPI_REAL8,icom,ierr)
         IF(MINLOC(sdist,DIM=1) == myid+1)THEN
           IF(l(k) == 0)ns(k,myid+1) = i+1
           l(k) = l(k) + 1
           xp(l(k),k) = p0x
           yp(l(k),k) = p0y
           zp(l(k),k) = p0z
         ENDIF
      ENDDO
    ENDIF
  ENDDO
  
  np = l
  CALL MPI_ALLGATHER(np,narc,MPI_INTEGER,interpInd,narc,MPI_INTEGER,icom,ierr)
  call MPI_ALLREDUCE(np,l,narc,MPI_INTEGER,MPI_SUM,icom,ierr)
  IF(myid == 0)THEN
    DO mm=1,narc
        WRITE(*,*)'Total interpolation points found on arc: ', mm, 'is: ', l(mm),', Total Requested: ',nnp(mm)
    ENDDO
  ENDIF
  CALL MPI_ALLGATHER(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,ns,narc,MPI_INTEGER,icom,ierr)

DO mm=1,narc
  
  order(mm,1:npro) = (/1:npro/)

  !--Stupid sorting (not important since number of procs is small)
  DO j = 1,npro-1
    DO k = j+1,npro
      IF(ns(mm,j) > ns(mm,k))THEN
        ll = ns(mm,j);ns(mm,j) = ns(mm,k);ns(mm,k) = ll
        ll = order(mm,j);order(mm,j) = order(mm,k);order(mm,k) = ll
      ENDIF
    ENDDO
  ENDDO

  DO i=2,npro
    interpInd(mm,order(mm,i)) = interpInd(mm,order(mm,i-1))+interpInd(mm,order(mm,i))
  ENDDO
ENDDO

END SUBROUTINE createInterpPoints

SUBROUTINE findMNearestPoints2D(x,y,z,xp,yp,zp,np,MNearest,indM,distM)
!--Very simple search algorithm to find nearest MNearest grid points to np interpolation point. 
USE MainVar3D,     ONLY : nr,lmx,lxi,let,lze,ijk,narc,snp
USE Subroutineso,  ONLY : indx3
IMPLICIT NONE

INTEGER,INTENT(IN)            :: np         !--number of interpolation points
INTEGER,INTENT(IN)            :: MNearest          !--Number of nearest point to be used in 
                                                   !--the interpolation for each point in {1..m}
REAL(KIND = nr),INTENT(IN)    :: x(0:lmx),y(0:lmx),z(0:lmx) !--coordinates of n data points
REAL(KIND = nr),INTENT(IN)    :: xp(snp),yp(snp),zp(snp) !--coordinates of mp interpolation points
INTEGER,INTENT(INOUT)         :: indM(MNearest,snp)    !--Index of mp nearest points in {1..MNearest} 
                                                      !--for each interpolation point in {1..m}
REAL(KIND = nr)               :: distM(MNearest,snp)   !--distance to mp nearest grid points
INTEGER                       :: m,l,p1,k
REAL(KIND = nr)               :: cdist

distM(:,:) = 1.D12
  DO m = 1,np
    DO l = 0,lmx
      cdist = SQRT( (x(l) - xp(m))**2 + (y(l) - yp(m))**2 + (z(l) - zp(m))**2 ) !--Calculate Distance
      DO p1 = 1,MNearest
        IF(cdist < distM(p1,m))THEN !--Find the position of the current point in the list
          distM(p1:MNearest,m)=EOSHIFT(distM(p1:MNearest,m),SHIFT = -1)
          indM(p1:MNearest,m)=EOSHIFT(indM(p1:MNearest,m),SHIFT = -1)
          distM(p1,m) = cdist  !--set the current value
          indM(p1,m) = l
          EXIT
        ENDIF 
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE findMNearestPoints2D

SUBROUTINE interpIDW(distM,varM,interpVar,MNearest,np,nvar)
!--INterpolation using Inverse Distance Weighting method
USE MainVar3D,     ONLY : nr,one,zero,narc,snp

IMPLICIT NONE

INTEGER,INTENT(IN)            :: np,nvar           !--number of interpolation points
INTEGER,INTENT(IN)            :: MNearest          !--Number of nearest point to be used in 
                                                   !--the interpolation for each point in {1..m}
REAL(KIND = nr)               :: distM(MNearest,snp),varM(nvar,MNearest,snp),interpVar(nvar,snp)

INTEGER         :: i,j
REAL(KIND = nr) :: dummy

interpVar = zero

  DO i = 1,np
    dummy = zero
    DO j = 1,MNearest
      IF(.NOT. ABS(distM(j,i)) > zero)THEN
        interpVar(:,i) = varM(:,j,i)
        dummy = one
        EXIT
      ENDIF
      interpVar(:,i) = interpVar(:,i) + varM(:,j,i)/distM(j,i)
      dummy = dummy + one/distM(j,i)
    ENDDO
    interpVar(:,i) = interpVar(:,i)/dummy
  ENDDO


END SUBROUTINE interpIDW

SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT) 
  USE MainVar3D,  ONLY : nr
  IMPLICIT NONE
  INTEGER         :: N
  REAL(KIND = nr) :: XX(N),YY(N),PX,PY 
  LOGICAL         :: MX,MY,NX,NY 
  REAL(KIND = nr) :: X(N),Y(N)
  INTEGER         :: I,J,INOUT
  !OUTPUT UNIT FOR PRINTED MESSAGES 
  
  DO 1 I=1,N 
    X(I)=XX(I)-PX 
  1 Y(I)=YY(I)-PY 
    INOUT=-1 
    DO 2 I=1,N 
      J=1+MOD(I,N) 
      MX=X(I).GE.0.0 
      NX=X(J).GE.0.0 
      MY=Y(I).GE.0.0 
      NY=Y(J).GE.0.0 
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2 
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3 
      INOUT=-INOUT 
      GO TO 2 
    3 IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5 
    4 INOUT=0 
      RETURN 
    5 INOUT=-INOUT 
    2 CONTINUE 
      RETURN 
END SUBROUTINE PNPOLY

ENDMODULE interpVar
