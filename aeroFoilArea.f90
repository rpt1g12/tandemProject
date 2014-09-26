MODULE aeroFoilArea

PRIVATE
PUBLIC :: setupFWH,writePressData
INTEGER,PARAMETER          :: permTen(3,3) = RESHAPE( (/1,2,3,2,3,1,3,1,2/), SHAPE(permTen) ) 
CONTAINS
!================================
SUBROUTINE calcArea 

USE MPI
USE MainVar3D,   ONLY :  mpro,nr,zero,yaco,etm,nbc,myid,ijk,&
                         one,two,half,nviscous,permTen,nbsize,&
                         icom,qo,qa,de
USE Subroutineso,ONLY :  indx3
IMPLICIT NONE
REAL(KIND = nr):: tmp,threedVect(3),coef,tmp2,g11,g33,g13
INTEGER        :: ip, nn, mm, l, i, ierr, np, j, k, jp


tmp = zero
tmp2= zero
DO nn=1,3
  DO ip=0,1
    IF(nbc(ip,nn)==20+5*nviscous)THEN
      i=ip*ijk(1,nn)
      tmp = zero
      !--Corner Treatment
      np = 1
      l = indx3(i,0,0,nn)
      coef = calc_coefc(nbc,permTen,nn,0,0)
      g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
      g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
      g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
      tmp2= tmp2 + coef*sqrt(g11*g33-g13**2)
      
      np = np + 1
      l = indx3(i,ijk(2,nn),0,nn)
      coef = calc_coefc(nbc,permTen,nn,1,0)
      g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
      g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
      g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
      tmp2= tmp2 + coef*sqrt(g11*g33-g13**2)      

      np = np + 1
      l = indx3(i,0,ijk(3,nn),nn)
      coef = calc_coefc(nbc,permTen,nn,0,1)
      g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
      g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
      g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
      tmp2= tmp2 + coef*sqrt(g11*g33-g13**2)      

      np = np + 1
      l = indx3(i,ijk(2,nn),ijk(3,nn),nn)
      coef = calc_coefc(nbc,permTen,nn,1,1)
      g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
      g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
      g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
      tmp2= tmp2 + coef*sqrt(g11*g33-g13**2)      

      !--Edge Treatment
      DO jp = 0,1; j = jp*ijk(2,nn)
        DO k = 1,ijk(3,nn)-1
          np = np + 1
          l = indx3(i,j,k,nn)
          coef = calc_coefe(nbc,jp,permTen(2,nn))
          g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
          g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
          g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
          tmp2= tmp2 + coef*sqrt(g11*g33-g13**2)          
        ENDDO
      ENDDO

      DO jp = 0,1; k = jp*ijk(3,nn)
        DO j = 1,ijk(2,nn)-1
          np = np + 1
          l = indx3(i,j,k,nn)
          coef = calc_coefe(nbc,jp,permTen(3,nn))
          g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
          g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
          g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
          tmp2= tmp2 + coef*sqrt(g11*g33-g13**2)          
        ENDDO
      ENDDO

      !--Internal node treatment
      DO k=1,ijk(3,nn)-1
        DO j=1,ijk(2,nn)-1
          np = np + 1
          l = indx3(i,j,k,nn)
          g11=qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
          g33=qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
          g13=qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
          tmp2= tmp2 + sqrt(g11*g33-g13**2) 
        ENDDO
      ENDDO
         
    ENDIF
  ENDDO
ENDDO

WRITE(*,'(A,F12.8,A,I4)')'Calculated surf area 2: ', tmp2, ', On processor: ', myid
CALL MPI_ALLREDUCE(MPI_IN_PLACE,tmp2,1,MPI_REAL8,MPI_SUM,icom,ierr)
IF(myid == 0)WRITE(*,'(A,F12.8)')'Calculated total surf area for FWH is: ', tmp2

END SUBROUTINE calcArea

FUNCTION calc_coefc(nbc,permTen,nn,ip1,ip2)

USE MainVar3D,    ONLY : nr,half,quarter,one

IMPLICIT NONE

INTEGER,INTENT(IN)  :: nbc(0:1,3),permTen(3,3),ip1,ip2,nn
REAL(nr)            :: calc_coefc

calc_coefc = one
SELECT CASE( nbc(ip1,permTen(2,nn)) )
  CASE(35)
    SELECT CASE( nbc(ip2,permTen(3,nn)) )
      CASE(35,45)
        calc_coefc = quarter
      CASE(40)
        calc_coefc = half
    END SELECT
  CASE(40)
    SELECT CASE( nbc(ip2,permTen(3,nn)) )
      CASE(35,45)
        calc_coefc = half
      CASE(40)
        calc_coefc = one
    END SELECT
  CASE(45)
    SELECT CASE( nbc(ip2,permTen(3,nn)) )
      CASE(35,45)
        calc_coefc = quarter
      CASE(40)
        calc_coefc = half
    END SELECT
  CASE DEFAULT
    WRITE(*,*)'Boundary not found. Setting coefficient to one...'
    calc_coefc = one
END SELECT

END FUNCTION calc_coefc

FUNCTION calc_coefe(nbc,ip1,ip2)

USE MainVar3D,    ONLY : nr,half,quarter,one

IMPLICIT NONE

INTEGER,INTENT(IN)  :: nbc(0:1,3),ip1,ip2
REAL(nr)            :: calc_coefe


SELECT CASE(nbc(ip1,ip2))
  CASE(35,45)
    calc_coefe = half
  CASE(40)
    calc_coefe = one
  CASE DEFAULT
    WRITE(*,*)'Boundary not found. Setting coefficient to one...'
    calc_coefe = one
END SELECT

END FUNCTION calc_coefe

END MODULE aeroFoilArea