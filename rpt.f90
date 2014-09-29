!**********************
!***** RPT MODULE *****
!**********************

module  rpt

use mainvar3d
use subroutineso
use mpi

contains

subroutine clComp(mode,ele,dir)
implicit none
   
   integer, intent(in) :: mode,ele,dir
   integer :: bblock,tblock
   logical :: flag
   real(nr) :: g11, g33, g13,coef,ninv
   real(nr) :: sumA,tsumA
   real(nr) :: dinp


   ! Define aerofoil blocks
   select case(ele)
   case(1)
   bblock = 1; tblock = 6
   case(2)
   bblock = 3; tblock = 8
   end select

   ! Find parallel grid position
   ip=mod(myid-mo(mb),npc(mb,1))
   jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
   kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

   ! Find master of the block
   mp = mo(mb)

   ! Find top or bottom
   if (mb==bblock) then
      coef=-one
   elseif (mb==tblock) then
      coef=one
   end if

   ! Initialise values
   cl(ele)=0;sumA=0;tsumA=0;flag=.false.
   g11=0;g33=g11;g13=g33

   if ((mb==bblock).AND.(jp==npc(mb,2)-1)) then
   j=ijk(1,2); flag=.true.
   elseif ((mb==tblock).AND.(jp==0)) then
   j=0; flag=.true.
   end if

   if (flag) then
   select case(mode)
   case(0) ! Compute Area
   do k=0,ijk(2,2)
   do i=0,ijk(3,2); l=indx3(j,k,i,2)
      g11 = qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
      g33 = qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
      g13 = qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
      dA(l) = sqrt(g11*g33-g13*g13)
      if ((ip==0).and.(i==0)) then
        dA(l) = dA(l)*half
      elseif ((ip==npc(mb,1)-1).and.(i==ijk(3,2))) then
        dA(l) = dA(l)*half
      end if
      if ((kp==0).and.(k==0)) then
        dA(l) = dA(l)*half
      elseif ((kp==npc(mb,3)-1).and.(k==ijk(2,2))) then
        dA(l) = dA(l)*half
      end if
   end do
   end do
   sumA = sum(dA)
   if (myid==mp) then
      do m = 1, (npc(mb,1)*npc(mb,3)-1)
      CALL MPI_RECV(tsumA,1,MPI_REAL8,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,ista,ierr)
      sumA=sumA+tsumA
      end do
      write(*,*) 'block',mb,'A =',sumA
   else
      CALL MPI_SEND(sumA,1,MPI_REAL8,mp,10,MPI_COMM_WORLD,ierr)
   end if
   case(1) ! Compute Cl
   ! Compute Dynamic pressure
   dinp=two/(umf(1)*umf(1)+umf(2)*umf(2)+umf(3)*umf(3))
   do k=0,ijk(2,2)
   do i=0,ijk(3,2); l=indx3(j,k,i,2)
      ninv = 1.0_nr/sqrt(etm(l,1)*etm(l,1)+etm(l,2)*etm(l,2)+etm(l,3)*etm(l,3))
      sumA = sumA + dinp*(-p(l))*coef*etm(l,dir)*ninv*dA(l)
   end do
   end do
   if (myid==mp) then
      do m = 1, (npc(mb,1)*npc(mb,3)-1)
      CALL MPI_RECV(tsumA,1,MPI_REAL8,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,ista,ierr)
      sumA=sumA+tsumA
      end do
      tsumA=sumA
   else
      CALL MPI_SEND(sumA,1,MPI_REAL8,mp,10,MPI_COMM_WORLD,ierr)
   end if
   end select
   end if

   CALL MPI_REDUCE(tsumA,cl(ele),1,MPI_REAL8,MPI_SUM,0,icom,ierr)
   
   if (myid==0) then
     if (mode==2) then
     close(9+ele)
     else
       if (n==0) then
          open(9+ele,file='misc/cl'//achar(ele+48)//'.dat')
          write(9+ele,*) 'variables=t,Cl'
       end if
          write(9+ele,'(2f8.5)') timo,cl(ele)
     end if
   end if

   
end subroutine clComp

end module rpt
