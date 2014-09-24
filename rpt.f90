!**********************
!***** RPT MODULE *****
!**********************

module  rpt

use mainvar3d
use subroutineso

contains

subroutine clComp(mode,ele)
   
   integer, intent(in) :: mode,ele
   integer :: bblock,tblock
   real(nr) :: g11, g33, g13
   real(nr) :: sumA

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

   if ((mb==bblock.and.jp==npc(mb,2)-1).or.(mb==tblock.and.jp==0)) then
   select case(mode)
   case(0) ! Compute Area
   if(jp==0) then 
      i==0
   elseif (jp==npc(mb,2)-1) then
      i==ijk(1,2)
   end if
   do k=0,ijk(3,2)
   do j=0,ijk(2,2)
      g11 = qo(l,1)*qo(l,1)+qa(l,1)qa(l,1)+de(l,1)*de(l,1)
      g33 = qo(l,3)*qo(l,3)+qa(l,3)qa(l,3)+de(l,3)*de(l,3)
      g13 = qo(l,1)*qo(l,3)+qa(l,1)qa(l,3)+de(l,1)*de(l,3)
      dA(l) = sqrt(g11*g33-g13*g13)
      if ((ip==0).or.(ip==npc(mb,1)-1)) then
        dA(l) = dA(l)*half
        if ((kp==0).or.(ip==npc(mb,3)-1)) then
           dA(l) = dA(l)*half
        end if
      end if
   end do
   end do
   sumA = sum(dA)
   CALL MPI_REDUCE(sumA,sumA,1,MPI_REAL8,MPI_SUM,mp,icom,ierr)
   if (myid==mp) then
     write(*,*) 'block',mp,'A =',sumA
   end if
   end select
   end if

end subroutine clComp

end module rpt
