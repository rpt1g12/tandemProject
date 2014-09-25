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
   logical :: flag
   real(nr) :: g11, g33, g13
   real(nr) :: sumA

   ! Define aerofoil blocks
   select case(ele)
   case(1)
   bblock = 2; tblock = 7
   case(2)
   bblock = 3; tblock = 8
   end select

   ! Find parallel grid position
   ip=mod(myid-mo(mb),npc(mb,1))
   jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
   kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

   ! Find master of the block
   mp = mo(mb)

   ! Initialise area
   sumA=0;g11=0;flag=.false.

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
      if ((ip==0).or.(ip==npc(mb,1)-1)) then
        if (i==0.or.i==ijk(3,2)) then
        dA(l) = dA(l)*half
        end if
        if ((kp==0).or.(ip==npc(mb,3)-1)) then
        if (k==0.or.k==ijk(3,2)) then
        dA(l) = dA(l)*half
        end if
        end if
      end if
   end do
   end do
   sumA = sum(dA)
   if (myid==mp) then
     write(*,*) 'block',mb,'A =',sumA
   end if
   end select
   end if

end subroutine clComp

end module rpt
