!**********************
!***** RPT MODULE *****
! Mode 0: Compute Area
! Mode 1: Compute Force coefficient
! Mode 2: Close file
!**********************

module  rpt

use mainvar3d
use subroutineso
use mpi

contains

subroutine clComp(mode,ele,dir)

use problemcase, only: span,delt1,delt2
implicit none
   
   integer, intent(in) :: mode,ele,dir
   integer :: bblock,tblock
   logical :: flag
   real(nr) :: g11, g33, g13,coef,normal
   real(nr) :: sumA,tsumA
   real(nr) :: dinp
   integer :: ddelt1,ddelt2
   integer :: wunit
   character :: char1,char2,char3
   character, dimension (4) :: cdelt1,cdelt2

   wunit=10*(dir)+(ele-1)
   ddelt1=delt1*1800.0_nr/pi
   ddelt2=delt2*1800.0_nr/pi


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



   ! Find top or bottom
   if (mb==bblock) then
      coef=-one
   elseif (mb==tblock) then
      coef=one
   end if

   ! Initialise values
   cl(ele,dir)=0;sumA=0;tsumA=0;flag=.false.
   g11=0;g33=g11;g13=g33

   if ((mb==bblock).AND.(jp==npc(mb,2)-1)) then
   ! Find master of the block
   mp = mo(mb) + npc(mb,1)*(npc(mb,2)-1)
   j=ijk(1,2); flag=.true.
   elseif ((mb==tblock).AND.(jp==0)) then
   j=0; flag=.true.
   ! Find master of the block
   mp = mo(mb)
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
   dinp=two/(amachoo*amachoo*span)
   do k=0,ijk(2,2)
   do i=0,ijk(3,2); l=indx3(j,k,i,2)
      normal = etm(l,dir)*coef/sqrt(etm(l,1)*etm(l,1)+etm(l,2)*etm(l,2)+etm(l,3)*etm(l,3))
      sumA = sumA + (p(l))*normal*dA(l)
   end do
   end do
   sumA = sumA * dinp
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

   CALL MPI_REDUCE(tsumA,cl(ele,dir),1,MPI_REAL8,MPI_SUM,0,icom,ierr)
   
   if (myid==0) then
     if (mode==2) then
     close(wunit)
     else
       if (n==0) then
          cdelt1(1)=achar((ddelt1/100)+48); 
          cdelt1(2)=achar((mod(ddelt1,100)/10)+48); 
          cdelt1(3)= '.'
          cdelt1(4)=achar((mod(ddelt1,10))+48); 
          cdelt2(1)=achar((int(ddelt2)/100)+48); 
          cdelt2(2)=achar((mod(ddelt2,100)/10)+48); 
          cdelt2(3)= '.'
          cdelt2(4)=achar((mod(ddelt2,10))+48); 
          select case(dir)
          case(2)
          open(wunit,file='misc/cl'//achar(ele+48)//'.dat')
          write(wunit,"('VARIABLES= t, C<sub>l',I1,'</sub>')") ele 
          write(wunit,"('ZONE T= Cl',I1,'A',4a,'_',4a)") ele,cdelt1,cdelt2
          case(1)
          open(wunit,file='misc/cd'//achar(ele+48)//'.dat')
          write(wunit,"('VARIABLES= t, C<sub>d',I1,'</sub>')") ele 
          write(wunit,"('ZONE T= Cd',I1,'A',4a,'_',4a)") ele,cdelt1,cdelt2
          end select
       end if
          write(wunit,'(f8.5,"   ",f10.5)') timo,cl(ele,dir)
     end if
   end if

   
end subroutine clComp

end module rpt
