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

   CALL MPI_ALLREDUCE(tsumA,cl(ele,dir),1,MPI_REAL8,MPI_SUM,icom,ierr)

   
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

subroutine post

use problemcase, only: finalout,ngridv
   
!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE

 if(dt==0.or.nout==2) then
    ! rpt-File is not closed so it can still be saved up to the last record
    !close(0,status='delete')
    if (myid==0) then
    write(*,*) "Overflow."
    end if
 ndata=ndati
 end if
 !else
 if(myid==0) then
    write(*,'("Simulation time was ",f6.2," hours")') wtime/(3600_nr*npro)
    write(*,*) "Writing Output files..."
 end if
    if (ndatp==1) then
    call finalout
    end if
 if(myid==mo(mb)) then
    open(9,file=coutput); close(9,status='delete')
 end if
    call MPI_BARRIER(icom,ierr)
    open(9,file=coutput,access='stream',shared)
    lh=0
 if(myid==mo(mb)) then
    write(9,pos=4*lh+1) '#!TDV112'; lh=lh+2
    write(9,pos=4*lh+1) 1; lh=lh+1 ! Header Section
    write(9,pos=4*lh+1) 0; lh=lh+1 ! File Type
    cinput='title'; call strio(9,lh,cinput) ! File Title
    write(9,pos=4*lh+1) int4(nwrec); lh=lh+1 ! Number of Variables
    cinput='x'; call strio(9,lh,cinput)
    cinput='y'; call strio(9,lh,cinput)
    cinput='z'; call strio(9,lh,cinput)
    if (ngridv==1) then
       cinput='xix'; call strio(9,lh,cinput)
       cinput='xiy'; call strio(9,lh,cinput)
       cinput='xiz'; call strio(9,lh,cinput)
       cinput='etax'; call strio(9,lh,cinput)
       cinput='etay'; call strio(9,lh,cinput)
       cinput='etaz'; call strio(9,lh,cinput)
       cinput='zetax'; call strio(9,lh,cinput)
       cinput='zetay'; call strio(9,lh,cinput)
       cinput='zetaz'; call strio(9,lh,cinput)
    end if
 do n=0,ndata
    no(2)=n/100; no(1)=mod(n,100)/10; no(0)=mod(n,10); cno=achar(no+48)
    cinput='r'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='u'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='v'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='w'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='p'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
 end do
 if (ndatp==1) then
    cinput='r'; call strio(9,lh,cinput)
    cinput='u'; call strio(9,lh,cinput)
    cinput='v'; call strio(9,lh,cinput)
    cinput='w'; call strio(9,lh,cinput)
    cinput='p'; call strio(9,lh,cinput)
 end if
    write(9,pos=4*lh+1) 299.0; lh=lh+1 ! Zone Marker
    cinput=czone; call strio(9,lh,cinput)
    write(9,pos=4*lh+1) -1; lh=lh+1 ! Parent Zone
    write(9,pos=4*lh+1) -2; lh=lh+1 ! Strand ID
    write(9,pos=4*lh+1) dble(0.0); lh=lh+2 ! Solution Time (Double)
    write(9,pos=4*lh+1) -1; lh=lh+1 ! (Not used. Set to -1.)
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Zone Type
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Specify Var Location
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Raw Local 1-to-1 Face Neighbours Suppliled
    write(9,pos=4*lh+1) 0; lh=lh+1 ! Number of Miscellaneous Face Neighbour Connections
    write(9,pos=4*lh+1) int4(lximb(mb)+1); lh=lh+1 ! IMax
    write(9,pos=4*lh+1) int4(letmb(mb)+1); lh=lh+1 ! JMax
    write(9,pos=4*lh+1) int4(lzemb(mb)+1); lh=lh+1 ! KMax
    write(9,pos=4*lh+1) 0; lh=lh+1 ! No Auxillary Data Pairs
    write(9,pos=4*lh+1) 357.0; lh=lh+1 ! End of Header Marker
    write(9,pos=4*lh+1) 299.0; lh=lh+1 ! Zone Marker
 do n=1,nwrec
    write(9,pos=4*lh+1) 1; lh=lh+1 ! 1 = Float / 2 = Double
 end do
    write(9,pos=4*lh+1) 0; lh=lh+1 ! No Passive Variables
    write(9,pos=4*lh+1) 0; lh=lh+1 ! No Variable Sharing
    write(9,pos=4*lh+1) -1; lh=lh+1 ! Zero Based Zone Number to Share
 do n=1,nwrec
    lh=lh+2 ! Minimum Value (Double) of Variables (to be filled)
    lh=lh+2 ! Maximum Value (Double) of Variables (to be filled)
 end do
    lhmb(mb)=lh
 end if

 ! rpt- COMPUTE MAX AND MIN
 do mm=0,mbk
    call MPI_BCAST(lhmb(mm),1,MPI_INTEGER,mo(mm),icom,ierr)
 end do
    ns=1; ne=nwrec; allocate(varmin(ns:ne),varmax(ns:ne))
    lp=lpos(myid)+lhmb(mb)
 do n=ns,ne; lq=(n-1)*ltomb
    read(0,rec=n) varr(:)
 do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
    write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
 end do; end do
    varmin(n)=minval(varr(:)); varmax(n)=maxval(varr(:))
 end do
    close(0,status='delete')
 do n=ns,ne
    res=varmin(n); call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MIN,icom,ierr); varmin(n)=fctr
    res=varmax(n); call MPI_ALLREDUCE(res,fctr,1,MPI_REAL8,MPI_MAX,icom,ierr); varmax(n)=fctr
 end do
 if(myid==mo(mb)) then
    l=0; lq=4*(nwrec)
 do n=ns,ne
    write(9,pos=4*(lp-lq+l)+1) dble(varmin(n)); l=l+2 ! 8-Bytes "Stream"
    write(9,pos=4*(lp-lq+l)+1) dble(varmax(n)); l=l+2 ! 8-Bytes "Stream"
 end do
 end if
    close(9)
 !end if

end subroutine post

end module rpt
