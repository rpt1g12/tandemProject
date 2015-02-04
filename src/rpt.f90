!**********************
!***** RPT MODULE *****
!**********************

module  rpt

use mainvar3d
use subroutineso
use mpi

contains

subroutine cpComp(ele)
use problemcase, only: span,delt1,delt2,wlew,wlea
implicit none

   integer, intent(in) :: ele
   logical :: flag
   integer :: nwlew
   integer :: bblock,tblock
   integer :: loc
   integer :: olxi
   integer :: dest,omp
   real(nr) :: wlew4,zcp,dinp
   real(nr),dimension(0:lxi-1) :: mycp
   real(nr),dimension(:,:),allocatable :: cp,tcp
   character,dimension(3) :: cloc
   character(4) :: cblock

   nwlew=4*int(span/wlew)
   wlew4=wlew*quarter
   flag=.false.

   ! Find parallel grid position
   ip=mod(myid-mo(mb),npc(mb,1))
   jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
   kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

   ! Define aerofoil blocks
   select case(ele)
   case(1)
   bblock = 6; tblock = 11
   cblock='fore'
   case(2)
   bblock = 8; tblock = 13
   cblock='aftr'
   end select


   if ((mb==bblock).AND.(jp==npc(mb,2)-1)) then
   ! Find master of the block
   mp = mo(mb) + npc(mb,1)*(npc(mb,2)-1)
   omp = mo(tblock)
   j=ijk(1,2);flag=.true. 
   elseif ((mb==tblock).AND.(jp==0)) then
   ! Find master of the block
   mp = mo(mb)
   omp = mo(bblock) + npc(bblock,1)*(npc(bblock,2)-1)
   j=0;flag=.true.
   end if

   if (flag) then
   if (myid==mp) then
      allocate(cp(0:lximb(mb)-1,0:nwlew))
      if (j==0) then
      allocate(tcp(0:lximb(mb)+lximb(bblock)-1,0:nwlew))
      end if
   end if
   do loc = 0, nwlew
      zcp=loc*wlew4
      k=int(zcp/(span/lzemb(mb)))
      l=indx3(0,j,k,1)
      mycp(:)=p(l:l+lxi)
      mycp(:)=two*(mycp(:)-poo)/(amachoo**2)
      if (myid==mp) then
         cp(0:lxi-1,loc)=mycp(:)
         do m = 0, npc(mb,1)-2
         dest=myid+m+1
         olxi=lxim(dest)
         CALL MPI_RECV(cp(lxi+m*olxi+1,loc),olxi,MPI_DOUBLE,dest,dest,MPI_COMM_WORLD,ista,ierr)
         end do
      else
         CALL MPI_SEND(mycp,lxi,MPI_DOUBLE,mp,myid,MPI_COMM_WORLD,ierr)
      end if

      if (myid==mp) then
         if (j==0) then
            CALL MPI_RECV(tcp(lximb(mb),loc),lximb(bblock),MPI_DOUBLE,omp,omp,MPI_COMM_WORLD,ista,ierr)
            cloc(1)=achar((loc/100)+48)
            cloc(2)=achar((mod(loc,100)/10)+48)
            cloc(3)=achar(mod(loc,10)+48)
            open(3,file='loc/'//cblock//cloc(1)//cloc(2)//cloc(3)//'.dat')
            write(3,"('VARIABLES= x,cp')") 
            write(3,"('ZONE T= ',a4,'loc',3a)") cblock,cloc
            tcp(0:lximb(mb)-1,loc)=cp(:,loc)
         do i = 0, lximb(mb)-1
           write(3,'(f10.5,"   ",f10.5)') tpwle(i,loc,1),tcp(i,loc)
         end do
         do i = lximb(mb)+lximb(bblock)-1, lximb(mb), -1
           write(3,'(f10.5,"   ",f10.5)') tpwle(i,loc,1),tcp(i,loc)
         end do
         close(3)
         else
            CALL MPI_SEND(cp(0,loc),lximb(mb),MPI_DOUBLE,omp,mp,MPI_COMM_WORLD,ierr)
         end if
      end if
   end do
   end if

end subroutine cpComp
!===============================================
subroutine spanLoc(ele)
use problemcase, only: span,delt1,delt2,wlew,wlea
implicit none

   integer, intent(in) :: ele
   logical :: flag
   integer :: nwlew
   integer :: bblock,tblock
   integer :: loc
   integer :: olxi
   integer :: dest,omp
   real(nr) :: wlew4,zcp
   real(nr),dimension(0:lxi-1) :: mypwle
   real(nr),dimension(:,:,:),allocatable :: pwle
   character,dimension(3) :: cloc

   nwlew=4*int(span/wlew)
   wlew4=wlew*quarter
   flag=.false.

   ! Find parallel grid position
   ip=mod(myid-mo(mb),npc(mb,1))
   jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
   kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))

   ! Define aerofoil blocks
   select case(ele)
   case(1)
   bblock = 6; tblock = 11
   case(2)
   bblock = 8; tblock = 13
   end select


   if ((mb==bblock).AND.(jp==npc(mb,2)-1)) then
   ! Find master of the block
   mp = mo(mb) + npc(mb,1)*(npc(mb,2)-1)
   omp = mo(tblock)
   j=ijk(1,2);flag=.true. 
   elseif ((mb==tblock).AND.(jp==0)) then
   ! Find master of the block
   mp = mo(mb)
   omp = mo(bblock) + npc(bblock,1)*(npc(bblock,2)-1)
   j=0;flag=.true.
   end if

   if (flag) then
   if (myid==mp) then
      allocate(pwle(0:lximb(mb)-1,0:nwlew,3))
      if (j==0) then
      allocate(tpwle(0:lximb(mb)+lximb(bblock)-1,0:nwlew,3))
      end if
   end if
   do loc = 0, nwlew
      zcp=loc*wlew4
      k=int(zcp/(span/lzemb(mb)))
      do nn=1,3
      l=indx3(0,j,k,1)
      mypwle(:)=ss(l:l+lxi,nn)
      if (myid==mp) then
         pwle(0:lxi-1,loc,nn)=mypwle(:)
         do m = 0, npc(mb,1)-2
         dest=myid+m+1
         olxi=lxim(dest)
         CALL MPI_RECV(pwle(lxi+m*olxi+1,loc,nn),olxi,MPI_DOUBLE,dest,dest,MPI_COMM_WORLD,ista,ierr)
         end do
      else
         CALL MPI_SEND(mypwle,lxi,MPI_DOUBLE,mp,myid,MPI_COMM_WORLD,ierr)
      end if
      if (myid==mp) then
         if (j==0) then
            tpwle(0:lximb(mb)-1,loc,nn)=pwle(:,loc,nn)
            CALL MPI_RECV(tpwle(lximb(mb),loc,nn),lximb(bblock),MPI_DOUBLE,omp,omp,MPI_COMM_WORLD,ista,ierr)
         else
            CALL MPI_SEND(pwle(0,loc,nn),lximb(mb),MPI_DOUBLE,omp,mp,MPI_COMM_WORLD,ierr)
         end if
      end if
      end do
   end do
   end if

end subroutine spanLoc
!===============================================
subroutine clComp(mode,ele,dir)
! Mode 0: Compute Area
! Mode 1: Compute Force coefficient
! Mode 2: Close file

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
   character, dimension (4) :: cdelt1,cdelt2

   wunit=10*(dir)+(ele-1)
   ddelt1=delt1*1800.0_nr/pi
   ddelt2=delt2*1800.0_nr/pi


   ! Define aerofoil blocks
   select case(ele)
   case(1)
   bblock = 6; tblock = 11
   case(2)
   bblock = 8; tblock = 13
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

!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE
subroutine post

use problemcase, only: finalout,ngridv

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
    if (nvarout==7) then
    cinput='dux'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='duy'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='duz'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='dvx'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='dvy'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='dvz'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='dwx'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='dwy'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    cinput='dwz'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
    end if
 end do
 if (ndatp==1) then
    cinput='r'; call strio(9,lh,cinput)
    cinput='u'; call strio(9,lh,cinput)
    cinput='v'; call strio(9,lh,cinput)
    cinput='w'; call strio(9,lh,cinput)
    cinput='p'; call strio(9,lh,cinput)
    if (nvarout==7) then
    cinput='dux'; call strio(9,lh,cinput)
    cinput='duy'; call strio(9,lh,cinput)
    cinput='duz'; call strio(9,lh,cinput)
    cinput='dvx'; call strio(9,lh,cinput)
    cinput='dvy'; call strio(9,lh,cinput)
    cinput='dvz'; call strio(9,lh,cinput)
    cinput='dwx'; call strio(9,lh,cinput)
    cinput='dwy'; call strio(9,lh,cinput)
    cinput='dwz'; call strio(9,lh,cinput)
    end if
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

! ====SET UP FORCING PARAMETERS
subroutine forceup

use problemcase, only: span

ra0=-log(0.0001_nr)/(rfor**2); ra1=2*pi/span
ra2=rfor**2
ll=-1
do l = 0, lmx
  rr(l,1)=(ss(l,1)-xfor)**2+(ss(l,2)-yfor)**2
  if (rr(l,1)-ra2<0) then
     ll=ll+1; de(ll,5)=l+sml
  end if
end do
lfor=ll
if (lfor.ne.-1) then
   allocate(lcfor(0:lfor),xafor(0:lfor,3),yafor(0:lfor,3),bfor(0:lfor,3))
do ll = 0, lfor; l=de(ll,5); lcfor(ll)=l
   bfor(ll,1)=cos(ra1*ss(l,3))
   bfor(ll,2)=cos(3*ra1*ss(l,3))
   bfor(ll,3)=cos(4*ra1*ss(l,3))
   ra2=rr(l,1)
   ra3=(qo(l,2)-qo(l,1))
   xafor(ll,1)=half*exp(-ra0*ra2)*bfor(ll,1)*ra3
   xafor(ll,2)=half*exp(-ra0*ra2)*bfor(ll,2)*ra3
   xafor(ll,3)=half*exp(-ra0*ra2)*bfor(ll,3)*ra3
   ra3=(qa(l,2)-qa(l,1))
   yafor(ll,1)=half*exp(-ra0*ra2)*bfor(ll,1)*ra3
   yafor(ll,2)=half*exp(-ra0*ra2)*bfor(ll,2)*ra3
   yafor(ll,3)=half*exp(-ra0*ra2)*bfor(ll,3)*ra3
end do
end if
   
end subroutine forceup

! ====FORCE IMPLEMENTATION
subroutine forcego

if ((timo+dtk>tsfor).and.(timo+dtk<tefor)) then
  ra0=(timo+dtk-tsfor)
  ra1=amfor*cos(ra0*48.76_nr*amachoo)/3
  ra2=amfor*cos(ra0*53.6_nr*amachoo)/3
  ra3=amfor*cos(ra0*53.6_nr*amachoo)/3
  do ll = 0, lfor; l=lcfor(ll)
    de(l,2)=de(l,2)+ra1*xafor(ll,1)+ra2*xafor(ll,2)+ra3*xafor(ll,3)
    de(l,3)=de(l,3)+ra1*yafor(ll,1)+ra2*yafor(ll,2)+ra3*yafor(ll,3)
  end do
end if
end subroutine forcego

! ====RECORD DATA FOR POST-PROCESSING
subroutine postDat
 call MPI_BARRIER(icom,ierr)
 open(3,file=cpostdat,access='stream',shared)
 lp=lpos(myid)
 do n=1,nwrec
    lq=(n-1)*ltomb
    read(0,rec=n) varr(:)
    do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
       write(3,pos=nr*(lp+lq+lio(j,k))+1) varr(l:l+lxi)
    end do; end do
 end do
 close(3)
end subroutine postDat

!=====COMPUTE CELL AREA OVER AEROFOILS
 subroutine wallArea
 implicit none
    integer :: bblock1,tblock1
    integer :: bblock2,tblock2
    real(nr) :: g11, g33, g13,coef
 
    bblock1 = 6; tblock1 = 11
    bblock2 = 8; tblock2 = 13
 
    ! Find parallel grid position
    ip=mod(myid-mo(mb),npc(mb,1))
    jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))
 
    ! Initialise values
    wflag=.false.
    g11=0;g33=g11;g13=g33
 
    if ((mb==bblock1).AND.(jp==npc(mb,2)-1)) then
    j=ijk(1,2); wflag=.true.
    elseif ((mb==tblock1).AND.(jp==0)) then
    j=0; wflag=.true.
    end if
    if ((mb==bblock2).AND.(jp==npc(mb,2)-1)) then
    j=ijk(1,2); wflag=.true.
    elseif ((mb==tblock2).AND.(jp==0)) then
    j=0; wflag=.true.
    end if
 
    ll=-1; lcwall=(nbsize(2)-1)
    if (wflag) then
    allocate(lwall(0:lcwall),area(0:lcwall))
    do k=0,ijk(2,2)
    do i=0,ijk(3,2); l=indx3(j,k,i,2)
       g11 = qo(l,1)*qo(l,1)+qa(l,1)*qa(l,1)+de(l,1)*de(l,1)
       g33 = qo(l,3)*qo(l,3)+qa(l,3)*qa(l,3)+de(l,3)*de(l,3)
       g13 = qo(l,1)*qo(l,3)+qa(l,1)*qa(l,3)+de(l,1)*de(l,3)
       ll=ll+1; lwall(ll)=l+sml
       area(ll)=sqrt(g11*g33-g13*g13)
       if ((ip==0).and.(i==0)) then
         area(ll) = area(ll)*half
       elseif ((ip==npc(mb,1)-1).and.(i==ijk(3,2))) then
         area(ll) = area(ll)*half
       end if
       if ((kp==0).and.(k==0)) then
         area(ll) = area(ll)*half
       elseif ((kp==npc(mb,3)-1).and.(k==ijk(2,2))) then
         area(ll) = area(ll)*half
       end if
    end do
    end do
    end if
 
 end subroutine wallArea

!=====COMPUTE WALL NORMAL VECTOR OVER AEROFOILS
 subroutine walldir
 implicit none
    
    integer :: bblock1,tblock1
    integer :: bblock2,tblock2
    real(nr) :: coef,tmp
 
    bblock1 = 6; tblock1 = 11
    bblock2 = 8; tblock2 = 13
 
    if (wflag) then
    ! Find top or bottom
    if ((mb==bblock1).or.(mb==bblock2)) then
       coef=-one
    elseif ((mb==tblock1).or.(mb==tblock2)) then
       coef=one
    end if
    allocate(wnor(0:lcwall,3))
    do ll = 0, lcwall; l=lwall(ll)
       tmp=coef/sqrt(etm(l,1)*etm(l,1)+etm(l,2)*etm(l,2)+etm(l,3)*etm(l,3))
       do m = 1, 3
       wnor(ll,m)=etm(l,m)*tmp
       end do
    end do
    end if
 
 end subroutine walldir

!=====COMPUTE LIFT COEFFICIENT OVER AEROFOILS
 subroutine clpost(ele,dir,nvar,post)
 
 use problemcase, only: span,delt1,delt2
 use subroutines3d, only: mpigo,deriv
 implicit none
    
    integer, intent(in) :: ele,nvar,dir
    integer :: bblock,tblock
    real(nr) :: dynp,clp,clv,tcl
    logical :: flag,post
 
    clp=0;clv=0;flag=.false.;tcl=0;
 
    ! Define aerofoil blocks
    select case(ele)
    case(1)
    bblock = 6; tblock = 11
    case(2)
    bblock = 8; tblock = 13
    end select
 
    if (mb==bblock) then
    ! Find master of the block
    mp = mo(mb) + npc(mb,1)*(npc(mb,2)-1)
    flag=.true.
    elseif (mb==tblock) then
    ! Find master of the block
    mp = mo(mb)
    flag=.true.
    end if
 
    if (wflag.and.flag) then
    if (post) then
       ! READ VARIABLES
       nread=nrec+(totVar*nvar)
       do nn = 1, 5
        nread=nread+1 
        call postread(nread)
        qo(:,nn)=varr(:)
       end do
       p(:)=qo(:,5)
    end if
    ! Compute Dynamic pressure
    dynp=two/(amachoo*amachoo*span)
    do ll = 0, lcwall; l=lwall(ll)
      clp=clp+(p(l)*wnor(ll,dir)*area(ll))
    end do
    if (nviscous==1) then
     if (post) then
     de(:,1)=1/qo(:,1)
     de(:,2)=qo(:,2)
     de(:,3)=qo(:,3)
     de(:,4)=qo(:,4)
     de(:,5)=gam*qo(:,5)*de(:,1)
     ss(:,1)=srefp1dre*de(:,5)**1.5_nr/(de(:,5)+srefoo)
     de(:,1)=ss(:,1)
 
     rr(:,1)=de(:,2)
     m=2; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     txx(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
     hzz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
     tzx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)
 
     rr(:,1)=de(:,3)
     m=3; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     txy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
     tyy(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
     hxx(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)
 
     rr(:,1)=de(:,4)
     m=4; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     hyy(:)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
     tyz(:)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
     tzz(:)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)
 
     fctr=2.0_nr/3
     rr(:,1)=de(:,1)*yaco(:)
     de(:,5)=fctr*(txx(:)+tyy(:)+tzz(:))
 
     txx(:)=yaco(:)*(2*txx(:)-de(:,5))
     tyy(:)=yaco(:)*(2*tyy(:)-de(:,5))
     tzz(:)=yaco(:)*(2*tzz(:)-de(:,5))
     txy(:)=yaco(:)*(txy(:)+hzz(:))
     tyz(:)=yaco(:)*(tyz(:)+hxx(:))
     tzx(:)=yaco(:)*(tzx(:)+hyy(:))
     end if
 
       if(.not.allocated(tw)) allocate(tw(0:lcwall,3))
       do ll = 0, lcwall; l=lwall(ll)
         tw(ll,1)=(txx(l)*wnor(ll,1)+txy(l)*wnor(ll,2)+tzx(l)*wnor(ll,3))/reoo
         tw(ll,2)=(txy(l)*wnor(ll,1)+tyy(l)*wnor(ll,2)+tyz(l)*wnor(ll,3))/reoo
         tw(ll,3)=(tzx(l)*wnor(ll,1)+tyz(l)*wnor(ll,2)+tzz(l)*wnor(ll,3))/reoo
         clv=clv+tw(ll,dir)*area(ll)
       end do
    end if
    tcl=(clp+clv)*dynp
    if (myid==mp) then
       do m = 1, (npc(mb,1)*npc(mb,3)-1)
       CALL MPI_RECV(clp,1,MPI_REAL8,MPI_ANY_SOURCE,10,MPI_COMM_WORLD,ista,ierr)
       tcl=tcl+clp
       end do
    else
       CALL MPI_SEND(tcl,1,MPI_REAL8,mp,10,MPI_COMM_WORLD,ierr)
    end if
    end if
 
    CALL MPI_ALLREDUCE(tcl,cl(ele,dir),1,MPI_REAL8,MPI_SUM,icom,ierr)
 
 end subroutine clpost

! ====READ DATA FOR POST-PROCESSING
subroutine postread(nread)
implicit none
integer, intent (in) :: nread
 n=nread
 lp=lpos(myid)
    lq=(n-1)*ltomb
    do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
       read(3,pos=nr*(lp+lq+lio(j,k))+1) varr(l:l+lxi)
    end do; end do
end subroutine postread

! ====WRITE DATA FOR POST-PROCESSING
subroutine postwrite(nread)
implicit none
integer, intent (in) :: nread
 call MPI_BARRIER(icom,ierr)
 n=nread
 lp=lpos(myid)
    lq=(n-1)*ltomb
    do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
       write(3,pos=nr*(lp+lq+lio(j,k))+1) varr(l:l+lxi)
    end do; end do
end subroutine postwrite

! ====WRITE DATA FOR POST-PROCESSING
subroutine datawrite
implicit none
integer :: n
 call MPI_BARRIER(icom,ierr)
    inquire(iolength=lh) varr
    open(0,file=cdata,access='direct',recl=lh)
    do n = 1, nwrec
       call postread(n)
       write(0,rec=nwrec) varr(:)
    end do
end subroutine datawrite
!===============================================
end module rpt
