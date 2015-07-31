!**********************
!***** RPT MODULE *****
!**********************

module  rpt

use mainvar3d
use subroutineso
use mpi

contains

!====================================================================================
!===== POST-PROCESSING & GENERATING PLOT3D DATA FILES
!====================================================================================
 subroutine plot3dgrid()
 integer :: n

      if (myid==0) then
        open(9,file='out/grid.xyz'); close(9,status='delete')
      end if
      CALL MPI_BARRIER(icom,ierr)
      open (unit=9, file='out/grid.xyz', access='stream',shared)
      lh=0
      if (myid==0) then
       write(9,pos=4*lh+1) mbk+1; lh=lh+1 ! Number of zones
       do mm = 0, mbk
          write(9,pos=4*lh+1) int4(lximb(mm)+1); lh=lh+1 ! IMax
          write(9,pos=4*lh+1) int4(letmb(mm)+1); lh=lh+1 ! JMax
          write(9,pos=4*lh+1) int4(lzemb(mm)+1); lh=lh+1 ! KMax
       end do
       lhmb(mb)=lh
       do mm = 0, mbk-1
          lhmb(mm+1)=lhmb(mm)+3*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
       end do
      end if
       call MPI_BCAST(lhmb,mbk+1,MPI_INTEGER,0,icom,ierr)
       lp=lpos(myid)+lhmb(mb)
       ns=1; ne=3
       do n=ns,ne; lq=(n-1)*ltomb
          varr(:)=ss(:,n)
       do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
          write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
       end do; end do
       end do
       close(9)
      CALL MPI_BARRIER(icom,ierr)
        if (myid==0) then
        write(*,*) 'Grid written!'
        end if
    
 end subroutine plot3dgrid

 subroutine plot3dsolution(ctime)
 character(8), intent (in) :: ctime
 integer :: n

       if (myid==0) then
         open(9,file='out/solT'//ctime//'.q'); close(9,status='delete')
       end if
       CALL MPI_BARRIER(icom,ierr)
       open (unit=9, file='out/solT'//ctime//'.q', access='stream',shared)
       lh=0
       if (myid==0) then
        write(9,pos=4*lh+1) mbk+1; lh=lh+1 ! Number of zones
        do mm = 0, mbk
           write(9,pos=4*lh+1) int4(lximb(mm)+1); lh=lh+1 ! IMax
           write(9,pos=4*lh+1) int4(letmb(mm)+1); lh=lh+1 ! JMax
           write(9,pos=4*lh+1) int4(lzemb(mm)+1); lh=lh+1 ! KMax
        end do
        lhmb(mb)=lh
        do mm = 0, mbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do
       end if
        call MPI_BCAST(lhmb,mbk+1,MPI_INTEGER,0,icom,ierr)
        if (myid==mo(mb)) then
           lh=lhmb(mb)
           write(9,pos=4*lh+1) real(amachoo,kind=4); lh=lh+1 ! Mach Number
           write(9,pos=4*lh+1) real(aoa,kind=4); lh=lh+1  ! AoA
           write(9,pos=4*lh+1) real(reoo,kind=4); lh=lh+1 ! Reynolds Number
           write(9,pos=4*lh+1) real(timo,kind=4); lh=lh+1 ! Time
        end if
        lp=lpos(myid)+lhmb(mb)+4
        ns=1; ne=5
        
        do n=ns,ne; lq=(n-ns)*ltomb
           selectcase(n)
           case(1); varr(:)=qa(:,n)
           case(2,3,4); !de(:,n)=qa(:,1)*umf(n-1)
           varr(:)=((qa(:,n)/qa(:,1))+umf(n-1))
           case(5); varr(:)=p(:)
           end select
        do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
          write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
        end do; end do
        end do
        close(9)
        if (myid==0) then
           write(*,"('Solution written! T= ',8a)") ctime 
        end if
    
 end subroutine plot3dsolution

 subroutine plot3d(gflag,sflag,bflag)

 integer, intent(in) :: gflag,sflag,bflag
 integer :: n
 character(8) :: ctime
 
 write(ctime,"(f8.4)") timo
    

    if (gflag==1) then
       call plot3dgrid
    end if

    if (sflag==1) then
       call plot3dsolution(ctime)
    end if
          
 end subroutine plot3d

!====================================================================================
!===== POST-PROCESSING & GENERATING TECPLOT DATA FILE
!====================================================================================
 subroutine post(average)
 use problemcase, only: finalout,ngridv
 logical, intent(in) :: average
 
  if(myid==mo(mb)) then
     if (ispost) then
     open(9,file=ctecout); close(9,status='delete')
     else
     open(9,file=coutput); close(9,status='delete')
     end if
  end if
     call MPI_BARRIER(icom,ierr)
     if (ispost) then
     open(9,file=ctecout,access='stream',shared)
     else
     open(9,file=coutput,access='stream',shared)
     end if
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
     if (average) then
     cinput='r'; call strio(9,lh,cinput)
     cinput='u'; call strio(9,lh,cinput)
     cinput='v'; call strio(9,lh,cinput)
     cinput='w'; call strio(9,lh,cinput)
     cinput='p'; call strio(9,lh,cinput)
     end if
  !do n=0,ndata
  !   no(2)=n/100; no(1)=mod(n,100)/10; no(0)=mod(n,10); cno=achar(no+48)
  !   cinput='Q'//cno(2)//cno(1)//cno(0); call strio(9,lh,cinput)
  !end do
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
     if (ispost) then
     call postread(n)   
     else
     read(0,rec=n) varr(:)
     end if
  do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
     write(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
  end do; end do
     varmin(n)=minval(varr(:)); varmax(n)=maxval(varr(:))
  end do
     if (ispost) then
        close(8)
     else
        if (nto==0) then
        close(0,status='delete')
        else
        close(0)
        end if
     end if
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
 
 end subroutine post

!====================================================================================
! ====SET UP FORCING PARAMETERS
!====================================================================================
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
    allocate(lcfor(0:lfor),xafor(0:lfor,3),bfor(0:lfor,3))
 do ll = 0, lfor; l=de(ll,5); lcfor(ll)=l
    bfor(ll,1)=cos(ra1*ss(l,3))
    bfor(ll,2)=cos(3*ra1*ss(l,3))
    bfor(ll,3)=cos(4*ra1*ss(l,3))
    ra1=half*exp(-ra0*rr(l,1))/yaco(l)
    !rpt- Xeta-Xxi
    !ra3=(qo(l,2)-qo(l,1))
    !ra1=half*exp(-ra0*rr(l,1))*(qo(l,2)-qo(l,1))
    !rpt- Yeta-Yxi
    !ra3=(qa(l,2)-qa(l,1))
    !ra2=half*exp(-ra0*rr(l,1))*(qa(l,2)-qa(l,1))
    xafor(ll,1)=ra1*bfor(ll,1)
    xafor(ll,2)=ra1*bfor(ll,2)
    xafor(ll,3)=ra1*bfor(ll,3)
 end do
 end if
 end subroutine forceup

!====================================================================================
! ====FORCE IMPLEMENTATION
!====================================================================================
 subroutine forcego
 if ((timo+dtk>tsfor).and.(timo+dtk<tefor)) then
   ra0=(timo+dtk-tsfor)
   ra1=amfor*cos(ra0*48.76_nr*amachoo)/3
   ra2=amfor*cos(ra0*53.6_nr*amachoo)/3
   ra3=amfor*cos(ra0*53.6_nr*amachoo)/3
   do ll = 0, lfor; l=lcfor(ll)
     de(l,2)=de(l,2)+qa(l,1)*(ra1*xafor(ll,1)+ra2*xafor(ll,2)+ra3*xafor(ll,3))
     de(l,3)=de(l,3)-qa(l,1)*(ra1*xafor(ll,1)+ra2*xafor(ll,2)+ra3*xafor(ll,3))
   end do
 end if
 end subroutine forcego

!====================================================================================
!=====COMPUTE CELL AREA OVER AEROFOILS
!====================================================================================
 subroutine wallArea
 implicit none
    integer :: bblock1,tblock1
    integer :: bblock2,tblock2
    real(nr) :: g11, g33, g13,coef
 
    select case(mbk)
    case(19)
    bblock1 = 6; tblock1 = 11
    bblock2 = 8; tblock2 = 13
    case(11)
    bblock1 = 4; tblock1 = 7
    bblock2 = 12; tblock2 = 13
    end select
 
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

!====================================================================================
!=====COMPUTE WALL NORMAL VECTOR OVER AEROFOILS
!====================================================================================
 subroutine walldir
 implicit none
    
    integer :: bblock1,tblock1
    integer :: bblock2,tblock2
    real(nr) :: coef,tmp
    real(nr), dimension(3) :: u,v,r
 
    select case(mbk)
    case(19)
    bblock1 = 6; tblock1 = 11
    bblock2 = 8; tblock2 = 13
    case(11)
    bblock1 = 4; tblock1 = 7
    bblock2 = 12; tblock2 = 13
    end select
    u=(/0,0,1/)
 
    if (wflag) then
    ! Find top or bottom
    if ((mb==bblock1).or.(mb==bblock2)) then
       coef=-one
    elseif ((mb==tblock1).or.(mb==tblock2)) then
       coef=one
    end if
    allocate(wnor(0:lcwall,3),wtan(0:lcwall,3))
    do ll = 0, lcwall; l=lwall(ll)
       tmp=coef/sqrt(etm(l,1)*etm(l,1)+etm(l,2)*etm(l,2)+etm(l,3)*etm(l,3))
       do m = 1, 3
       wnor(ll,m)=etm(l,m)*tmp
       end do
       v=wnor(ll,:)
       r=cross(u,v)
       wtan(ll,:)=r(:)*coef
    end do
    end if
 
 end subroutine walldir

!====================================================================================
!=====COMPUTE LIFT COEFFICIENT OVER AEROFOILS
!====================================================================================
 subroutine clpost(ele,nvar)
 
 use problemcase, only: span,delt1,delt2
 implicit none
    
    integer, intent(in) :: ele,nvar
    integer :: bblock,tblock,m,ll,dir
    real(nr) :: dynp,clp,clv,tcl
    logical :: flag
 
    clp=0;clv=0;flag=.false.;tcl=0;

 
    ! Define aerofoil blocks
    select case(mbk)
    case(19)
    select case(ele)
    case(1)
    bblock = 6; tblock = 11
    case(2)
    bblock = 8; tblock = 13
    end select
    case(11)
    select case(ele)
    case(1)
    bblock = 4; tblock = 7
    case(2)
    bblock = 12; tblock = 13
    end select
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
 
    do dir = 1, 2
    clp=0;clv=0;tcl=0;
       if (flag) then
          ! Compute Dynamic pressure
          dynp=two/(amachoo*amachoo*span)
          if (dir==1) then
             if (ispost) then
                call gettw(nvar)
             else
                call gettwrun
             end if
          end if
          if (wflag) then
             do ll = 0, lcwall; l=lwall(ll)
               clp=clp+(p(l)*wnor(ll,dir)*area(ll))
               if (nviscous==1) then
               clv=clv+tw(ll,dir)*area(ll)
               end if
             end do
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
       end if
    CALL MPI_ALLREDUCE(tcl,cl(ele,dir),1,MPI_REAL8,MPI_SUM,icom,ierr)
    end do
 
 end subroutine clpost

!====================================================================================
!=====COMPUTE WALL SHEAR STRESS AT RUNTIME
!====================================================================================
 subroutine gettwrun
 use subroutines3d, only: mpigo,deriv
 implicit none
 integer :: nn,ll


   if (nviscous==1) then
     de(:,1)=1/qa(:,1)
     de(:,2)=qa(:,2)*de(:,1)
     de(:,3)=qa(:,3)*de(:,1)
     de(:,4)=qa(:,4)*de(:,1)
     de(:,5)=gam*p(:)*de(:,1)
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

     rr(:,1)=de(:,5)
     m=5; call mpigo(ntdrv,nrone,n45no,m); call deriv(3,1); call deriv(2,1); call deriv(1,1)
     ss(:,1)=xim(:,1)*rr(:,1)+etm(:,1)*rr(:,2)+zem(:,1)*rr(:,3)
     ss(:,2)=xim(:,2)*rr(:,1)+etm(:,2)*rr(:,2)+zem(:,2)*rr(:,3)
     ss(:,3)=xim(:,3)*rr(:,1)+etm(:,3)*rr(:,2)+zem(:,3)*rr(:,3)

     fctr=2.0_nr/3
     rr(:,1)=-de(:,1)*yaco(:)
     rr(:,2)=gamm1prndtli*rr(:,1)
     de(:,5)=fctr*(txx(:)+tyy(:)+tzz(:))

     txx(:)=rr(:,1)*(2*txx(:)-de(:,5))
     tyy(:)=rr(:,1)*(2*tyy(:)-de(:,5))
     tzz(:)=rr(:,1)*(2*tzz(:)-de(:,5))
     txy(:)=rr(:,1)*(txy(:)+hzz(:))
     tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
     tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
 
      if (wflag) then
      if(.not.allocated(tw)) allocate(tw(0:lcwall,3))
      do ll = 0, lcwall; l=lwall(ll)
        tw(ll,1)=qa(l,1)*(txx(l)*wnor(ll,1)+txy(l)*wnor(ll,2)+tzx(l)*wnor(ll,3))!/reoo
        tw(ll,2)=qa(l,1)*(txy(l)*wnor(ll,1)+tyy(l)*wnor(ll,2)+tyz(l)*wnor(ll,3))!/reoo
        tw(ll,3)=qa(l,1)*(tzx(l)*wnor(ll,1)+tyz(l)*wnor(ll,2)+tzz(l)*wnor(ll,3))!/reoo
      end do
      end if
   end if
    
 end subroutine gettwrun

!====================================================================================
!=====COMPUTE WALL SHEAR STRESS FROM WRITTEN DATA
!====================================================================================
 subroutine gettw(nvar)
 use subroutines3d, only: mpigo,deriv
 implicit none
 integer, intent(in) :: nvar
 integer :: nn,ll

    if (wflag) then
    ! READ VARIABLES
       selectcase(output)
       case(0)
       nread=nrec+(totVar*nvar)
       do nn = 1, 5
        if (tecplot) then
        nread=nread+1; call tpostread(nread,lsta)
        else
        nread=nread+1; call postread(nread)
        end if
        qo(:,nn)=varr(:)
       end do
       case(1)
          call p3dread(gsflag=0,nout=nvar)
       end select
       p(:)=qo(:,5)
    if (nviscous==1) then
       de(:,1)=1/qo(:,1)
       de(:,2)=qo(:,2)
       de(:,3)=qo(:,3)
       de(:,4)=qo(:,4)
       de(:,5)=gam*p(:)*de(:,1)
       de(:,1)=srefp1dre*de(:,5)**1.5_nr/(de(:,5)+srefoo)
 
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
     rr(:,1)=-de(:,1)*yaco(:)
     de(:,5)=fctr*(txx(:)+tyy(:)+tzz(:))
 
     txx(:)=rr(:,1)*(2*txx(:)-de(:,5))
     tyy(:)=rr(:,1)*(2*tyy(:)-de(:,5))
     tzz(:)=rr(:,1)*(2*tzz(:)-de(:,5))
     txy(:)=rr(:,1)*(txy(:)+hzz(:))
     tyz(:)=rr(:,1)*(tyz(:)+hxx(:))
     tzx(:)=rr(:,1)*(tzx(:)+hyy(:))
 
       if(.not.allocated(tw)) allocate(tw(0:lcwall,3))
       do ll = 0, lcwall; l=lwall(ll)
         tw(ll,1)=qo(l,1)*(txx(l)*wnor(ll,1)+txy(l)*wnor(ll,2)+tzx(l)*wnor(ll,3))!/reoo
         tw(ll,2)=qo(l,1)*(txy(l)*wnor(ll,1)+tyy(l)*wnor(ll,2)+tyz(l)*wnor(ll,3))!/reoo
         tw(ll,3)=qo(l,1)*(tzx(l)*wnor(ll,1)+tyz(l)*wnor(ll,2)+tzz(l)*wnor(ll,3))!/reoo
       end do
    end if
    end if
    
 end subroutine gettw

!====================================================================================
! ====READ DATA FOR POST-PROCESSING FROM TECPLOT FILE
!====================================================================================
 subroutine tpostread(num,lsta)
 implicit none
 integer, intent (in) :: num,lsta
  lp=lpos(myid)+lsta
     lq=(num-1)*ltomb
     do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
        read(9,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi)
     end do; end do
 end subroutine tpostread

!====================================================================================
! ====READ DATA FOR POST-PROCESSING
!====================================================================================
 subroutine postread(num)
 implicit none
 integer, intent (in) :: num
 integer :: lp,lq,l,k,j
  lp=lpos(myid)
     lq=(num-1)*ltomb
     do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
        read(8,pos=nr*(lp+lq+lio(j,k))+1) varr(l:l+lxi)
     end do; end do
 end subroutine postread

!====================================================================================
!=====READ RESULTS PLOT3D
!====================================================================================
subroutine p3dread(gsflag,nout)

integer, intent(in) :: gsflag,nout
integer :: n,nfile
real(4) :: res
character(8) :: ctime

nfile=5

   selectcase(gsflag)
   case(1)
     open (unit=nfile, file='out/grid.xyz', access='stream',shared)
     lh=0
      read(nfile,pos=4*lh+1) mbk; mbk=mbk-1; lh=lh+1 ! Number of zones
      do mm = 0, mbk
         read(nfile,pos=4*lh+1) lximb(mm);lximb(mm)=lximb(mm)-1; lh=lh+1 ! IMax
         read(nfile,pos=4*lh+1) letmb(mm);letmb(mm)=letmb(mm)-1; lh=lh+1 ! JMax
         read(nfile,pos=4*lh+1) lzemb(mm);lzemb(mm)=lzemb(mm)-1; lh=lh+1 ! KMax
      end do
      lhmb(0)=lh
      do mm = 0, mbk-1
         lhmb(mm+1)=lhmb(mm)+3*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
      end do
      lp=lpos(myid)+lhmb(mb)
      ns=1; ne=3
      do n=ns,ne; lq=(n-1)*ltomb
      do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         read(nfile,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
      end do; end do
      ss(:,n)=varr(:)
      end do
      close(nfile)
     CALL MPI_BARRIER(icom,ierr)
       if (myid==0) then
       write(*,*) 'Grid read!'
       end if
   case(0)
      if (nout>ndata) then
         open (unit=nfile, file='out/solA.qa', access='stream',shared)
      else
         open (unit=nfile, file=ofiles(nout), access='stream',shared)
      end if
      lh=0
      read(nfile,pos=4*lh+1) mbk; mbk=mbk-1; lh=lh+1 ! Number of zones
      do mm = 0, mbk
         read(nfile,pos=4*lh+1) lximb(mm);lximb(mm)=lximb(mm)-1; lh=lh+1 ! IMax
         read(nfile,pos=4*lh+1) letmb(mm);letmb(mm)=letmb(mm)-1; lh=lh+1 ! JMax
         read(nfile,pos=4*lh+1) lzemb(mm);lzemb(mm)=lzemb(mm)-1; lh=lh+1 ! KMax
      end do
      lhmb(0)=lh
       do mm = 0, mbk-1
          lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
       end do
          lh=lhmb(mb)
          read(nfile,pos=4*lh+1) res; amachoo=res; lh=lh+1 ! Mach Number
          read(nfile,pos=4*lh+1) res; aoa=res; lh=lh+1  
          read(nfile,pos=4*lh+1) res; reoo=res; lh=lh+1 ! Reynolds Number
          read(nfile,pos=4*lh+1) res; timo=res; lh=lh+1 ! Time
       lp=lpos(myid)+lhmb(mb)+4
       ns=1; ne=5
       do n=ns,ne; lq=(n-ns)*ltomb
       do k=0,lze; do j=0,let; l=indx3(0,j,k,1)
         read(nfile,pos=4*(lp+lq+lio(j,k))+1) varr(l:l+lxi) ! 4-Bytes "Stream"
       end do; end do
       if (nout>ndata) then
       selectcase(n)
          case(1); qo(:,n)=varr(:)
          case(2,3,4); qo(:,n)=varr(:)
          case(5); qo(:,n)=varr(:)
       end select
       else
       selectcase(n)
          case(1); qo(:,n)=varr(:)
          ! If the data was written in conservatives variables change this
          case(2,3,4); qo(:,n)=varr(:)!/qo(:,1)
          case(5); qo(:,n)=varr(:)!&
          !gamm1*(varr(:)-half*qo(:,1)*(qo(:,2)*qo(:,2)+qo(:,3)*qo(:,3)+qo(:,4)*qo(:,4)))
       end select
       end if
       end do
       close(nfile)
       write(ctime,"(f8.4)") times(min(nout,ndata))
       if (myid==0) then
          if (nout>ndata) then
          write(*,"('Averaged Solution read! T= ',8a)") ctime 
          else
          write(*,"('Solution read! T= ',8a)") ctime 
          end if
       end if
   end select
         
end subroutine p3dread

!====================================================================================
! ====CROSS PRODUCT OF TWO VECTORS 
!====================================================================================
 function cross(u,v) result(r)
 real(nr), dimension(3), intent(in) :: u,v
 real(nr), dimension(3) :: r

    r(1)=(u(2)*v(3)-u(3)*v(2))
    r(2)=(u(3)*v(1)-u(1)*v(3))
    r(3)=(u(1)*v(2)-u(2)*v(1))
 end function cross

end module rpt
