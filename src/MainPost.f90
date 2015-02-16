!*****
!***** 3D PARALLEL SOLVER
!*****

 program mainpost

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 use rpt
 use rptpost
 implicit none

!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

tecplot=.false.
call setup
 !===== INPUT RECORDING PARAMETERS
  open(9,file='data/post.dat',shared)
  read(9,*) cinput,ngridv
  read(9,*) cinput,ndata
  read(9,*) cinput,nrec
  read(9,*) cinput,nwrec
  close(9)
    totVar=5

  if (tecplot) then
  lsta=17+ngridv*45+(ndata+1)*25+24+nwrec+3+4*nwrec
  end if
 !===== INQUIRE RECORD LENGTH
  open(3,file=cpostdat,access='stream',shared)
  open(9,file=coutput,access='stream',shared)
  nread=0


!===== COMPUTE AVERAGE VALUES IF NOT AVAILABLE YET
! ltz=-1
! call average

!!===== SAVE AVERAGE VALUES INTO ARRAY
! nread=nwrec-totVar
! do nn=1,5
!    nread=nread+1
!    if (tecplot) then
!    call tpostread(nread,lsta)
!    else 
!    call postread(nread)
!    end if
!    qa(:,nn)=varr(:)
! end do

!===COMPUTE FORCE COEFFICIENT
!call clpost(1,2,1,.true.)

!==COMPUTE SPACE AVERAGE
!if (myid==11) then
!varr(:)=qa(:,5)
!call spcavg(2,3,0)
!end if

!==WRITE DATA SO IT CAN BE READ BY POST SUBROUTINE
if (tecplot) then
   call tdatawrite(lsta)
else
   call datawrite
end if
close(3)

!!===== WRITE TECPLOT FILE
 call post

!==COMPUTE WALL DISTANCES
call getwplus(1)
if (myid==11) then
    open(7,file='out/wplus.dat')
    write(7,"('x y z x+ y+ z+')") 
 do nn = 0, lcwall;l=lwall(nn)
 write(7,"(f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5,' ',f10.5)")&
       xyz(l,1),xyz(l,2),xyz(l,3),wplus(nn,1),wplus(nn,2),wplus(nn,3)
 end do
 close(7)
allocate(wvarr(0:lcwall))
 wvarr=wplus(:,2)
 call wavg(2,wall=.true.)
end if

!===== END OF JOB
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program mainpost
