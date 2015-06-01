!*****
!***** 3D INTERPOLATION SCHEME
!*****

 program maininter

 use mpi
 use subroutineso
 use subroutines3d
 use problemcase
 use rpt
 use rptinter
 implicit none
 logical :: flag
 real(nr) :: tol


!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL
    wts=MPI_WTIME()
    flag=.true.

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

   call setup(210,315,210,210,90,100)
lxii=lxi;leti=let;lzei=lze
lxiio=lxio;letio=leto;lzeio=lzeo
allocate(qb(0:lmx,5))
call prepareArrays
     allocate(xyz2(0:lmx,3),ixis(0:lmx,3))
if (.true.) then
call getGrid
     do i = 1, 3
        xyz2(:,i)=ss(:,i)
     end do
call writeGrid
else
call readGrid
end if

call deallocateArrays

if (flag) then
   call setup(160,320,160,160,80,100)
   color=mb;ncom=mb
   call prepareArrays
   call getGrid
   call MPI_COMM_SPLIT(icom,color,myid,ncom,ierr)
        allocate(xyz(0:ltomb-1,3))
        do i = 1, 3
           varr=ss(:,i);call joinBlock
           xyz(:,i)=lvarr
        end do
   call getMetrics
   ra0 = minval(yaco,1)
   CALL MPI_ALLREDUCE(ra0,tol,1,MPI_REAL8,MPI_MIN,ncom,ierr)
   tol=abs(1.0_nr/tol**(1.0_nr/3.0_nr))
   write(*,*) tol,myid
   call readRestart
   
   do n = 1, 5
      call getDeri(n)
      call interpolate(n,tol)
   end do
   
   
   call deallocateArrays
   call setup(210,315,210,210,90,100)
   call writeRestart
    if(myid==mo(mb)) then
       write(*,*) "Finished",mb
    end if

   wte=MPI_WTIME(); res=wte-wts
   call MPI_ALLREDUCE(res,wtime,1,MPI_REAL8,MPI_SUM,icom,ierr)
   if(myid==0) then
      if (wtime/(3600_nr*npro)<1_nr) then
         if (wtime/(60_nr*npro)<1_nr) then
            write(*,'("Interpolation time was ",f6.2," seconds")') wtime/(1_nr*npro)
         else
            write(*,'("Interpolation time was ",f6.2," minutes")') wtime/(60_nr*npro)
         end if
      else
         write(*,'("Simulation time was ",f6.2," hours")') wtime/(3600_nr*npro)
      end if
   end if
end if


!===== END OF JOB
 close(8)
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program maininter
