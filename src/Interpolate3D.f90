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
 real(nr) :: tol


!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL
    wts=MPI_WTIME()
    !iflag=.true.;igflag=.true.

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))

    ll=max(npro,12); allocate(ista(MPI_STATUS_SIZE,ll),ireq(ll))

    inquire(iolength=ll) real(1.0,kind=ieee32); nrecs=ll
    inquire(iolength=ll) real(1.0,kind=ieee64); nrecd=ll


   call interSetUp
   call setup(1)
lxii=lxi;leti=let;lzei=lze
write(*,*) myid,lxi,let,lze
lxiio=lxio;letio=leto;lzeio=lzeo
allocate(qb(0:lmx,5))
call prepareArrays
   allocate(xyz2(0:lmx,3),ixis(0:lmx,3))
   if (igflag) then
       call getGrid
        do i = 1, 3
           xyz2(:,i)=ss(:,i)
        end do
       !call wrIGrid
       call writeGrid
       call wrP3dG
   else
      !call rdIGrid
      call readGrid
!        do i = 1, 3
!           ss(:,i)=xyz2(:,i)
!        end do
!      call wrP3dG
   end if
   call rdRsta
   ndati=ndata;call wrP3dS

call deallocateArrays

if (iflag) then
   call setup(0)
   write(*,*) myid,lxi,let,lze
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
   call rdIRsta
   if (myid==0) then
   write(*,"('Last time step was at:',f7.4)") timo
   end if

   ! Get boundary values
   do n = 1, 2
   ra0=minval(xyz(:,n),1)
   ra1=maxval(xyz(:,n),1)
   CALL MPI_ALLREDUCE(ra0,bounds(0,n),1,MPI_REAL8,MPI_MIN,icom,ierr)
   CALL MPI_ALLREDUCE(ra1,bounds(1,n),1,MPI_REAL8,MPI_MAX,icom,ierr)
   bounds(0,n)=bounds(0,n)-2*sml
   bounds(1,n)=bounds(1,n)+2*sml
   end do
   bounds(0,3)=-half*span;bounds(1,3)=half*span
   ! Set outside values
   outside(1)=rhooo
   outside(2:4)=0.0_nr
   outside(5)=poo*hamm1/rhooo
   !do n = 1, 5
   !   call getDeri(n)
   !   call interpolate(n,tol)
   !end do
   call spanCopy
   
   
   call deallocateArrays
   call setup(1)
   call prepareArrays
       allocate(fout(0:lmx,5))
       fout(:,:)=qb(:,:)
       call wrP3dF('inter',0,5)
       deallocate(fout)
   call wrIRsta
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
         write(*,'("Interpolation time was ",f6.2," hours")') wtime/(3600_nr*npro)
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
