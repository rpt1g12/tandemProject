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


!===== PREPARATION FOR PARALLEL COMPUTING

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)

    mpro=npro-1; icom=MPI_COMM_WORLD; info=MPI_INFO_NULL

    allocate(lxim(0:mpro),letm(0:mpro),lzem(0:mpro),lpos(0:mpro),vmpi(0:mpro))
    allocate(ista(MPI_STATUS_SIZE,12))

call setup(250,250,250,240,80,15)
lxii=lxi;leti=let;lzei=lze
allocate(qb(0:lmx,5))
call prepareArrays
call getGrid
     allocate(xyz2(0:lmx,3),ixis(0:lmx,3))
     do i = 1, 3
        xyz2(:,i)=ss(:,i)
     end do

call deallocateArrays

call setup(200,200,200,200,40,15)
call prepareArrays
call getGrid
     allocate(xyz(0:lmx,3))
     if (myid==0) then
      write(*,*) lxi,let,lze,lmx  
     end if
     do i = 1, 3
        xyz(:,i)=ss(:,i)
     end do
call getMetrics
call readRestart
do n = 1, 5
call getDeri(n)
do i = 0, lxii
   do j = 0, leti
      do k = 0, lzei;l2=indx4(i,j,k,1)
         if (n==1) then
         xs(:)=(/xyz2(l2,1),xyz2(l2,2),xyz2(l2,3)/)
         start(:)=(/nint(real(lxi*i/lxii,nr)),nint(real(let*j/leti,nr)),nint(real(lze*k/lzei,nr))/)
         start(:)=max(start(:),(/0.0_nr,0.0_nr,0.0_nr/));
         start(:)=min(start(:),(/real(lxi,nr),real(let,nr),real(lze,nr)/))
         xin(:)=(/start(1),start(2),start(3)/)
         hxi(:)=(/0,0,0/)
         limit=sqrt(real(lxi**2+let**2+lze**2))/&
               sqrt(real(lxii**2+leti**2+lzei**2))
         err = 2.0_nr;err1=1
         do while(err1>sml)
         xin=xin+hxi
         xin(:)=max(xin(:),(/0.0_nr,0.0_nr,0.0_nr/));
         xin(:)=min(xin(:),(/real(lxi,nr),real(let,nr),real(lze,nr)/))
         thisxi=xin(1);thiset=xin(2);thisze=xin(3)
         varr=xyz(:,1);xn=trilinr(thisxi,thiset,thisze)
         varr=xyz(:,2);yn=trilinr(thisxi,thiset,thisze)
         varr=xyz(:,3);zn=trilinr(thisxi,thiset,thisze)
         hxn=xs(:)-(/xn,yn,zn/)
         varr=xxi;xxin=trilinr(thisxi,thiset,thisze)
         varr=xet;xetn=trilinr(thisxi,thiset,thisze)
         varr=xze;xzen=trilinr(thisxi,thiset,thisze)
         varr=yxi;yxin=trilinr(thisxi,thiset,thisze)
         varr=yet;yetn=trilinr(thisxi,thiset,thisze)
         varr=yze;yzen=trilinr(thisxi,thiset,thisze)
         varr=zxi;zxin=trilinr(thisxi,thiset,thisze)
         varr=zet;zetn=trilinr(thisxi,thiset,thisze)
         varr=zze;zzen=trilinr(thisxi,thiset,thisze)
         jaco(1,:)=(/xxin,yxin,zxin/)
         jaco(2,:)=(/xetn,yetn,zetn/)
         jaco(3,:)=(/xzen,yzen,zzen/)
         hxi(1)=sum(jaco(1,:)*hxn(:))
         hxi(2)=sum(jaco(2,:)*hxn(:))
         hxi(3)=sum(jaco(3,:)*hxn(:))
         tol=(xxin*yetn*zzen+xetn*yzen*zxin+xzen*yxin*zzen)-&
             (zxin*yetn*xzen+zetn*yzen*xxin+zzen*yxin*xetn)
         tol=tol**(1.0_nr/3.0_nr)
         err1=sqrt(hxn(1)**2+hxn(2)**2+hxn(3)**2)
         err=err1/tol
         end do
         ixis(l2,:)=(/thisxi,thiset,thisze/)
         end if
         qb(l2,n)=htrilinr(ixis(l2,1),ixis(l2,2),ixis(l2,3))
      end do
 if (myid==0) then
 write(*,*) i,j
 end if
   end do
end do
end do
CALL MPI_BARRIER(icom,ierr)

call deallocateArrays
call setup(250,250,250,240,80,15)
call writeRestart




!===== END OF JOB
 if(myid==0) then
    write(*,*) "Finished."
 end if
 call MPI_FINALIZE(ierr)

end program maininter
