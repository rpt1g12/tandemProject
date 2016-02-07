!**********************
!***** RPT MODULE *****
!**********************

module  subsets

use mainvar3d
use subroutineso
use mpi

contains
!====================================================================================
!=====  PLOT3D XYZ FILES WRITE
!====================================================================================
  subroutine wrP3dG_Sub(mblkin)
     integer, intent(in),optional :: mblkin
     character(len=*),parameter :: fname='SubGrid'
     character(len=:),allocatable :: lfname
     character(2) :: cout
     character(len=*),parameter :: cext='.xyz',cpath='out/'
     integer :: n,l,i,lh,iolen,foper,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: fh,amode,garr
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf

     ! rpt- Set default option to Multiblock
     if(present(mblkin)) then
        mblk=mblkin
     else
        mblk=1
     end if

     selectcase(mblk);
     case(1)
        cout=''
        foper=0
        wrcom=icom
        nbk=mbk
     case(0)
        write(cout,"(i2)") mb
        do i = 0, 1
           l=scan(cout,' ')
           if (l==0) exit
           cout(l:l)='0'
        end do
        foper=mo(mb)
        wrcom=bcom
        nbk=0
     case default
        if(myid==0) write(*,*) 'Wrong multiblock option! Aborting...'
        CALL MPI_ABORT(icom,err,ierr)
     end select
     l=len(cpath)+len(fname)+len(trim(cout))+len(cext)
     allocate(character(len=l) :: lfname)
     lfname=cpath//trim(fname)//trim(cout)//cext
     if(myid==foper) CALL MPI_FILE_DELETE(lfname,info,ierr)

     wrlen=3*(lmx+1)
     amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

     CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
     gsizes(:)=(/mbijkl(:),3/)
     lsizes(:)=(/mpijkl(:),3/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL4,garr,ierr) 
     CALL MPI_TYPE_COMMIT(garr,ierr)
     
     CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,fh,ierr)

     lh=0
     if (myid==foper) then
      ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
      CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
      do l = 0, nbk
         mm=l+(1-mblk)*mb
         ibuf=lximb(mm)+1; offset=lh*iolen ! IMax
         CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         ibuf=letmb(mm)+1; offset=lh*iolen ! JMax
         CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         ibuf=lzemb(mm)+1; offset=lh*iolen ! KMax
         CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
      end do
     end if
     l=(1-mblk)*mb
     lhmb(l)=1+(nbk+1)*3
     do mm = 0, nbk-1
        lhmb(mm+1)=lhmb(mm)+3*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
     end do
     disp=lhmb(mb)*iolen
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,garr,'native',info,ierr)
     CALL MPI_FILE_WRITE_ALL_BEGIN(fh,xyz4,wrlen,MPI_REAL4,ierr)
     CALL MPI_FILE_WRITE_ALL_END(fh,xyz4,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(garr,ierr)
     if (myid==foper) then
        write(*,"('Grid written!')")
     end if
  end subroutine wrP3dG

end module subsets
