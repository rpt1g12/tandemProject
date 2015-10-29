!**********************
!***** RPT MODULE *****
!**********************

module  rpt

use mainvar3d
use subroutineso
use mpi

contains
!====================================================================================
!=====  PLOT3D F FILES WRITE
!====================================================================================
  subroutine wrP3dF(fname,nout,ndim,mblkin)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nout,ndim
     character(len=*), intent(in) :: fname
     character(len=:),allocatable :: lfname
     character(3) :: cout,cblk
     character(len=*),parameter :: cext='.f',cpath='out/'
     integer :: n,l,i,lh,iolen,foper,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: fh,amode,farr
     integer(k4) :: ibuf
     integer, dimension (4) :: gsizes,lsizes,starts

     ! rpt- Set default option to Multiblock
     if(present(mblkin)) then
        mblk=mblkin
     else
        mblk=1
     end if

     selectcase(mblk);
     case(1)
        cblk=''
        foper=0
        wrcom=icom
        nbk=mbk
     case(0)
        write(cblk,"(i2,a)") mb,'n'
        do i = 0, 1
           l=scan(cblk,' ')
           if (l==0) exit
           cblk(l:l)='0'
        end do
        foper=mo(mb)
        wrcom=bcom
        nbk=0
     case default
        if(myid==0) write(*,*) 'Wrong multiblock option! Aborting...'
        CALL MPI_ABORT(icom,err,ierr)
     end select
     write(cout,"(i3)") nout
     do i = 0, 2
        l=scan(cout,' ')
        if (l==0) exit
        cout(l:l)='0'
     end do
     l=len(cpath)+len(fname)+len(trim(cblk))+len(cout)+len(cext)
     allocate(character(len=l) :: lfname)
     lfname=cpath//trim(fname)//trim(cblk)//cout//cext
     if(myid==foper) CALL MPI_FILE_DELETE(lfname,info,ierr)

     wrlen=ndim*(lmx+1)
     amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

     CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
     gsizes(:)=(/mbijkl(:),ndim/)
     lsizes(:)=(/mpijkl(:),ndim/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL4,farr,ierr) 
     CALL MPI_TYPE_COMMIT(farr,ierr)
     
     if(myid==foper) CALL MPI_FILE_DELETE(lfname,info,ierr)
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
         ibuf=ndim; offset=lh*iolen        ! #Dimensions
         CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
      end do
     end if
     l=(1-mblk)*mb
     lhmb(l)=1+(nbk+1)*4
     do mm = 0, nbk-1
        lhmb(mm+1)=lhmb(mm)+ndim*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
     end do
     disp=lhmb(mb)*iolen
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,farr,'native',info,ierr)
     CALL MPI_FILE_WRITE_ALL_BEGIN(fh,fout,wrlen,MPI_REAL4,ierr)
     CALL MPI_FILE_WRITE_ALL_END(fh,fout,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(farr,ierr)
     if (myid==foper) then
        write(*,"(a,' funtion written!')") trim(fname)//cout
     end if
  end subroutine wrP3dF
    

!====================================================================================
!=====  PLOT3D Q FILES WRITE
!====================================================================================
  subroutine wrP3dS(mblkin)
     integer, intent(in),optional :: mblkin
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout
     character(8) :: ctime
     character(len=*),parameter :: cext='.q',cpath='out/'
     integer :: n,l,i,lh,iolen,foper,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

     if (wrsfg) then
        CALL MPI_FILE_WRITE_ALL_END(q4fh,q4,ista,ierr)
        CALL MPI_FILE_CLOSE(q4fh,ierr)
        wrsfg=.false.
     end if
     if(.not.wrsfg) then
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
           write(cout,"(a,i2)") 'b',mb
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
        write(ctime,"(f8.4)") timo
        do i = 0, 8
        l=scan(ctime,' ')
        if (l==0) exit
        ctime(l:l)='0'
        end do
        l=len(cpath)+len(fname)+len(trim(adjustl(ctime)))+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(adjustl(ctime))//trim(cout)//cext
        if(myid==foper) CALL MPI_FILE_DELETE(lfname,info,ierr)

        wrlen=5*(lmx+1)
        if(.not.allocated(q4)) allocate(q4(0:lmx,5))

           q4(:,1)=qa(:,1)
        do i = 2, 4
           q4(:,i)=((qa(:,i)/qa(:,1))+umf(i-1))
        end do
           q4(:,5)=p(:)

        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.q4flag) then
           gsizes(:)=(/mbijkl(:),5/)
           lsizes(:)=(/mpijkl(:),5/)
           starts(:)=(/mpijks(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                           MPI_ORDER_FORTRAN,MPI_REAL4,q4arr,ierr) 
           CALL MPI_TYPE_COMMIT(q4arr,ierr)
           q4flag=.true.
    end if

        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,q4fh,ierr)

       lh=0
        if (myid==foper) then
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            ibuf=lximb(mm)+1; offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=letmb(mm)+1; offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=lzemb(mm)+1; offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(q4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do
        if (myid==mo(mb)) then
           lh=lhmb(mb)
            rbuf=amachoo; offset=lh*iolen ! Mach Number
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=aoa    ; offset=lh*iolen ! AoA
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=reoo   ; offset=lh*iolen ! Reynolds Number
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=timo   ; offset=lh*iolen ! Time
            CALL MPI_FILE_WRITE_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
        end if
        disp=(lhmb(mb)+4)*iolen
        CALL MPI_FILE_SET_VIEW(q4fh,disp,MPI_REAL4,q4arr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL_BEGIN(q4fh,q4,wrlen,MPI_REAL4,ierr)
        wrsfg=.true.
        if (myid==foper) then
           write(*,"('Solution written! T= ',8a)") ctime 
        end if
        if (ndati.ge.ndata) then
           CALL MPI_FILE_WRITE_ALL_END(q4fh,q4,ista,ierr)
           CALL MPI_FILE_CLOSE(q4fh,ierr)
           CALL MPI_TYPE_FREE(q4arr,ierr)
           wrsfg=.false.
        end if
     end if !wrsfg
  end subroutine wrP3dS

!====================================================================================
!=====  PLOT3D XYZ FILES WRITE
!====================================================================================
  subroutine wrP3dG(mblkin)
     integer, intent(in),optional :: mblkin
     character(len=*),parameter :: fname='grid'
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

!===================================================================================
!=====  PLOT3D XYZ FILES READ
!===================================================================================
  subroutine rdP3dG(mblkin)
     integer, intent(in),optional :: mblkin
     character(len=*),parameter :: fname='grid'
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

     wrlen=3*(lmx+1)
     amode=MPI_MODE_RDONLY
     allocate(xyz4(0:lmx,3))

     CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
     gsizes(:)=(/mbijkl(:),3/)
     lsizes(:)=(/mpijkl(:),3/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL4,garr,ierr) 
     CALL MPI_TYPE_COMMIT(garr,ierr)
     
     if (fflag) then
        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,fh,ierr)

          lh=0
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            offset=lh*iolen ! IMax
            CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            lximb(mm)=ibuf-1;
            offset=lh*iolen ! JMax
            CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            letmb(mm)=ibuf-1;
            offset=lh*iolen ! KMax
            CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            lzemb(mm)=ibuf-1;
         end do
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+3*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
        end do
        disp=lhmb(mb)*iolen
        CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,garr,'native',info,ierr)
        CALL MPI_FILE_READ_ALL(fh,xyz4,wrlen,MPI_REAL4,ista,ierr)
        CALL MPI_FILE_CLOSE(fh,ierr)
        CALL MPI_TYPE_FREE(garr,ierr)
        do i = 1, 3
           ss(:,i)=xyz4(:,i)
        end do
        if (myid==foper) then
           write(*,"('Grid',i3,' Read!')") mb
        end if
     else
        ss(:,:)=0
     end if
  end subroutine rdP3dG

!===================================================================================
!=====  PLOT3D Q FILES READ
!===================================================================================
  subroutine rdP3dS(nout,mblkin)
     integer, intent(in) :: nout
     integer, intent(in),optional :: mblkin
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout
     character(8) :: ctime
     character(len=*),parameter :: cext='.q',cpath='out/'
     integer :: n,l,i,lh,iolen,foper,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

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
           write(cout,"(a,i2)") 'b',mb
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

        if (nout==ndata+1) then
           l=len('out/solA'//trim(cout)//'.qa')
           allocate(character(len=l) :: lfname)
           lfname='out/solA'//trim(cout)//'.qa'
        elseif (nout==ndata+2) then
           l=len('out/solRMS'//trim(cout)//'.qa')
           allocate(character(len=l) :: lfname)
           lfname='out/solRMS'//trim(cout)//'.qa'
        else 
           l=len(ofiles(nout))
           allocate(character(len=l) :: lfname)
           lfname=ofiles(nout)
        end if

        wrlen=5*(lmx+1)
        if(.not.allocated(q4)) allocate(q4(0:lmx,5))

        amode=MPI_MODE_RDONLY

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.q4flag) then
           gsizes(:)=(/mbijkl(:),5/)
           lsizes(:)=(/mpijkl(:),5/)
           starts(:)=(/mpijks(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                           MPI_ORDER_FORTRAN,MPI_REAL4,q4arr,ierr) 
           CALL MPI_TYPE_COMMIT(q4arr,ierr)
           q4flag=.true.
        end if

        if (fflag) then
           CALL MPI_FILE_OPEN(wrcom,trim(lfname) ,amode ,info ,q4fh,ierr)

           l=(1-mblk)*mb
           lhmb(l)=1+(nbk+1)*3
           do mm = 0, nbk-1
              lhmb(mm+1)=lhmb(mm)+4+5*(lximb(mm)+1)*(letmb(mm)+1)*(lzemb(mm)+1)
           end do
           lh=lhmb(mb)
           offset=lh*iolen ! Mach Number
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           amachoo=rbuf;
           offset=lh*iolen ! AoA
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           aoa=rbuf;
           offset=lh*iolen ! Reynolds Number
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           reoo=rbuf;
           offset=lh*iolen ! Time
           CALL MPI_FILE_READ_AT(q4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
           timo=rbuf;

           disp=(lhmb(mb)+4)*iolen
           CALL MPI_FILE_SET_VIEW(q4fh,disp,MPI_REAL4,q4arr,'native',info,ierr)
           CALL MPI_FILE_READ_ALL(q4fh,q4,wrlen,MPI_REAL4,ista,ierr)
           if (myid==foper) then
              if (nout==ndata+1) then
                 write(*,"('AVG Solution read!')") 
              else if(nout==ndata+2) then
                 write(*,"('RMS Solution read!')") 
              else
                 write(*,"('Solution read! T= ',f8.4)") times(nout)
              end if
           end if

           CALL MPI_FILE_CLOSE(q4fh,ierr)
           qo(:,:)=q4(:,:)
        else
           qo(:,:)=0
        end if
  end subroutine rdP3dS
!====================================================================================
!=====  WRITE RAW RESTART
!====================================================================================
  subroutine wrRsta()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp,offset
     integer :: amode,iolen
     integer, dimension (4) :: gsizes,lsizes,starts
     real(k8) :: rbuf
     integer(k4) :: ibuf

     if (wrrfg) then
        CALL MPI_FILE_WRITE_ALL_END(qfh,q8,ista,ierr)
        CALL MPI_FILE_CLOSE(qfh,ierr)
        wrrfg=.false.
     end if
     if (.not.wrrfg) then
        wrlen=5*(lmx+1)
        if(.not.allocated(q8)) allocate(q8(0:lmx,5))
        q8(:,:)=qa(:,:)
        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
        CALL MPI_TYPE_EXTENT(MPI_REAL8,iolen,ierr)

        if (.not.qflag) then
        gsizes(:)=(/mbijkl(:),5/)
        lsizes(:)=(/mpijkl(:),5/)
        starts(:)=(/mpijks(:),0/)
        CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                         MPI_ORDER_FORTRAN,MPI_REAL8,qarr,ierr) 
        CALL MPI_TYPE_COMMIT(qarr,ierr)
        qflag=.true.
       end if
        

        if(myid==mo(mb)) CALL MPI_FILE_DELETE(crestart,info,ierr)
        CALL MPI_FILE_OPEN(bcom,crestart,amode,info,qfh,ierr)
       lh=0
       if (myid==mo(mb)) then
            ibuf=n; offset=lh*iolen ! Iteration Number
            CALL MPI_FILE_WRITE_AT(qfh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ndt; offset=lh*iolen ! ?
            CALL MPI_FILE_WRITE_AT(qfh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            rbuf=dt; offset=lh*iolen ! Timestep
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
            rbuf=dts; offset=lh*iolen ! ?
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
            rbuf=dte; offset=lh*iolen ! ?
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
            rbuf=timo; offset=lh*iolen ! time
            CALL MPI_FILE_WRITE_AT(qfh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1
        else
            lh=lh+6
       end if
        disp=lh*iolen
        CALL MPI_FILE_SET_VIEW(qfh,disp,MPI_REAL8,qarr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL_BEGIN(qfh,q8,wrlen,MPI_REAL8,ierr)
        wrrfg=.true.
       if (myid==0) then
           write(*,"('Restart file written!')") 
       end if
        if (ndati==ndata) then
           CALL MPI_FILE_WRITE_ALL_END(qfh,q4,ista,ierr)
           CALL MPI_FILE_CLOSE(qfh,ierr)
           CALL MPI_TYPE_FREE(qarr,ierr)
           wrrfg=.false.
    end if
     end if !wrrfg
  end subroutine wrRsta

!====================================================================================
!=====  READ RAW RESTART
!====================================================================================
  subroutine rdRsta()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp,offset
     integer :: fh,amode,qarr,iolen
     integer, dimension (4) :: gsizes,lsizes,starts
     real(k8) :: rbuf
     integer(k4) :: ibuf

     if (myid==0) then
        write(*,"('Reading restart file..')") 
     end if
          
     wrlen=5*(lmx+1)
     amode=MPI_MODE_RDONLY
     CALL MPI_TYPE_EXTENT(MPI_REAL8,iolen,ierr)

     gsizes(:)=(/mbijkl(:),5/)
     lsizes(:)=(/mpijkl(:),5/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,qarr,ierr) 
     CALL MPI_TYPE_COMMIT(qarr,ierr)
     

     CALL MPI_FILE_OPEN(bcom,crestart,amode,info,fh,ierr)
     lh=0

         offset=lh*iolen ! Iteration Number
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1; n=ibuf
         offset=lh*iolen ! 10*(n/10)+1
         CALL MPI_FILE_READ_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1; ndt=ibuf
         offset=lh*iolen ! Timestep
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; dt=rbuf
         offset=lh*iolen ! ?
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; dts=rbuf
         offset=lh*iolen ! ?
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; dte=rbuf
         offset=lh*iolen ! time
         CALL MPI_FILE_READ_AT(fh,offset,rbuf,1,MPI_REAL8,ista,ierr); lh=lh+1; timo=rbuf

     disp=lh*iolen
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL8,qarr,'native',info,ierr)
     CALL MPI_FILE_READ_ALL(fh,qa,wrlen,MPI_REAL8,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(qarr,ierr)


  end subroutine rdRsta

!====================================================================================
!=====  READ RAW GRID
!====================================================================================
  subroutine rdGrid()
     integer(kind=MPI_OFFSET_KIND) :: wrlen,disp
     integer :: fh,amode,garr
     integer, dimension (4) :: gsizes,lsizes,starts

     wrlen=3*(lmx+1)
     amode=MPI_MODE_RDONLY

     gsizes(:)=(/mbijkl(:),3/)
     lsizes(:)=(/mpijkl(:),3/)
     starts(:)=(/mpijks(:),0/)
     CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL8,garr,ierr) 
     CALL MPI_TYPE_COMMIT(garr,ierr)
     
     disp=0

     CALL MPI_FILE_OPEN(bcom,cgrid,amode,info,fh,ierr)
     CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL8,garr,'native',info,ierr)
     CALL MPI_FILE_READ_ALL(fh,ss,wrlen,MPI_REAL8,ista,ierr)
     CALL MPI_FILE_CLOSE(fh,ierr)
     CALL MPI_TYPE_FREE(garr,ierr)

     if(myid==mo(mb)) CALL MPI_FILE_DELETE(cgrid,info,ierr)
  end subroutine rdGrid

!====================================================================================
! ====SET UP FORCING PARAMETERS
!====================================================================================
 subroutine forceup
 use problemcase, only: span
 
 ra0=-log(0.0001_k8)/(rfor**2); ra1=2*pi/span
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
   ra1=amfor*cos(ra0*48.76_k8*amachoo)/3
   ra2=amfor*cos(ra0*53.6_k8*amachoo)/3
   ra3=amfor*cos(ra0*53.6_k8*amachoo)/3
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
    integer(k4) :: bct,bcb,bcw,color
    real(k8) :: g11, g33, g13,coef
 
    ! Find parallel grid position
    ip=mod(myid-mo(mb),npc(mb,1))
    jp=mod((myid-mo(mb))/npc(mb,1),npc(mb,2))
    kp=mod((myid-mo(mb))/(npc(mb,1)*npc(mb,2)),npc(mb,3))
 
    ! Initialise values
    wflag=.false.
    g11=0;g33=g11;g13=g33
    color=MPI_UNDEFINED
 
    ! Find processors in contact with wall
    bct=nbc(1,2); bcb=nbc(0,2); bcw=20+5*nviscous
    if (bcb==bcw) then
    j=0; wflag=.true.
    elseif (bct==bcw) then
    j=ijk(1,2); wflag=.true.
    end if
 
    ll=-1; lcwall=(nbsize(2)-1)
    if(myid==0) color=1
    if (wflag) then
       color=1
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
 
    call MPI_COMM_SPLIT(icom,color,myid,wcom,ierr)

    if (color==1) then
       call MPI_COMM_SPLIT(wcom,mb,myid,bwcom,ierr)
    end if
 
 end subroutine wallArea

!====================================================================================
!=====COMPUTE WALL NORMAL VECTOR OVER AEROFOILS
!====================================================================================
 subroutine walldir
 implicit none
    
    integer(k4) :: bblock1,tblock1
    integer(k4) :: bblock2,tblock2
    integer(k4) :: bct,bcb,bcw
    real(k8) :: coef,tmp
    real(k8), dimension(3) :: u,r
 
    u=(/0.0_k8,0.0_k8,1.0_k8/)
 
    if (wflag) then
    ! Find top or bottom
       bct=nbc(1,2); bcb=nbc(0,2); bcw=20+5*nviscous
       if (bcb==bcw) then
       coef=one
       elseif (bct==bcw) then
       coef=-one
    end if
    allocate(wnor(0:lcwall,3),wtan(0:lcwall,3))
    do ll = 0, lcwall; l=lwall(ll)
       tmp=coef/sqrt(etm(l,1)*etm(l,1)+etm(l,2)*etm(l,2)+etm(l,3)*etm(l,3))
       do m = 1, 3
       wnor(ll,m)=etm(l,m)*tmp
       end do
       r=cross(u,wnor(ll,1:3))
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
    
    integer(k4), intent(in) :: ele,nvar
    integer(k4) :: bct,bcb,bcw,m,ll,dir
    real(k8) :: dynp,clp,clv,tcl
 
    clp=0;clv=0;;tcl=0;

    do dir = 1, 2
    clp=0;clv=0;tcl=0;
    if (ispost.and.(dir==1)) call gettw(nvar)
    if (wflag) then
      ! Compute Dynamic pressure
      dynp=two/(amachoo*amachoo*span)
      do ll = 0, lcwall; l=lwall(ll)
        clp=clp+(p(l)*wnor(ll,dir)*area(ll))
      end do
      if (nviscous==1) then
        if ((.not.ispost).and.(dir==1)) then
           if(.not.allocated(tw)) allocate(tw(0:lcwall,3))
           do ll = 0, lcwall; l=lwall(ll)
              tw(ll,1)=qa(l,1)*(txx(l)*wnor(ll,1)+txy(l)*wnor(ll,2)+tzx(l)*wnor(ll,3))
              tw(ll,2)=qa(l,1)*(txy(l)*wnor(ll,1)+tyy(l)*wnor(ll,2)+tyz(l)*wnor(ll,3))
           end do
        end if
        do ll = 0, lcwall; l=lwall(ll)
          clv=clv+tw(ll,dir)*area(ll)
        end do
      end if
      tcl=(clp+clv)*dynp
    end if
    if(wflag.or.(myid==0)) then
       if (ispost) then
          CALL MPI_REDUCE(tcl,cl(ele,dir),1,MPI_REAL8,MPI_SUM,0,bwcom,ierr)
       else
          CALL MPI_REDUCE(tcl,cl(ele,dir),1,MPI_REAL8,MPI_SUM,0,wcom,ierr)
       end if
    end if
    end do
 
 end subroutine clpost

!====================================================================================
!=====COMPUTE WALL SHEAR STRESS FROM WRITTEN DATA
!====================================================================================
 subroutine gettw(nvar)
 use subroutines3d, only: mpigo,deriv
 implicit none
 integer(k4), intent(in) :: nvar
 integer(k4) :: nn,ll

    if (wflag) then
    ! READ VARIABLES
       call rdP3dS(nvar,fmblk)
       !call p3dread(gsflag=0,nout=nvar)
       p(:)=qo(:,5)
    if (nviscous==1) then
       de(:,1)=1/qo(:,1)
       de(:,2)=qo(:,2)
       de(:,3)=qo(:,3)
       de(:,4)=qo(:,4)
       de(:,5)=gam*p(:)*de(:,1)
       de(:,1)=srefp1dre*de(:,5)**1.5_k8/(de(:,5)+srefoo)
 
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
 
     fctr=2.0_k8/3
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
         tw(ll,1)=qo(l,1)*(txx(l)*wnor(ll,1)+txy(l)*wnor(ll,2)+tzx(l)*wnor(ll,3))
         tw(ll,2)=qo(l,1)*(txy(l)*wnor(ll,1)+tyy(l)*wnor(ll,2)+tyz(l)*wnor(ll,3))
         tw(ll,3)=qo(l,1)*(tzx(l)*wnor(ll,1)+tyz(l)*wnor(ll,2)+tzz(l)*wnor(ll,3))
       end do
    end if
    end if
    
 end subroutine gettw

!====================================================================================
!=====READ RESULTS PLOT3D
!====================================================================================
subroutine p3dread(gsflag,nout)

integer(k4), intent(in) :: gsflag,nout
integer(k4) :: n,nfile
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
      if (nout==ndata+1) then
         open (unit=nfile, file='out/solA.qa', access='stream',shared)
      elseif (nout==ndata+2) then
         open (unit=nfile, file='out/solRMS.qa', access='stream',shared)
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
 real(k8), dimension(3), intent(in) :: u,v
 real(k8), dimension(3) :: r

    r(1)=(u(2)*v(3)-u(3)*v(2))
    r(2)=(u(3)*v(1)-u(1)*v(3))
    r(3)=(u(1)*v(2)-u(2)*v(1))
 end function cross

end module rpt
