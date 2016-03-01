!**************************
!***** SUBSETS MODULE *****
!**************************

module  subsets

use mainvar3d
use subroutineso
use mpi

contains
!====================================================================================
!=====  SUBSETS SETUP
!====================================================================================
  subroutine ssSetUp
     implicit none
     integer :: n,i,j,k,m,ll
     character(3) :: cnum
     character(10) :: chstr

     open(unit=9, file='inputs.dat')
     read(9,*) cinput,tss
     do n = 1, tss
           read(9,*) cinput,cinput,ssFreq 
        do nn = 1, 2+mb
           read(9,*) 
        end do
           read(9,*) cinput,ssRange(1,:),ssRange(2,:),ssRange(3,:)
        do nn = 0, mbk-mb-1
           read(9,*) 
        end do
     end do

     do i = 1, 3
        ssSize(i)=ssRange(i,2)-ssRange(i,1)+1
     end do

     ssFlag=.true.
     do m = 1, 2
     if (ssFlag) then
        if (((mpijke(m).ge.ssRange(m,1)).and.&
            (mpijks(m).lt.ssRange(m,1))).or.&
           ((mpijks(m).ge.ssRange(m,1)).and.&
            (mpijke(m).le.ssRange(m,2))).or.&
           ((mpijks(m).le.ssRange(m,2)).and.&
            (mpijke(m).gt.ssRange(m,2)))) then
            ssFlag=.true.
            color=1
        else
            ssFlag=.false.
            color=MPI_UNDEFINED
        end if
     end if
     end do

     ! rpt- Create SubSet communicator 
      CALL MPI_COMM_SPLIT(icom,color,myid,sscom,ierr)   
     if (ssFlag) then
        allocate(ssGSzs(0:mbk,3))
        ssGSzs(mb,:)=ssSize(:)
        ssGStr=(/max(mpijks(1),ssRange(1,1)),&
                   max(mpijks(2),ssRange(2,1)),&
                   max(mpijks(3),ssRange(3,1))/)
        ssGEnd=(/min(mpijke(1),ssRange(1,2)),&
                 min(mpijke(2),ssRange(2,2)),&
                 min(mpijke(3),ssRange(3,2))/)
        ssStr=ssGStr-ssRange(:,1)
        ssEnd=ssGEnd-ssRange(:,1)
        ssLSize=ssGEnd(:)-ssGStr(:)+1
        sslmx=(ssLSize(1))*(ssLSize(2))*(ssLSize(3))-1
        ! rpt- Rank and Sizes for SubSet Communicator
        call MPI_COMM_RANK(sscom,ssid,ierr)
        call MPI_COMM_SIZE(sscom,ssnp,ierr)
        ! rpt- Create SubSet block communicator 
        CALL MPI_COMM_SPLIT(sscom,mb,myid,ssbcom,ierr)   
        call MPI_COMM_RANK(ssbcom,bssid,ierr)
        if ((bssid==0)) then
           if (ssid==0) then
              do m = 1, mbk
               CALL MPI_RECV(ssGSzs(m,:),3,MPI_INTEGER4,MPI_ANY_SOURCE,m,sscom,ista,ierr)
              end do
           else
              CALL MPI_SEND(ssGSzs(mb,:),3,MPI_INTEGER4,0,mb,sscom,ierr)
           end if
        end if
        CALL MPI_BCAST(ssGSzs,3*(mbk+1),MPI_INTEGER4,0,sscom,ierr)
     end if

     do i = 1, tss
        write(cnum,"(i3,a)") i
        do ii = 1, 3
           l=scan(cnum,' ')
           if (l==0) exit
           cnum(l:l)='0'
        end do
        chstr='out/ss'//cnum//'/'
        if (ssid==0) call system('mkdir -p '//chstr)
     end do


  end subroutine ssSetUp
!====================================================================================
!=====  SUBSET PLOT3D XYZ FILES WRITE
!====================================================================================
  subroutine wrP3dG_ss(mblkin,nssin)
     integer, intent(in),optional :: mblkin
     integer, intent(in),optional :: nssin
     character(len=*),parameter :: fname='grid'
     character(len=:),allocatable :: lfname
     character(2) :: cout
     character(3) :: cnum
     character(10) :: cpath
     character(len=*),parameter :: cext='.xyz'
     integer :: n,l,i,lh,iolen,comid,wrcom,nbk,err,mblk,nss
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: fh,amode,garr
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf

     if (ssFlag) then
        ! rpt- Set default option to Multiblock
        if(present(mblkin)) then
           mblk=mblkin
        else
           mblk=1
        end if
        if(present(nssin)) then
           nss=nssin
        else
           nss=1
        end if
        write(cnum,"(i3,a)") nss
        do ii = 1, 3
           l=scan(cnum,' ')
           if (l==0) exit
           cnum(l:l)='0'
        end do
        cpath='out/ss'//cnum//'/'

        selectcase(mblk);
        case(1)
           cout=''
           wrcom=sscom
           nbk=mbk
        case(0)
           write(cout,"(i2)") mb
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
           wrcom=ssbcom
           nbk=0
        case default
           if(ssid==0) write(*,*) 'Wrong multiblock option! Aborting...'
           CALL MPI_ABORT(icom,err,ierr)
        end select

        ! rpt- Ranks and Size in SubSet communicator
        call MPI_COMM_RANK(wrcom,comid,ierr)

        l=len(cpath)+len(fname)+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(cout)//cext
        if(comid==0) CALL MPI_FILE_DELETE(lfname,info,ierr)

        wrlen=3*(sslmx+1)
        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        gsizes(:)=(/ssSize(:),3/)
        lsizes(:)=(/ssLSize(:),3/)
        starts(:)=(/ssStr(:),0/)
        CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,MPI_ORDER_FORTRAN,&
                                      MPI_REAL4,garr,ierr) 
        CALL MPI_TYPE_COMMIT(garr,ierr)
        
        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,fh,ierr)

        lh=0
        if (comid==0) then
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            ibuf=ssGSzs(mm,1); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+3*(ssGSzs(mm,1))*(ssGSzs(mm,2))*(ssGSzs(mm,3))
        end do
        disp=lhmb(mb)*iolen
        CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,garr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL_BEGIN(fh,ssxyz4,wrlen,MPI_REAL4,ierr)
        CALL MPI_FILE_WRITE_ALL_END(fh,ssxyz4,ista,ierr)
        CALL MPI_FILE_CLOSE(fh,ierr)
        CALL MPI_TYPE_FREE(garr,ierr)
        if (comid==0) then
           write(*,"('SSGrid written!')")
        end if
     end if !ssFlag
  end subroutine wrP3dG_ss
!====================================================================================
!=====  SUBSETS PLOT3D Q FILES WRITE
!====================================================================================
  subroutine wrP3dS_ss(mblkin,nssin)
     integer, intent(in),optional :: mblkin
     integer, intent(in),optional :: nssin
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout,cnum
     character(8) :: ctime
     character(10) :: cpath
     character(len=*),parameter :: cext='.q'
     integer :: n,l,ll,ii,jj,kk,i,lh,iolen,comid,bcomid,wrcom,nbk,err,mblk,nss
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

     if (ssFlag) then
     if (wrsfg) then
        CALL MPI_FILE_WRITE_ALL_END(ssq4fh,ssq4,ista,ierr)
        CALL MPI_FILE_CLOSE(ssq4fh,ierr)
        sswrsfg=.false.
     end if
     if(.not.sswrsfg) then
        ! rpt- Set default option to Multiblock
        if(present(mblkin)) then
           mblk=mblkin
        else
           mblk=1
        end if
        if(present(nssin)) then
           nss=nssin
        else
           nss=1
        end if
        write(cnum,"(i3,a)") nss
        do ii = 1, 3
           l=scan(cnum,' ')
           if (l==0) exit
           cnum(l:l)='0'
        end do
        cpath='out/ss'//cnum//'/'

        selectcase(mblk);
        case(1)
           cout=''
           wrcom=sscom
           nbk=mbk
        case(0)
           write(cout,"(a,i2)") 'b',mb
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
           wrcom=ssbcom
           nbk=0
        case default
           if(ssid==0) write(*,*) 'Wrong multiblock option! Aborting...'
           CALL MPI_ABORT(icom,err,ierr)
        end select

        ! rpt- Ranks and Size in SubSet communicator
        call MPI_COMM_RANK(wrcom,comid,ierr)
        bcomid=bssid

        write(ctime,"(f8.4)") timo
        do i = 0, 8
        l=scan(ctime,' ')
        if (l==0) exit
        ctime(l:l)='0'
        end do
        l=len(cpath)+len(fname)+len(trim(adjustl(ctime)))+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(adjustl(ctime))//trim(cout)//cext
        if(comid==0) CALL MPI_FILE_DELETE(lfname,info,ierr)

        wrlen=5*(sslmx+1)
        if(.not.allocated(ssq4)) allocate(ssq4(0:sslmx,5))

        do ll = 0, sslmx; l=lss(ll)
        ssq4(ll,1)=qa(l,1)
        do i = 2, 4
           ssq4(ll,i)=((qa(l,i)/qa(l,1))+umf(i-1))
        end do
        ssq4(ll,5)=p(l)
        end do

        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.ssq4flag) then
           gsizes(:)=(/ssSize(:),5/)
           lsizes(:)=(/ssLSize(:),5/)
           starts(:)=(/ssStr(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                           MPI_ORDER_FORTRAN,MPI_REAL4,ssq4arr,ierr) 
           CALL MPI_TYPE_COMMIT(ssq4arr,ierr)
           ssq4flag=.true.
        end if
        
        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,ssq4fh,ierr)

        lh=0
        if (comid==0) then
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_WRITE_AT(ssq4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            ibuf=ssGSzs(mm,1); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*ssGSzs(mm,1)*ssGSzs(mm,2)*ssGSzs(mm,3)
        end do
        if (bcomid==0) then
           lh=lhmb(mb)
            rbuf=amachoo; offset=lh*iolen ! Mach Number
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=aoa    ; offset=lh*iolen ! AoA
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=reoo   ; offset=lh*iolen ! Reynolds Number
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=timo   ; offset=lh*iolen ! Time
            CALL MPI_FILE_WRITE_AT(ssq4fh,offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
        end if
        disp=(lhmb(mb)+4)*iolen
        CALL MPI_FILE_SET_VIEW(ssq4fh,disp,MPI_REAL4,ssq4arr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL_BEGIN(ssq4fh,ssq4,wrlen,MPI_REAL4,ierr)
        sswrsfg=.true.
        if (comid==0) then
           write(*,"('Subset Solution written! T= ',8a)") ctime 
        end if
        if (ndati.ge.ndata) then
           CALL MPI_FILE_WRITE_ALL_END(ssq4fh,ssq4,ista,ierr)
           CALL MPI_FILE_CLOSE(ssq4fh,ierr)
           CALL MPI_TYPE_FREE(ssq4arr,ierr)
           sswrsfg=.false.
        end if
     end if !sswrsfg
     end if !ssFlag
  end subroutine wrP3dS_ss
end module subsets
