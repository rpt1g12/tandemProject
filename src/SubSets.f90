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

     ! Allocate SubSet (SS) Flags
     allocate(ssFlag(tss)) 
     ! Allocate Communicators, ids and # processors
     allocate(sscom(tss),ssbcom(tss),ssid(tss),bssid(tss),ssnp(tss)) 
     ! Allocate SS lmx array and SS frequencies
     allocate(sslmx(tss),ssFreq(tss))
     ! Allocate SubArray Types and FileHandlers arrays
     allocate(ssq4arr(tss),ssq4fh(tss))
     ! Allocate SubArray Types flag array
     allocate(ssq4flag(tss)) 
     ! Allocate array of start and end indices
     allocate(lss0(tss),lssn(tss))


     do nss = 1, tss
           read(9,*) cinput,cinput,ssFreq(nss) 
        do nn = 1, 2+mb
           read(9,*) 
        end do
           read(9,*) cinput,ssRange(1,:),ssRange(2,:),ssRange(3,:)
        do nn = 0, mbk-mb-1
           read(9,*) 
        end do
     end do

     do nss = 1, tss
     if (tss.ge.1) then
        do i = 1, 3
           ssSize(i)=ssRange(i,2)-ssRange(i,1)+1
        end do

        ssFlag(nss)=.true.
        ssq4flag(nss)=.false.
        do m = 1, 3
        if (ssFlag(nss)) then
           if (((mpijke(m).ge.ssRange(m,1)).and.&
               (mpijks(m).lt.ssRange(m,1))).or.&
              ((mpijks(m).ge.ssRange(m,1)).and.&
               (mpijke(m).le.ssRange(m,2))).or.&
              ((mpijks(m).le.ssRange(m,2)).and.&
               (mpijke(m).gt.ssRange(m,2)))) then
               ssFlag(nss)=.true.
               color=1
           else
               ssFlag(nss)=.false.
               color=MPI_UNDEFINED
           end if
        end if
        end do

        ! rpt- Create SubSet communicator 
         CALL MPI_COMM_SPLIT(icom,color,myid,sscom(nss),ierr)   
        if (ssFlag(nss)) then
           if(.not.allocated(ssGSzs)) allocate(ssGSzs(0:mbk,3))
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
           sslmx(nss)=(ssLSize(1))*(ssLSize(2))*(ssLSize(3))-1
           ! rpt- Rank and Sizes for SubSet Communicator
           call MPI_COMM_RANK(sscom(nss),ssid(nss),ierr)
           call MPI_COMM_SIZE(sscom(nss),ssnp(nss),ierr)
           ! rpt- Create SubSet block communicator 
           CALL MPI_COMM_SPLIT(sscom(nss),mb,myid,ssbcom(nss),ierr)   
           call MPI_COMM_RANK(ssbcom(nss),bssid(nss),ierr)
           if ((bssid(nss)==0)) then
              if (ssid(nss)==0) then
                 do m = 1, mbk
                  CALL MPI_RECV(ssGSzs(m,:),3,MPI_INTEGER4,MPI_ANY_SOURCE,m,sscom(nss),ista,ierr)
                 end do
              else
                 CALL MPI_SEND(ssGSzs(mb,:),3,MPI_INTEGER4,0,mb,sscom(nss),ierr)
              end if
           end if
           CALL MPI_BCAST(ssGSzs,3*(mbk+1),MPI_INTEGER4,0,sscom(nss),ierr)
           do i = 1, tss
              write(cnum,"(i3,a)") i
              do ii = 1, 3
                 l=scan(cnum,' ')
                 if (l==0) exit
                 cnum(l:l)='0'
              end do
              chstr='out/ss'//cnum//'/'
              if (ssid(nss)==0) call system('mkdir -p '//chstr)
           end do
           if(.not.allocated(lss)) allocate(lss(0:sslmx(nss)))
           ll=0
           do kk = ssGStr(3)-mpijks(3),ssGEnd(3)-mpijks(3)
              do jj = ssGStr(2)-mpijks(2),ssGEnd(2)-mpijks(2)
                 do ii = ssGStr(1)-mpijks(1),ssGEnd(1)-mpijks(1)
                    l=indx3(ii,jj,kk,1)
                    lss(ll)=l
                    ll=ll+1;
                 end do
              end do
           end do
           lssn(nss)=ll-1;lss0(nss)=lssn(nss)-sslmx(nss)
           if (ssid(nss)==0) write(*,*) 'SubSets Ready to use!'
        end if
     else
       ssFlag(nss)=.false.
     end if !tss
     end do !nss=1,tss

  end subroutine ssSetUp
!====================================================================================
!=====  SUBSET Checking
!====================================================================================
subroutine ssCheck()
integer :: nss,n
   
   do nss = 1, tss
      if (ssFlag(nss)) then
         do n = 0, ssnp(nss)-1
         if (ssid(nss)==n) then
            write(*,"(a,i3,a,i3,a,i6,a,i6)") &
            'id:',myid,' ssid:',ssid(nss),'lss0:',lss0(nss),'lssn:',lssn(nss)
         end if
         end do
      end if
   end do
end subroutine ssCheck
!====================================================================================
!=====  SUBSET PLOT3D XYZ FILES WRITE
!====================================================================================
  subroutine wrP3dG_ss(mblkin,nss)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nss
     character(len=*),parameter :: fname='grid'
     character(len=:),allocatable :: lfname
     character(2) :: cout
     character(3) :: cnum
     character(10) :: cpath
     character(len=*),parameter :: cext='.xyz'
     integer :: n,l,i,lh,iolen,comid,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: fh,amode,garr
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf

     if (ssFlag(nss)) then
        ! rpt- Set default option to Multiblock
        if(present(mblkin)) then
           mblk=mblkin
        else
           mblk=1
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
           wrcom=sscom(nss)
           nbk=mbk
        case(0)
           write(cout,"(i2)") mb
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
           wrcom=ssbcom(nss)
           nbk=0
        case default
           if(ssid(nss)==0) write(*,*) 'Wrong multiblock option! Aborting...'
           CALL MPI_ABORT(icom,err,ierr)
        end select

        ! rpt- Ranks and Size in SubSet communicator
        call MPI_COMM_RANK(wrcom,comid,ierr)

        l=len(cpath)+len(fname)+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(cout)//cext
        if(comid==0) CALL MPI_FILE_DELETE(lfname,info,ierr)

        wrlen=3*(sslmx(nss)+1)
        if(.not.allocated(ssxyz4)) allocate(ssxyz4(0:sslmx(nss),3))
        do ll = 0, sslmx(nss); l=lss(ll)
           ssxyz4(ll,:)=ss(l,:)
        end do
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
        CALL MPI_FILE_WRITE_ALL(fh,ssxyz4,wrlen,MPI_REAL4,ista,ierr)
        CALL MPI_FILE_CLOSE(fh,ierr)
        CALL MPI_TYPE_FREE(garr,ierr)
        deallocate(ssxyz4)
        if (comid==0) then
           write(*,"('SSGrid ',i3,' written!')") nss
        end if
     end if !ssFlag
  end subroutine wrP3dG_ss
!====================================================================================
!=====  SUBSETS PLOT3D Q FILES WRITE
!====================================================================================
  subroutine wrP3dS_ss(mblkin,nss)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nss
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout,cnum
     character(8) :: ctime
     character(10) :: cpath
     character(len=*),parameter :: cext='.q'
     integer :: n,l,ll,ii,jj,kk,i,lh,iolen,comid,bcomid,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

     if (ssFlag(nss)) then
        ! rpt- Set default options
        if(present(mblkin)) then
           mblk=mblkin
        else
           mblk=1
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
           wrcom=sscom(nss)
           nbk=mbk
        case(0)
           write(cout,"(a,i2)") 'b',mb
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
           wrcom=ssbcom(nss)
           nbk=0
        case default
           if(ssid(nss)==0) write(*,*) 'Wrong multiblock option! Aborting...'
           CALL MPI_ABORT(icom,err,ierr)
        end select

        ! rpt- Ranks and Size in SubSet communicator
        call MPI_COMM_RANK(wrcom,comid,ierr)
        bcomid=bssid(nss)

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

        wrlen=5*(sslmx(nss)+1)
        if(.not.allocated(ssq4)) allocate(ssq4(0:sslmx(nss),5))

        do ll = 0, sslmx(nss); l=lss(ll)
        ssq4(ll,1)=qa(l,1)
        do i = 2, 4
           ssq4(ll,i)=((qa(l,i)/qa(l,1))+umf(i-1))
        end do
        ssq4(ll,5)=p(l)
        end do

        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.ssq4flag(nss)) then
           gsizes(:)=(/ssSize(:),5/)
           lsizes(:)=(/ssLSize(:),5/)
           starts(:)=(/ssStr(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                           MPI_ORDER_FORTRAN,MPI_REAL4,ssq4arr(nss),ierr) 
           CALL MPI_TYPE_COMMIT(ssq4arr(nss),ierr)
           ssq4flag(nss)=.true.
        end if
        
        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,ssq4fh(nss),ierr)

        lh=0
        if (comid==0) then
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            ibuf=ssGSzs(mm,1); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
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
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=aoa    ; offset=lh*iolen ! AoA
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=reoo   ; offset=lh*iolen ! Reynolds Number
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=timo   ; offset=lh*iolen ! Time
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
        end if
        disp=(lhmb(mb)+4)*iolen
        CALL MPI_FILE_SET_VIEW(ssq4fh(nss),disp,MPI_REAL4,ssq4arr(nss),'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL(ssq4fh(nss),ssq4,wrlen,MPI_REAL4,ista,ierr)
        CALL MPI_FILE_CLOSE(ssq4fh(nss),ierr)
        if (comid==0) then
           write(*,"('Subset Solution written! T= ',8a)") ctime 
        end if
        if (ndati.ge.ndata) then
           CALL MPI_TYPE_FREE(ssq4arr(nss),ierr)
           ssq4flag(nss)=.false.
        end if
     end if !ssFlag
  end subroutine wrP3dS_ss
!====================================================================================
!=====  SUBSETS PLOT3D Q FILES WRITE
!====================================================================================
  subroutine wrP3dP_ss(nout,mblkin,cname,nss)
     integer, intent(in),optional :: mblkin
     character(*), intent(in),optional :: cname
     integer, intent(in) :: nout,nss
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout,cnum,ncout
     character(8) :: ctime
     character(10) :: cpath
     character(len=*),parameter :: cext='.q'
     integer :: n,l,ll,ii,jj,kk,i,lh,iolen,comid,bcomid,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

     if (ssFlag(nss)) then
        ! rpt- Set default options
        if(present(mblkin)) then
           mblk=mblkin
        else
           mblk=1
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
           wrcom=sscom(nss)
           nbk=mbk
        case(0)
           write(cout,"(a,i2)") 'b',mb
           do i = 0, 1
              l=scan(cout,' ')
              if (l==0) exit
              cout(l:l)='0'
           end do
           wrcom=ssbcom(nss)
           nbk=0
        case default
           if(ssid(nss)==0) write(*,*) 'Wrong multiblock option! Aborting...'
           CALL MPI_ABORT(icom,err,ierr)
        end select

        if(present(cname)) then
           write(ncout,"(i3)") nout
           do i = 0, 2
              l=scan(ncout,' ')
              if (l==0) exit
              ncout(l:l)='0'
           end do
           ctime=trim(cname)//ncout
        else
           if (nout.le.ndata) then
              write(ctime,"(f8.4)") times(nout)
              do i = 0, 8
                 l=scan(ctime,' ')
                 if (l==0) exit
                 ctime(l:l)='0'
              end do
           else if (nout==ndata+1) then
              ctime='A'
           else if (nout==ndata+2) then
              ctime='RMS'
           end if
        end if

        ! rpt- Ranks and Size in SubSet communicator
        call MPI_COMM_RANK(wrcom,comid,ierr)
        bcomid=bssid(nss)

        l=len(cpath)+len(fname)+len(trim(adjustl(ctime)))+len(trim(cout))+len(cext)
        allocate(character(len=l) :: lfname)
        lfname=cpath//trim(fname)//trim(adjustl(ctime))//trim(cout)//cext
        if(comid==0) CALL MPI_FILE_DELETE(lfname,info,ierr)

        wrlen=5*(sslmx(nss)+1)
        if(.not.allocated(ssq4)) allocate(ssq4(0:sslmx(nss),5))

        do ll = 0, sslmx(nss); l=lss(ll)
           ssq4(ll,:)=qo(l,:)
        end do

        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.ssq4flag(nss)) then
           gsizes(:)=(/ssSize(:),5/)
           lsizes(:)=(/ssLSize(:),5/)
           starts(:)=(/ssStr(:),0/)
           CALL MPI_TYPE_CREATE_SUBARRAY(4,gsizes,lsizes,starts,&
                           MPI_ORDER_FORTRAN,MPI_REAL4,ssq4arr(nss),ierr) 
           CALL MPI_TYPE_COMMIT(ssq4arr(nss),ierr)
           ssq4flag(nss)=.true.
        end if
        
        CALL MPI_FILE_OPEN(wrcom,lfname ,amode ,info ,ssq4fh(nss),ierr)

        lh=0
        if (comid==0) then
         ibuf=nbk+1; offset=lh*iolen          ! Number of blocks
         CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         do l = 0, nbk
            mm=l+(1-mblk)*mb
            ibuf=ssGSzs(mm,1); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
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
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=aoa    ; offset=lh*iolen ! AoA
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=reoo   ; offset=lh*iolen ! Reynolds Number
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
            rbuf=timo   ; offset=lh*iolen ! Time
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,rbuf,1,MPI_REAL4,ista,ierr); lh=lh+1
        end if
        disp=(lhmb(mb)+4)*iolen
        CALL MPI_FILE_SET_VIEW(ssq4fh(nss),disp,MPI_REAL4,ssq4arr(nss),'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL(ssq4fh(nss),ssq4,wrlen,MPI_REAL4,ista,ierr)
        if (comid==0) then
           write(*,"('Subset Solution written! T= ',8a)") trim(ctime)//trim(cout) 
        end if
        CALL MPI_FILE_CLOSE(ssq4fh(nss),ierr)
     end if !ssFlag
  end subroutine wrP3dP_ss
end module subsets
