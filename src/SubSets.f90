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
     read(9,*) cinput,tss ! Read total number of SubSets

     ! Allocate SubSets' ranges
     allocate(ssRange(3,2,tss))
     ! Allocate SubSets' sizes
     allocate(ssSize(3,tss),ssLSize(3,tss))
     allocate(ssGSzs(0:mbk,3,tss))
     ! Allocate SubSets' Starts and Ends
     allocate(ssGStr(3,tss),ssGEnd(3,tss))
     allocate(ssStr(3,tss),ssEnd(3,tss))
     ! Allocate Subset Flags
     allocate(ssFlag(tss),ssq4flag(tss),sswrsfg(tss))
     ssFlag(:)=.false.;ssq4flag(:)=.False.;sswrsfg(:)=.False.
     ! Allocate Comunicators array
     allocate(sscom(tss),ssbcom(tss))
     ! Allocate Ids
     allocate(ssid(tss),bssid(tss),ssnp(tss))
     ! Allocate maximum SubSets maximum sizes
     allocate(sslmx(tss))
     ! Allocate SubSets' frequency
     allocate(ssFreq(tss))
     ! Allocate output arrays
     allocate(nout_ss(tss),ndati_ss(tss),ssq4arr(tss),ssq4fh(tss))

     do n = 1, tss
           read(9,*) cinput,cinput,ssFreq(n)
        do nn = 1, 2+mb
           read(9,*) 
        end do
           read(9,*) cinput,ssRange(1,:,n),ssRange(2,:,n),ssRange(3,:,n)
        do nn = 0, mbk-mb-1
           read(9,*) 
        end do
     end do

     if (tss.ge.1) then
     do n = 1, tss
        do i = 1, 3
           ssSize(i,n)=ssRange(i,2,n)-ssRange(i,1,n)+1
        end do

        ssFlag(n)=.true.
        do m = 1, 3
        if (ssFlag(n)) then
           if (((mpijke(m).ge.ssRange(m,1,n)).and.&
               (mpijks(m).lt.ssRange(m,1,n))).or.&
              ((mpijks(m).ge.ssRange(m,1,n)).and.&
               (mpijke(m).le.ssRange(m,2,n))).or.&
              ((mpijks(m).le.ssRange(m,2,n)).and.&
               (mpijke(m).gt.ssRange(m,2,n)))) then
               ssFlag(n)=.true.
               color=1
           else
               ssFlag(n)=.false.
               color=MPI_UNDEFINED
           end if
        end if
        end do

        ! rpt- Create SubSet communicator 
        CALL MPI_COMM_SPLIT(icom,color,myid,sscom(n),ierr)   
        if (ssFlag(n)) then
           ! Global sizes
           ssGSzs(mb,:,n)=ssSize(:,n)
           ! Global starts
           ssGStr(:,n)=(/max(mpijks(1),ssRange(1,1,n)),&
                    max(mpijks(2),ssRange(2,1,n)),&
                    max(mpijks(3),ssRange(3,1,n))/)
           ! Global ends
           ssGEnd(:,n)=(/min(mpijke(1),ssRange(1,2,n)),&
                    min(mpijke(2),ssRange(2,2,n)),&
                    min(mpijke(3),ssRange(3,2,n))/)
           ssStr(:,n)=ssGStr(:,n)-ssRange(:,1,n)  ! Start Points
           ssEnd(:,n)=ssGEnd(:,n)-ssRange(:,1,n)  ! End points
           ssLSize(:,n)=ssGEnd(:,n)-ssGStr(:,n)+1 ! Local Sizes
           ! Total points
           sslmx(n)=(ssLSize(1,n))*(ssLSize(2,n))*(ssLSize(3,n))-1
           ! rpt- Rank and Sizes for SubSet Communicator
           call MPI_COMM_RANK(sscom(n),ssid(n),ierr) ! Id in SubSet
           call MPI_COMM_SIZE(sscom(n),ssnp(n),ierr) ! # pcs in SubSet
           ! rpt- Create SubSet block communicator 
           CALL MPI_COMM_SPLIT(sscom(n),mb,myid,ssbcom(n),ierr)   
           call MPI_COMM_RANK(ssbcom(n),bssid(n),ierr) ! Id in Block-SubSet
           if ((bssid(n)==0)) then
              if (ssid(n)==0) then
                 do m = 1, mbk
                  CALL MPI_RECV(ssGSzs(m,:,n),3,MPI_INTEGER4,MPI_ANY_SOURCE,m,sscom(n),ista,ierr)
                 end do
              else
                 CALL MPI_SEND(ssGSzs(mb,:,n),3,MPI_INTEGER4,0,mb,sscom(n),ierr)
              end if
           end if
           CALL MPI_BCAST(ssGSzs(0,1,n),3*(mbk+1),MPI_INTEGER4,0,sscom(n),ierr)
        end if
     end do
     else
       ssFlag=.false.
     end if !tss

     ! Allocate grid & index arrays
     allocate(lss(0:maxval(sslmx(:)),tss))
     if (tss.ge.1) then
     do n = 1, tss
        if (ssFlag(n)) then
           ! Create SubSet Folder
                write(cnum,"(i3,a)") n
           do ii = 1, 3
              l=scan(cnum,' ')
              if (l==0) exit
              cnum(l:l)='0'
           end do
                chstr='out/ss'//cnum//'/'
                if (ssid(n)==0) call system('mkdir -p '//chstr)

           ll=0
           do kk = ssGStr(3,n)-mpijks(3),ssGEnd(3,n)-mpijks(3)
              do jj = ssGStr(2,n)-mpijks(2),ssGEnd(2,n)-mpijks(2)
                 do ii = ssGStr(1,n)-mpijks(1),ssGEnd(1,n)-mpijks(1)
                    l=indx3(ii,jj,kk,1)
                    lss(ll,n)=l
                    ll=ll+1;
                 end do
              end do
           end do
           if (ssid(n)==0) write(*,"('SubSet ',i3,' ready to use!')") n
        end if
     end do
     end if !tss

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

     if (ssFlag(nss)) then
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
        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        if(allocated(ssxyz4)) then
           deallocate(ssxyz4)
           allocate(ssxyz4(sslmx(nss),3))
        end if

        do ll = 0, sslmx(nss); l=lss(ll,nss)
           ssxyz4(ll,:)=xyz4(l,:)
        end do

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        gsizes(:)=(/ssSize(:,nss),3/)
        lsizes(:)=(/ssLSize(:,nss),3/)
        starts(:)=(/ssStr(:,nss),0/)
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
            ibuf=ssGSzs(mm,1,nss); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2,nss); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3,nss); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+3*(ssGSzs(mm,1,nss))*(ssGSzs(mm,2,nss))*(ssGSzs(mm,3,nss))
        end do
        disp=lhmb(mb)*iolen
        CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,garr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL_BEGIN(fh,ssxyz4,wrlen,MPI_REAL4,ierr)
        CALL MPI_FILE_WRITE_ALL_END(fh,ssxyz4,ista,ierr)
        CALL MPI_FILE_CLOSE(fh,ierr)
        CALL MPI_TYPE_FREE(garr,ierr)
        if (comid==0) then
           write(*,"('SSGrid ',i3,' written!')") nss
        end if
     end if !ssFlag
  end subroutine wrP3dG_ss
!====================================================================================
!=====  SUBSETS PLOT3D Q FILES WRITE
!====================================================================================
  subroutine wrP3dS_ss(mblkin,nssin)
     integer, intent(in),optional :: mblkin,nssin
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

     ! rpt- Set default options
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

     if (ssFlag(nss)) then

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
        if(allocated(ssq4)) then
           deallocate(ssq4)
           allocate(ssq4(0:sslmx(nss),5))
        else
           allocate(ssq4(0:sslmx(nss),5))
        end if

        do ll = 0, sslmx(nss); l=lss(ll,nss)
        ssq4(ll,1)=qa(l,1)
        do i = 2, 4
           ssq4(ll,i)=((qa(l,i)/qa(l,1))+umf(i-1))
        end do
        ssq4(ll,5)=p(l)
        end do

        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.ssq4flag(nss)) then
           gsizes(:)=(/ssSize(:,nss),5/)
           lsizes(:)=(/ssLSize(:,nss),5/)
           starts(:)=(/ssStr(:,nss),0/)
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
            ibuf=ssGSzs(mm,1,nss); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2,nss); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3,nss); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*ssGSzs(mm,1,nss)*ssGSzs(mm,2,nss)*ssGSzs(mm,3,nss)
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
           write(*,"('Subset ',i3,' Solution written! T= ',8a)") nss,ctime 
        end if
        CALL MPI_FILE_CLOSE(ssq4fh(nss),ierr)
        if (ndati.ge.ndata) then
           CALL MPI_TYPE_FREE(ssq4arr(nss),ierr)
           ssq4flag(nss)=.false.
        end if
     end if !ssFlag
  end subroutine wrP3dS_ss
!====================================================================================
!=====  SUBSETS PLOT3D Q FILES WRITE
!====================================================================================
  subroutine wrP3dP_ss(nout,mblkin,cname,nssin)
     integer, intent(in),optional :: mblkin,nssin
     integer, intent(in) :: nout
     character(*), intent(in),optional :: cname
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout,cnum,ncout
     character(8) :: ctime
     character(10) :: cpath
     character(len=*),parameter :: cext='.q'
     integer :: n,l,ll,ii,jj,kk,i,lh,iolen,comid,bcomid,wrcom,nbk,err,mblk,nss
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode
     integer, dimension (4) :: gsizes,lsizes,starts
     integer(k4) :: ibuf
     real   (k4) :: rbuf

     ! rpt- Set default options
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

     if (ssFlag(nss)) then

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
        if(allocated(ssq4)) then
           deallocate(ssq4)
           allocate(ssq4(0:sslmx(nss),5))
        else
           allocate(ssq4(0:sslmx(nss),5))
        end if

        do ll = 0, sslmx(nss); l=lss(ll,nss)
           ssq4(ll,:)=qo(l,:)
        end do

        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

        CALL MPI_TYPE_EXTENT(MPI_INTEGER4,iolen,ierr)
        if (.not.ssq4flag(nss)) then
           gsizes(:)=(/ssSize(:,nss),5/)
           lsizes(:)=(/ssLSize(:,nss),5/)
           starts(:)=(/ssStr(:,nss),0/)
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
            ibuf=ssGSzs(mm,1,nss); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,2,nss); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(mm,3,nss); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*ssGSzs(mm,1,nss)*ssGSzs(mm,2,nss)*ssGSzs(mm,3,nss)
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
           write(*,"('Subset ',i3,' Solution written! T= ',8a)") nss,trim(ctime)//trim(cout) 
        end if
        CALL MPI_FILE_CLOSE(ssq4fh(nss),ierr)
     end if !ssFlag
  end subroutine wrP3dP_ss
end module subsets
