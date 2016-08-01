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
     integer :: n,i,j,k,m,ll,idum,gsize,qsize,l,ii,jj,kk
     character(3) :: cnum
     character(10) :: chstr

     open(unit=9, file='inputs.dat')
     read(9,*) cinput,tss

     ! Allocate SubSet (SS) Flags
     allocate(ssFlag(tss)) 
     ! Allocate SS ranges
     allocate(ssRange(3,2,tss))
     ! Allocate SS Sizes, Starts, and Ends
     allocate(ssGSzs(3,0:mbk,tss))
     allocate(ssSize(3,tss),ssLSize(3,tss),ssGStr(3,tss),ssGEnd(3,tss),ssStr(3,tss))
     ! Allocate Communicators, ids and # processors
     allocate(sscom(tss),ssbcom(tss),ssmcom(tss),ssid(tss),bssid(tss),ssnp(tss),ssmb(tss)) 
     ! Allocate SS lmx array and SS frequencies
     allocate(sslmx(tss),ssFreq(tss))
     sslmx(:)=0;ssFreq(:)=0
     ! Allocate SubArray Types and FileHandlers arrays
     allocate(ssq4arr(tss),ssq4fh(tss))
     ! Allocate SubArray Types flag array
     allocate(ssq4flag(tss)); ssq4flag(:)=.False.
     ! Allocate array of start and end indices
     allocate(lss0(tss),lssn(tss))


     do nss = 1, tss
           read(9,*) cinput,cinput,ssFreq(nss) 
        do nn = 1, 2+mb
           read(9,*) 
        end do
           read(9,*) cinput,ssRange(1,:,nss),ssRange(2,:,nss),ssRange(3,:,nss)
        do nn = 0, mbk-mb-1
           read(9,*) 
        end do
     end do

     do nss = 1, tss
     if (tss.ge.1) then
        idum=1
        do i = 1, 3
           ssSize(i,nss)=ssRange(i,2,nss)-ssRange(i,1,nss)+1
           if (ssSize(i,nss)<0) ssSize(i,nss)=0
           idum=idum*ssSize(i,nss)
        end do

        if (idum==0) then
           ssFlag(nss)=.false.
           color=MPI_UNDEFINED
        else
           ssFlag(nss)=.true.
        end if

        do m = 1, 3
        if (ssFlag(nss)) then
           if (((mpijke(m).ge.ssRange(m,1,nss)).and.&
               (mpijks(m).lt.ssRange(m,1,nss))).or.&
              ((mpijks(m).ge.ssRange(m,1,nss)).and.&
               (mpijke(m).le.ssRange(m,2,nss))).or.&
              ((mpijks(m).le.ssRange(m,2,nss)).and.&
               (mpijke(m).gt.ssRange(m,2,nss)))) then
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
           ssGStr(:,nss)=(/max(mpijks(1),ssRange(1,1,nss)),&
                      max(mpijks(2),ssRange(2,1,nss)),&
                      max(mpijks(3),ssRange(3,1,nss))/)
           ssGEnd(:,nss)=(/min(mpijke(1),ssRange(1,2,nss)),&
                    min(mpijke(2),ssRange(2,2,nss)),&
                    min(mpijke(3),ssRange(3,2,nss))/)
           ssStr(:,nss)=ssGStr(:,nss)-ssRange(:,1,nss)
           ssEnd(:,nss)=ssGEnd(:,nss)-ssRange(:,1,nss)
           ssLSize(:,nss)=ssGEnd(:,nss)-ssGStr(:,nss)+1
           sslmx(nss)=(ssLSize(1,nss))*(ssLSize(2,nss))*(ssLSize(3,nss))-1
           ! rpt- Rank and Sizes for SubSet Communicator
           call MPI_COMM_RANK(sscom(nss),ssid(nss),ierr)
           call MPI_COMM_SIZE(sscom(nss),ssnp(nss),ierr)
           ! rpt- Create SubSet block communicator 
           CALL MPI_COMM_SPLIT(sscom(nss),mb,myid,ssbcom(nss),ierr)   
           call MPI_COMM_RANK(ssbcom(nss),bssid(nss),ierr)
           if (bssid(nss)==0) then
              color=1
           else 
              color=MPI_UNDEFINED
           end if
           ! rpt- Create SubSet masters communicator 
           CALL MPI_COMM_SPLIT(sscom(nss),color,myid,ssmcom(nss),ierr)   
           ! rpt- Number of blocks involved in SubSet
           call MPI_COMM_SIZE(ssmcom(nss),idum,ierr)
           ssmb(nss)=idum-1
           if ((bssid(nss)==0)) then
              if (ssid(nss)==0) then
                 ssGSzs(:,mb,nss)=ssSize(:,nss)
                 do m = 1, ssmb(nss)
                    CALL MPI_RECV(ssGSzs(1,m,nss),3,MPI_INTEGER4,MPI_ANY_SOURCE,&
                                  m,sscom(nss),ista,ierr)
                    do i = 1, 3
                       if (ssGSzs(i,m,nss)<0) ssGSzs(i,m,nss)=0
                    end do
                 end do
              else
                 CALL MPI_SEND(ssSize(1,nss),3,MPI_INTEGER4,0,mb,sscom(nss),ierr)
              end if
              end if
           CALL MPI_BCAST(ssGSzs(:,:,nss),3*(ssmb(nss)+1),MPI_INTEGER4,0,sscom(nss),ierr)
           end if
     else
       ssFlag(nss)=.false.
     end if !tss
     end do !nss=1,tss


     idum=sum(sslmx(:)+1)-1
     gsize=idum*3-1
     qsize=idum*5-1
     allocate(lss(0:idum),ssxyz4(0:gsize),ssq4(0:qsize))
     ll=0
     do nss = 1, tss
        if (ssFlag(nss)) then
           write(cnum,"(i3,a)") nss
              do ii = 1, 3
              i=scan(cnum,' ')
              if (i==0) exit
              cnum(i:i)='0'
              end do
              chstr='out/ss'//cnum//'/'
           if (ssid(nss)==0) call system('mkdir -p '//chstr)

           lss0(nss)=ll
           do kk = ssGStr(3,nss)-mpijks(3),ssGEnd(3,nss)-mpijks(3)
              do jj = ssGStr(2,nss)-mpijks(2),ssGEnd(2,nss)-mpijks(2)
                 do ii = ssGStr(1,nss)-mpijks(1),ssGEnd(1,nss)-mpijks(1)
                    l=indx3(ii,jj,kk,1)
                    lss(ll)=l
                    ll=ll+1;
                 end do
              end do
           end do
           lssn(nss)=ll-1
           if (ssid(nss)==0) write(*,"(a,i3,a)") 'SubSet ',nss,' Ready to use!'
        end if
     end do

  end subroutine ssSetUp
!====================================================================================
!=====  SUBSET Checking
!====================================================================================
subroutine ssCheck()
integer :: nss,n,m
   
   do nss = 1, tss
      if (ssFlag(nss)) then
         if (bssid(nss)==0) then
            write(*,*) &
            'ss:',nss,' id:',myid,' bk:',mb
         end if
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
     integer :: n,l,llss,i,ii,lh,iolen,comid,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: fh,amode,garr,idum
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
           nbk=ssmb(nss)
           if (nbk==0) mblk=0
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
        idum=0
        do i = 1, nss-1
           idum=idum+3*(sslmx(i)+1)
        end do

        do i = 1, 3
           do ll = 0, sslmx(nss); l=lss(ll+lss0(nss))
              ii=ll+(i-1)*(sslmx(nss)+1)+idum
              ssxyz4(ii)=ss(l,i)
           end do
        end do
        amode=IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)

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
            ibuf=ssGSzs(1,mm,nss); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(2,mm,nss); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(3,mm,nss); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(fh,offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+3*(ssGSzs(1,mm,nss))*(ssGSzs(2,mm,nss))*(ssGSzs(3,mm,nss))
        end do
        disp=lhmb(mb)*iolen
        CALL MPI_FILE_SET_VIEW(fh,disp,MPI_REAL4,garr,'native',info,ierr)
        CALL MPI_FILE_WRITE_ALL(fh,ssxyz4(idum),wrlen,MPI_REAL4,ista,ierr)
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
  subroutine wrP3dS_ss(mblkin,nss)
     integer, intent(in),optional :: mblkin
     integer, intent(in) :: nss
     character(len=*),parameter :: fname='solT'
     character(len=:),allocatable :: lfname
     character(3) :: cout,cnum
     character(8) :: ctime
     character(10) :: cpath
     character(len=*),parameter :: cext='.q'
     integer :: n,l,ll,llss,ii,jj,kk,i,lh,iolen,comid,bcomid,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode,idum
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
           nbk=ssmb(nss)
           if (nbk==0) mblk=0
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

        idum=0
        do i = 1, nss-1
           idum=idum+5*(sslmx(i)+1)
        end do

        do ll = 0, sslmx(nss); l=lss(ll+lss0(nss))
           ii=ll+idum
           ssq4(ii)=qa(l,1)
        end do
        do i = 2, 4
           do ll = 0, sslmx(nss); l=lss(ll+lss0(nss))
              ii=ll+(i-1)*(sslmx(nss)+1)+idum
              ssq4(ii)=((qa(l,i)/qa(l,1))+umf(i-1))
        end do
        end do
        do ll = 0, sslmx(nss); l=lss(ll+lss0(nss))
           ii=ll+4*(sslmx(nss)+1)+idum
           ssq4(ii)=p(l)
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
            ibuf=ssGSzs(1,mm,nss); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(2,mm,nss); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(3,mm,nss); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*ssGSzs(1,mm,nss)*ssGSzs(2,mm,nss)*ssGSzs(3,mm,nss)
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
        CALL MPI_FILE_WRITE_ALL(ssq4fh(nss),ssq4(idum),wrlen,MPI_REAL4,ista,ierr)
        CALL MPI_FILE_CLOSE(ssq4fh(nss),ierr)
        if (comid==0) then
           write(*,"('Subset Solution',i3,' written! T= ',8a)") nss,ctime 
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
     integer :: n,l,ll,llss,ii,jj,kk,i,lh,iolen,comid,bcomid,wrcom,nbk,err,mblk
     integer(kind=MPI_OFFSET_KIND) :: wrlen,offset,disp
     integer :: amode,idum
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
           nbk=ssmb(nss)
           if (nbk==0) mblk=0
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

        idum=0
        do i = 1, nss-1
           idum=idum+5*(sslmx(i)+1)
        end do

        do i = 1, 5
           do ll = 0, sslmx(nss); l=lss(ll+lss0(nss))
              ii=ll+(i-1)*(sslmx(nss)+1)+idum
              ssq4(ii)=qo(l,i)
           end do
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
            ibuf=ssGSzs(1,mm,nss); offset=lh*iolen ! IMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(2,mm,nss); offset=lh*iolen ! JMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
            ibuf=ssGSzs(3,mm,nss); offset=lh*iolen ! KMax
            CALL MPI_FILE_WRITE_AT(ssq4fh(nss),offset,ibuf,1,MPI_INTEGER4,ista,ierr); lh=lh+1
         end do
        end if
        l=(1-mblk)*mb
        lhmb(l)=1+(nbk+1)*3
        do mm = 0, nbk-1
           lhmb(mm+1)=lhmb(mm)+4+5*ssGSzs(1,mm,nss)*ssGSzs(2,mm,nss)*ssGSzs(3,mm,nss)
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
        CALL MPI_FILE_WRITE_ALL(ssq4fh(nss),ssq4(idum),wrlen,MPI_REAL4,ista,ierr)
        if (comid==0) then
           write(*,"('Subset Solution',i3,' written! T= ',8a)") nss,trim(ctime)//trim(cout) 
        end if
        CALL MPI_FILE_CLOSE(ssq4fh(nss),ierr)
     end if !ssFlag
  end subroutine wrP3dP_ss
end module subsets
