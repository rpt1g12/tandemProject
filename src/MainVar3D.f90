!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

 module mainvar3d

 implicit none

!===== CONSTANT PARAMETERS
!
 integer,parameter :: icray=1

 integer,parameter :: int32=selected_int_kind(9),int64=selected_int_kind(18)
 integer,parameter :: ieee32=selected_real_kind(6,37),ieee64=selected_real_kind(15,307)
 integer,parameter :: ni=int32,nr=ieee64
 integer,parameter :: k4=int32,k8=ieee64

 integer(kind=ni),parameter :: ntdrv=0,ntflt=1,nrall=0,nrone=1,n45no=0,n45go=1
 integer(kind=ni),parameter :: lmd=11,lmf=8,lmp=max(lmd,lmf),mbci=3
 integer(kind=ni),parameter :: liofs=16,liofl=24

 character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

 real(kind=nr),parameter :: quarter=0.25_nr
 real(kind=nr),parameter :: half=0.5_nr,zero=0,one=1,two=2,three=3,four=4,five=5
 real(kind=nr),parameter :: one_three=one/three,two_three=two/three
 real(kind=nr),parameter :: pi=acos(-one),halfpi=pi/two,twopi=two*pi,sqrt2=sqrt(two),sqrt2i=one/sqrt2
 real(kind=nr),parameter :: sml=1.0e-6_nr,free=1.0e+6_nr
 real(kind=nr),parameter :: gam=1.4_nr,gamm1=gam-one,ham=one/gam,hamm1=one/gamm1
 real(kind=nr),parameter :: prndtl=0.71_nr,gamm1prndtli=one/(gamm1*prndtl)
 ! rpt- LES Constants
 real(k8),parameter :: tprndtl=0.99_k8,tgamm1prndtli=1/(gamm1*tprndtl)

 real(kind=nr),parameter :: alpha=0.5862704032801503_nr
 real(kind=nr),parameter :: beta=0.09549533555017055_nr
 real(kind=nr),parameter :: aa=0.6431406736919156_nr
 real(kind=nr),parameter :: ab=0.2586011023495066_nr
 real(kind=nr),parameter :: ac=0.007140953479797375_nr

 real(kind=nr),parameter :: beta20=0.03250008295108466_nr
 real(kind=nr),parameter :: alpha21=0.3998040493524358_nr
 real(kind=nr),parameter :: alpha23=0.7719261277615860_nr
 real(kind=nr),parameter :: beta24=0.1626635931256900_nr
 real(kind=nr),parameter :: a20=-0.1219006056449124_nr
 real(kind=nr),parameter :: a21=-0.6301651351188667_nr
 real(kind=nr),parameter :: a23=0.6521195063966084_nr
 real(kind=nr),parameter :: a24=0.3938843551210350_nr
 real(kind=nr),parameter :: a25=0.01904944407973912_nr
 real(kind=nr),parameter :: a26=-0.001027260523947668_nr

 real(kind=nr),parameter :: alpha10=0.08360703307833438_nr
 real(kind=nr),parameter :: alpha12=2.058102869495757_nr
 real(kind=nr),parameter :: beta13=0.9704052014790193_nr
 real(kind=nr),parameter :: a10=-0.3177447290722621_nr
 real(kind=nr),parameter :: a12=-0.02807631929593225_nr
 real(kind=nr),parameter :: a13=1.593461635747659_nr
 real(kind=nr),parameter :: a14=0.2533027046976367_nr
 real(kind=nr),parameter :: a15=-0.03619652460174756_nr
 real(kind=nr),parameter :: a16=0.004080281419108407_nr

 real(kind=nr),parameter :: alpha01=5.912678614078549_nr
 real(kind=nr),parameter :: beta02=3.775623951744012_nr
 real(kind=nr),parameter :: a01=-3.456878182643609_nr
 real(kind=nr),parameter :: a02=5.839043358834730_nr
 real(kind=nr),parameter :: a03=1.015886726041007_nr
 real(kind=nr),parameter :: a04=-0.2246526470654333_nr
 real(kind=nr),parameter :: a05=0.08564940889936562_nr
 real(kind=nr),parameter :: a06=-0.01836710059356763_nr

!===== ALLOCATABLE MAIN ARRAYS

 integer(kind=ni),dimension(:,:),allocatable :: npc,lio
 integer(kind=ni),dimension(:),allocatable :: li,lcsz,lctz,mxc
 integer(kind=ni),dimension(:),allocatable :: lxim,letm,lzem,lpos
 integer(kind=ni),dimension(:),allocatable :: lximb,letmb,lzemb,mo
 integer(kind=int64),dimension(:),allocatable :: lhmb

 real(kind=nr),dimension(:,:),allocatable :: qo,qa,qb,de
 real(kind=nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

 real(kind=nr),dimension(:,:),allocatable :: xim,etm,zem
 real(kind=nr),dimension(:),allocatable :: p,yaco
 real(kind=nr),dimension(:),allocatable :: sbcc

 real(kind=nr),dimension(:,:),allocatable :: rr,ss

 real(kind=nr),dimension(:,:),allocatable :: xu,yu
 real(kind=nr),dimension(:,:),allocatable :: xl,yl
 real(kind=nr),dimension(:),allocatable :: sa,sb

 real(kind=nr),dimension(:,:),allocatable :: ran,sit,ait
 real(kind=nr),dimension(:),allocatable :: xit,yit,zit
 real(kind=nr),dimension(:),allocatable :: asz,bsz,atz
 real(kind=nr),dimension(:),allocatable :: times,vmpi,rpex

 real(kind=nr),dimension(:,:,:),pointer :: drva,drvb,send,recv,cm

 real(kind=nr),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
 real(kind=nr),dimension(:,:,:),allocatable,target :: drvb1,drvb2,drvb3
 real(kind=nr),dimension(:,:,:),allocatable,target :: send1,send2,send3
 real(kind=nr),dimension(:,:,:),allocatable,target :: recv1,recv2,recv3
 real(kind=nr),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3

 real(kind=ieee32),dimension(:,:),allocatable :: varm
 real(kind=ieee32),dimension(:),allocatable :: varr,vart,vara,varb,vmean,varmin,varmax

 character(13),dimension(:),allocatable :: ctecplt,cthead
 character(4),dimension(:),allocatable :: cfilet
 character(7),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

 integer(kind=ni),dimension(0:1,0:1,3) :: ndf
 integer(kind=ni),dimension(3,3) :: ijk
 integer(kind=ni),dimension(0:1,3) :: nbc,ncd,nsz
 integer(kind=ni),dimension(0:4) :: no
 integer(kind=ni),dimension(-2:2) :: ilag
 integer(kind=ni),dimension(3) :: nbsize,nbcs,nbce,ncds,ncde,ms,me
 integer(kind=ni) :: lxio,leto,lzeo,lxi,let,lze,lmx,lim,nrecs,nrecd
 integer(kind=ni) :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,kp,l,lh,ll,lp,lq,ltomb
 integer(kind=ni) :: m,ma,mb,mm,mp,mq,mbk,mps,mpe,n,ndt,nn,nk,ns,ne,np,nt,nz,ndati,nsigi,nout,nfile
 integer(kind=ni) :: nts,nscrn,nsgnl,ndata,nviscous,nkrk,nsmf,nfskp,nrestart,nextrabc,nextgcic
 integer(kind=int64) :: nlmx,llmb,llmo,lis,lie,ljs,lje

 real(kind=nr),dimension(0:lmp,0:1,0:1) :: pbci,pbco
 real(kind=nr),dimension(-2:2,0:2,0:1) :: albed,albef
 real(kind=nr),dimension(mbci,mbci) :: cbca,cbcs
 real(kind=nr),dimension(5,5) :: xt
 real(kind=nr),dimension(0:5,0:2) :: abc
 real(kind=nr),dimension(0:1,0:1) :: pbcot
 real(kind=nr),dimension(0:lmp) :: sap
 real(kind=nr),dimension(mbci) :: rbci,sbci
 real(kind=nr),dimension(5) :: cha,dha
 real(kind=nr),dimension(-2:2) :: alag,blag,tlag
 real(kind=nr),dimension(3) :: ve,dm,rv,uoo,umf,dudtmf
 real(kind=nr),dimension(0:2) :: fam,fbm,fcm
 real(kind=nr) :: alphf,betf,fa,fb,fc
 real(kind=nr) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
 real(kind=nr) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,fltk,fltkbc,dto
 real(kind=nr) :: rhooo,poo,aoo,amachoo,srefoo,srefp1dre
 real(kind=nr) :: dt,dts,dte,dtk,dtko,dtsum,timo,tsam,wts,wte,wtime
 real(kind=nr) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

 character(1),dimension(0:4) :: cno
 character(3) :: cnzone,cndata
 character(5) :: cnnode
 character(7) :: czone
 character(17) :: coutput
 character(16) :: cinput
 character(16) :: cgrid
 character(18) :: cdata,cturb
 character(19) :: crestart
 character(256) :: cstring

!===== INTEGER VARIABLES FOR MPI COMMANDS

 integer(kind=ni),dimension(:,:),allocatable :: ista
 integer(kind=ni),dimension(:),allocatable :: ireq
 integer(kind=ni) :: ir,mpro,npro,myid,itag,info,icom,ierr

!===== INTEGER VARIABLES FOR RECORDING BY RPT

 integer(k4) :: ngrec,nwrec,nread,totVar

!===== VARIABLES FOR FORCING BY RPT

 real(k8)  ::  xfor,yfor,rfor,amfor,tsfor,tefor
 integer(k4)   ::  lfor
 integer(k4) ,allocatable,dimension(:)  ::  lcfor
 real(k8),allocatable,dimension(:,:)  ::  xafor,yafor,bfor

!===== MPI-IO VARIABLES BY RPT
 integer(k4), dimension (3) :: mpc,mbijkl,mpijkl,mpijks,mpijke
 integer(k4), dimension (:), allocatable :: ibegin,jbegin,kbegin
 integer :: bcom
 integer :: q4arr,q4fh,qarr,qfh
 logical :: qflag=.false.,gflag=.false.,q4flag=.false.,faflag=.false.
 logical :: wrsfg=.false.,wrrfg=.false.
 real(k8), dimension(:,:), allocatable :: q8
 real(k4), dimension(:,:), allocatable :: fout,xyz4,q4
!===== VARIABLES FOR WALL OPERATIONS BY RPT
 integer(k4) :: lcwall
 integer(k4), dimension(:), allocatable ::lwall
 real(k8), dimension(:), allocatable ::area
 real(k8), dimension(:,:), allocatable ::wnor,wtan,tw,pna
 logical :: wflag
 real(k8), dimension(:,:), allocatable,target :: xyz
 integer(k4) :: wcom,bwcom

!===== SUBSETS VARIABLES
 integer (k4), dimension (:,:,:), allocatable :: ssRange
 integer (k4), dimension (:,:), allocatable :: ssSize,ssLSize,ssGStr,ssGEnd,ssStr,ssEnd
 integer (k4), dimension (:,:,:), allocatable :: ssGSzs
 integer (k4), dimension (:), allocatable :: lss
 integer (k4), dimension (:), allocatable :: lss0,lssn
 integer (k4), dimension (:), allocatable :: sslmx,ssFreq
 integer (k4), dimension (:), allocatable :: sscom,ssbcom,ssid,bssid,ssnp,ssmbk,ssmb
 integer (k4) :: nss,tss
 integer (k4) :: color
 integer (k4), dimension (:), allocatable :: nout_ss,ndati_ss
 integer, dimension(:), allocatable :: ssq4arr,ssq4fh
 logical, dimension(:), allocatable :: ssFlag
 logical, dimension(:), allocatable :: ssq4flag
 integer(k4), dimension(:,:), allocatable :: ssblks
 real(k4), dimension(:), allocatable :: ssxyz4,ssq4

!===== POST-PROCESSING VARIABLES BY RPT

 integer(k4) :: lsta
 logical :: tecplot,ispost
 real(k8), dimension(:,:), allocatable :: wplus
 real(k8), dimension(:), allocatable :: wvarr
 character(:),dimension(:),allocatable :: ofiles 
 real(k8), dimension(2,2) :: cl
 real(k8), dimension(2,2,2) :: clh
 real(k8), dimension(2) :: clrng

 !=== ipost.dat input flags
 integer :: fparallel,fmblk
 integer :: favg
 integer :: fcoef
 integer :: floc
 integer :: fwss,fcf,fwplus
 integer :: fcurl
 integer :: fstrip
 integer :: fprobcirc
 integer :: fijkmax

 real(k8),dimension(:),allocatable :: delt
 real(k8), dimension (:,:), allocatable :: svarr
! real(k8), dimension (:,:,:,:), allocatable :: qxyz
 logical :: fflag

 ! Cylinder probes variables
 real(k8), dimension(:,:), allocatable :: xyprob
 integer, dimension(:,:), allocatable :: nlprob
 real(k8), dimension(:,:,:), allocatable :: nklprob
 integer, dimension(:), allocatable :: mprob,probcom,lprob
 real(k8), dimension(:), allocatable :: aprob,vprob
 integer :: lcprob,nprob
 real(k8) :: orprob
 real(k8), dimension(2) :: sprob,eprob
 logical, dimension(:), allocatable :: probflag

 logical :: intgflag
 integer :: intgcom
 integer(k4) :: lcintg
 real(k8),dimension(:),allocatable :: vintg,aintg
 integer(k4),dimension(:),allocatable :: lintg
 real :: rdis,xpos,ypos
 integer :: fintg,atk

 ! Maximum probe variables
 real(k4), dimension(:,:,:,:), allocatable :: mval
 integer(k4) , dimension(:,:,:), allocatable :: maxpos
 real(k4) , dimension(:,:,:), allocatable :: maxxyz
 real(k4) :: rout
 integer(k4), dimension(:,:), allocatable :: nose
 integer(k4), dimension(:), allocatable :: xarr,zarr
 integer(k4), dimension(2) :: msizes
 integer(k4) :: mblock,mxst,mxend,mxsz,mzst,mzend,mzsz

 ! Dummy variables
 integer(k4) :: idum
 real(k4) :: rdum
 real(k8) :: ddum
 character(128) :: cdum
 character(1) :: onechr

!===== ADITIONAL INPUT VARIABLES BY RPT
 integer(k4)  :: nto,iwrec
 integer(k4)  :: forcing,LES
 real(k8) :: tgustd,tguste
 real(k8) :: tps,tp,tpe,aoa0,aoa1,aoa,raoa
 real(k8) :: smago1,smago2
 integer(k4)  :: output,ogrid,osol,oblock

 integer(k4) :: bkx,bky,bkz
 integer(k4), dimension(:),allocatable:: lxibk,letbk,lzebk


!=====

 end module mainvar3d

!*****
