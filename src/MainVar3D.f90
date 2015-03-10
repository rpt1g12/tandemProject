!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

 module mainvar3d

 implicit none

!===== CONSTANT PARAMETERS

 integer,parameter :: nr=kind(0.0d0),ntdrv=0,ntflt=1,nrall=0,nrone=1,n45no=0,n45go=1
 integer,parameter :: lmd=11,lmf=8,lmp=max(lmd,lmf),mfbi=3,mbci=3
 integer,parameter :: liofs=16,liofl=24

 character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

 real(nr),parameter :: pi=3.141592653589793_nr
 real(nr),parameter :: zero=0,one=1,half=0.5_nr,sqrt2=sqrt(2.0_nr),sqrt2i=1/sqrt2,sml=1.0e-6_nr
 real(nr),parameter :: two=2,quarter=0.25_nr
 real(nr),parameter :: gam=1.4_nr,gamm1=gam-1,ham=1/gam,hamm1=1/gamm1
 real(nr),parameter :: prndtl=0.71_nr,gamm1prndtli=1/(gamm1*prndtl)

 real(nr),parameter :: alpha=0.5862704032801503_nr
 real(nr),parameter :: beta=0.09549533555017055_nr
 real(nr),parameter :: aa=0.6431406736919156_nr
 real(nr),parameter :: ab=0.2586011023495066_nr
 real(nr),parameter :: ac=0.007140953479797375_nr

 real(nr),parameter :: beta20=0.03250008295108466_nr
 real(nr),parameter :: alpha21=0.3998040493524358_nr
 real(nr),parameter :: alpha23=0.7719261277615860_nr
 real(nr),parameter :: beta24=0.1626635931256900_nr
 real(nr),parameter :: a20=-0.1219006056449124_nr
 real(nr),parameter :: a21=-0.6301651351188667_nr
 real(nr),parameter :: a23=0.6521195063966084_nr
 real(nr),parameter :: a24=0.3938843551210350_nr
 real(nr),parameter :: a25=0.01904944407973912_nr
 real(nr),parameter :: a26=-0.001027260523947668_nr

 real(nr),parameter :: alpha10=0.08360703307833438_nr
 real(nr),parameter :: alpha12=2.058102869495757_nr
 real(nr),parameter :: beta13=0.9704052014790193_nr
 real(nr),parameter :: a10=-0.3177447290722621_nr
 real(nr),parameter :: a12=-0.02807631929593225_nr
 real(nr),parameter :: a13=1.593461635747659_nr
 real(nr),parameter :: a14=0.2533027046976367_nr
 real(nr),parameter :: a15=-0.03619652460174756_nr
 real(nr),parameter :: a16=0.004080281419108407_nr

 real(nr),parameter :: alpha01=5.912678614078549_nr
 real(nr),parameter :: beta02=3.775623951744012_nr
 real(nr),parameter :: a01=-3.456878182643609_nr
 real(nr),parameter :: a02=5.839043358834730_nr
 real(nr),parameter :: a03=1.015886726041007_nr
 real(nr),parameter :: a04=-0.2246526470654333_nr
 real(nr),parameter :: a05=0.08564940889936562_nr
 real(nr),parameter :: a06=-0.01836710059356763_nr

 real(nr),parameter,dimension(3) :: fex=(/45,-9,1/)/(mfbi*30.0_nr)
 real(nr),parameter,dimension(0:5,0:2) :: abc=(/a01,a02,a03,a04,a05,a06,&
                                                a10,a12,a13,a14,a15,a16,&
                                                a20,a21,a23,a24,a25,a26/)

!===== ALLOCATABLE MAIN ARRAYS

 integer,dimension(:,:),allocatable :: npc,lio
 integer,dimension(:),allocatable :: li,lcsz,lctz,mxc
 integer,dimension(:),allocatable :: lxim,letm,lzem,lpos
 integer,dimension(:),allocatable :: lximb,letmb,lzemb,lhmb,mo

 real(nr),dimension(:,:),allocatable :: qo,qa,qb,de
 real(nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

 real(nr),dimension(:,:),allocatable :: xim,etm,zem
 real(nr),dimension(:),allocatable :: p,yaco

 real(nr),dimension(:,:),allocatable :: rr,ss

 real(nr),dimension(:,:),allocatable :: xu,yu
 real(nr),dimension(:,:),allocatable :: xl,yl
 real(nr),dimension(:),allocatable :: sa,sb

 real(nr),dimension(:,:),allocatable :: ran,sit,ait
 real(nr),dimension(:),allocatable :: xit,yit,zit
 real(nr),dimension(:),allocatable :: asz,bsz,csz,dsz
 real(nr),dimension(:),allocatable :: times,varmin,varmax,vmpi
 real(4),dimension(:),allocatable :: varr

 real(nr),dimension(:,:,:),pointer :: drva,drvb,send,recv,cm
 real(nr),dimension(:,:),pointer :: cmm

 real(nr),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
 real(nr),dimension(:,:,:),allocatable,target :: drvb1,drvb2,drvb3
 real(nr),dimension(:,:,:),allocatable,target :: send1,send2,send3
 real(nr),dimension(:,:,:),allocatable,target :: recv1,recv2,recv3
 real(nr),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3
 real(nr),dimension(:,:),allocatable,target :: cmm1,cmm2,cmm3

 character(13),dimension(:),allocatable :: cfilet
 character(7),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

 integer,dimension(0:1,0:1,3) :: ndf
 integer,dimension(3,3) :: ijk
 integer,dimension(0:1,3) :: nbc,ncd,nsz
 integer,dimension(0:4) :: no
 integer,dimension(-2:2) :: ilag
 integer,dimension(3) :: nbsize,nbcs,nbce,ncds,ncde,ms,me
 integer(8) :: lp,lq,ltomb
 integer :: lxio,leto,lzeo,lxi,let,lze,lmx,lim
 integer :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,kp,l,lh,ll
 integer :: m,ma,mb,mh,mm,mp,mbk,n,ndt,nn,nk,ns,ne,np,nt,nz,ndati,nsigi,nout,nfile
 integer :: nts,nscrn,nsgnl,ndata,nviscous,nkrk,nsmf,nfskp,nrestart

 real(nr),dimension(0:lmp,0:1,0:1) :: pbci,pbco
 real(nr),dimension(-2:2,0:2,0:1) :: albed,albef
 real(nr),dimension(mbci,mbci) :: cbca,cbcs
 real(nr),dimension(5,5) :: xt
 real(nr),dimension(0:1,0:1) :: pbcot
 real(nr),dimension(0:lmp) :: sap
 real(nr),dimension(mbci) :: rbci,sbci
 real(nr),dimension(5) :: cha,dha
 real(nr),dimension(-2:2) :: alag,blag,tlag
 real(nr),dimension(3) :: ve,dm,rv,uoo,umf,dudtmf
 real(nr),dimension(0:2) :: fam,fbm,fcm
 real(nr) :: alphf,betf,fa,fb,fc
 real(nr) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
 real(nr) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,fltk,fltkbc,dto
 real(nr) :: rhooo,poo,aoo,amachoo,srefoo,srefp1dre
 real(nr) :: dt,dts,dte,dtk,dtko,dtsum,timo,tsam,wts,wte,wtime
 real(nr) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

 character(1),dimension(0:4) :: cno
 character(3) :: cnzone,cndata
 character(5) :: cnnode
 character(7) :: czone
 character(17) :: coutput
 character(17) :: ctecout
 character(16) :: cinput
 character(16) :: cgrid
 character(18) :: cdata,cturb
 character(19) :: crestart,cpostdat

!===== INTEGER VARIABLES FOR MPI COMMANDS

 integer(4),dimension(:,:),allocatable :: ista
 integer(4),dimension(12) :: ireq
 integer(4) :: ir,mpro,npro,myid,itag,info,icom,ierr

!===== INTEGER VARIABLES FOR RECORDING BY RPT

 integer :: nrec,nwrec,nread,totVar

!===== VARIABLES FOR RECORDING BY RPT

 real(nr)  ::  xfor,yfor,rfor,amfor,tsfor,tefor
 integer   ::  lfor
 integer ,allocatable,dimension(:)  ::  lcfor
 real(nr),allocatable,dimension(:,:)  ::  xafor,yafor,bfor
 integer :: lcwall
 integer, dimension(:), allocatable ::lwall
 real(nr), dimension(:), allocatable ::area
 real(nr), dimension(:,:), allocatable ::wnor,wtan,tw
 real(nr), dimension(:,:), allocatable :: xyz
 logical :: wflag

!===== POST-PROCESSING VARIABLES BY RPT

 integer :: lsta
 logical :: tecplot,ispost
 real(nr), dimension(:,:), allocatable :: wplus
 real(nr), dimension(:), allocatable :: wvarr

 real(nr), dimension(2,2) :: cl
!=====

 end module mainvar3d

!*****
