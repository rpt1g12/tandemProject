!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

 module mainvar3d

 implicit none

!===== CONSTANT PARAMETERS

 integer,parameter :: nr=kind(0.0d0),ng0=0,ng1=1
 integer,parameter :: lmd=11,lmf=8,lmp=max(lmd,lmf)
 integer,parameter :: liofs=16,liofl=24

 character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

 real(nr),parameter :: pi=3.141592653589793_nr
 real(nr),parameter :: zero=0,one=1,half=0.5_nr,sqrt2=sqrt(2.0_nr),sqrt2i=1/sqrt2,sml=1.0e-6_nr
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

! real(nr),parameter :: beta20=0.04127253978047144_nr
! real(nr),parameter :: alpha21=0.4708395755079016_nr
! real(nr),parameter :: alpha23=0.5713690208719099_nr
! real(nr),parameter :: beta24=0.06287995158522702_nr
! real(nr),parameter :: a20=-0.1534532664885535_nr
! real(nr),parameter :: a21=-0.6866311200147498_nr
! real(nr),parameter :: a23=0.7176431952789228_nr
! real(nr),parameter :: a24=0.2186728528907302_nr
! real(nr),parameter :: a25=-0.001419994100359792_nr
! real(nr),parameter :: a26=0.0005236289985873258_nr
!
! real(nr),parameter :: alpha10=0.09486703622867607_nr
! real(nr),parameter :: alpha12=1.852980118858077_nr
! real(nr),parameter :: beta13=0.7841681122699989_nr
! real(nr),parameter :: a10=-0.3469447847494813_nr
! real(nr),parameter :: a12=0.1652135357932134_nr
! real(nr),parameter :: a13=1.379421330446014_nr
! real(nr),parameter :: a14=0.1789691155384915_nr
! real(nr),parameter :: a15=-0.02142195128295235_nr
! real(nr),parameter :: a16=0.001958948887672967_nr
!
! real(nr),parameter :: alpha01=5.590226531590711_nr
! real(nr),parameter :: beta02=3.911115464821060_nr
! real(nr),parameter :: a01=-3.320861355280472_nr
! real(nr),parameter :: a02=5.452259004221430_nr
! real(nr),parameter :: a03=1.150275611660523_nr
! real(nr),parameter :: a04=-0.1839611359673221_nr
! real(nr),parameter :: a05=0.05771607628595115_nr
! real(nr),parameter :: a06=-0.01431288821544212_nr

!===== ALLOCATABLE MAIN ARRAYS

 integer,dimension(:,:),allocatable :: npc,lio
 integer,dimension(:),allocatable :: li,lcsz,lctz,mxc
 integer,dimension(:),allocatable :: lxim,letm,lzem,lpos
 integer,dimension(:),allocatable :: lximb,letmb,lzemb,lhmb,mo

 real(nr),dimension(:,:),allocatable :: qo,qa,de
 real(nr),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

 real(nr),dimension(:,:),allocatable :: xim,etm,zem
 real(nr),dimension(:),allocatable :: p,yaco

 real(nr),dimension(:,:),allocatable :: rr,ss

 real(nr),dimension(:,:),allocatable :: xu,yu
 real(nr),dimension(:,:),allocatable :: xl,yl
 real(nr),dimension(:),allocatable :: sa,sb

 real(nr),dimension(:),allocatable :: asz,bsz,csz,dsz
 real(nr),dimension(:),allocatable :: ran,xit,yit,zit,sit,ait,bit,cit
 real(nr),dimension(:),allocatable :: times,varmin,varmax
 real(4),dimension(:),allocatable :: varr

 real(nr),dimension(:,:,:),pointer :: cm,send,recv,drva,drvb

 real(nr),dimension(:,:,:),allocatable,target :: drva1,drvb1
 real(nr),dimension(:,:,:),allocatable,target :: drva2,drvb2
 real(nr),dimension(:,:,:),allocatable,target :: drva3,drvb3

 real(nr),dimension(:,:,:),allocatable,target :: send1,recv1
 real(nr),dimension(:,:,:),allocatable,target :: send2,recv2
 real(nr),dimension(:,:,:),allocatable,target :: send3,recv3

 real(nr),dimension(:,:,:),allocatable,target :: cm1
 real(nr),dimension(:,:,:),allocatable,target :: cm2
 real(nr),dimension(:,:,:),allocatable,target :: cm3

!===== CONSTANT-SIZED MAIN VARIABLES

 integer,dimension(3,3) :: ijk
 integer,dimension(0:1,3) :: nbc,ncd,ndf,nsz
 integer,dimension(0:4) :: no
 integer,dimension(-2:2) :: ilag
 integer,dimension(3) :: nbsize,nbcs,nbce,ncds,ncde,ms,me
 integer(8) :: lp,lq,ltomb
 integer :: lxio,leto,lzeo,lxi,let,lze,lmx,lim
 integer :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,kp,l,lh,ll
 integer :: m,ma,mb,mh,mm,mp,mbk,n,nn,nn1,nn2,nn3,nk,ns,ne,np,nt,npg,ndati,nout,nfile
 integer :: nts,nscrn,ndata,ndatp,nviscous,nkrk,nsmf,nfbco,nfbcn,nfskp,nrestart

 real(nr),dimension(0:lmp,0:1,0:1) :: pbci,pbco
 real(nr),dimension(0:5,0:2) :: abc
 real(nr),dimension(0:8,0:1) :: albed,albef
 real(nr),dimension(5,5) :: xt
 real(nr),dimension(0:1,0:1) :: pbcot
 real(nr),dimension(0:4) :: fbc
 real(nr),dimension(0:lmp) :: sap
 real(nr),dimension(5) :: cha,dha
 real(nr),dimension(-2:2) :: alag,blag,tlag
 real(nr),dimension(3) :: ve,dm,rv,uoo,umf,dudtmf
 real(nr) :: alphf,betf,fa,fb,fc
 real(nr) :: alphf01,betf02,alphf10,alphf12,betf13
 real(nr) :: betf20,alphf21,alphf23,betf24,f20,f21,f23,f24,f25,fex1,fex2,fex3
 real(nr) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
 real(nr) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,fltk,fltkbc,dto
 real(nr) :: rhooo,poo,aoo,amachoo,srefoo,srefp1dre
 real(nr) :: dt,dtk,dtko,dtwi,timo,tsam,wts,wte,wtime
 real(nr) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

 character(1),dimension(0:4) :: cno
 character(7) :: czone
 character(13) :: coutput
 character(16) :: cinput
 character(16) :: cgrid
 character(18) :: cdata
 character(19) :: crestart

!===== INTEGER VARIABLES FOR MPI COMMANDS

 integer(4),dimension(:,:),allocatable :: ista
 integer(4),dimension(12) :: ireq
 integer(4) :: ir,mpro,npro,myid,itag,info,icom,ierr

!===== INTEGER VARIABLES FOR MPI COMMANDS by RPT

 integer :: nrec
 integer :: stnrec,etnrec
!=====
 end module mainvar3d

!*****
