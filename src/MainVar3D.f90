!*****
!***** MAIN VARIABLES & DATA FOR 3D NAVIER-STOKES/EULER SOLVER
!*****

 module mainvar3d

 implicit none

!===== CONSTANT PARAMETERS

 integer,parameter :: k4=kind(1.0e0),k8=kind(1.0d0)
 integer(k4),parameter :: ntdrv=0,ntflt=1,nrall=0,nrone=1,n45no=0,n45go=1
 integer(k4),parameter :: lmd=11,lmf=8,lmp=max(lmd,lmf),mfbi=2,mbci=3
 integer(k4),parameter :: liofs=16,liofl=24

 character(len=*),parameter :: fmts='es15.8',fmtl='es23.16',fmtsa=fmts//',a',fmtla=fmtl//',a'

 real(k8),parameter :: pi=3.141592653589793_k8
 real(k8),parameter :: zero=0,one=1,half=0.5_k8,sqrt2=sqrt(2.0_k8),sqrt2i=1/sqrt2,sml=1.0e-6_k8,free=1.0d+6
 real(k8),parameter :: two=2,quarter=0.25_k8
 real(k8),parameter :: gam=1.4_k8,gamm1=gam-1,ham=1/gam,hamm1=1/gamm1
 real(k8),parameter :: prndtl=0.71_k8,gamm1prndtli=1/(gamm1*prndtl)
 real(k8),parameter :: tprndtl=0.99_k8,tgamm1prndtli=1/(gamm1*tprndtl)

 real(k8),parameter :: alpha=0.5862704032801503_k8
 real(k8),parameter :: beta=0.09549533555017055_k8
 real(k8),parameter :: aa=0.6431406736919156_k8
 real(k8),parameter :: ab=0.2586011023495066_k8
 real(k8),parameter :: ac=0.007140953479797375_k8

 real(k8),parameter :: beta20=0.03250008295108466_k8
 real(k8),parameter :: alpha21=0.3998040493524358_k8
 real(k8),parameter :: alpha23=0.7719261277615860_k8
 real(k8),parameter :: beta24=0.1626635931256900_k8
 real(k8),parameter :: a20=-0.1219006056449124_k8
 real(k8),parameter :: a21=-0.6301651351188667_k8
 real(k8),parameter :: a23=0.6521195063966084_k8
 real(k8),parameter :: a24=0.3938843551210350_k8
 real(k8),parameter :: a25=0.01904944407973912_k8
 real(k8),parameter :: a26=-0.001027260523947668_k8

 real(k8),parameter :: alpha10=0.08360703307833438_k8
 real(k8),parameter :: alpha12=2.058102869495757_k8
 real(k8),parameter :: beta13=0.9704052014790193_k8
 real(k8),parameter :: a10=-0.3177447290722621_k8
 real(k8),parameter :: a12=-0.02807631929593225_k8
 real(k8),parameter :: a13=1.593461635747659_k8
 real(k8),parameter :: a14=0.2533027046976367_k8
 real(k8),parameter :: a15=-0.03619652460174756_k8
 real(k8),parameter :: a16=0.004080281419108407_k8

 real(k8),parameter :: alpha01=5.912678614078549_k8
 real(k8),parameter :: beta02=3.775623951744012_k8
 real(k8),parameter :: a01=-3.456878182643609_k8
 real(k8),parameter :: a02=5.839043358834730_k8
 real(k8),parameter :: a03=1.015886726041007_k8
 real(k8),parameter :: a04=-0.2246526470654333_k8
 real(k8),parameter :: a05=0.08564940889936562_k8
 real(k8),parameter :: a06=-0.01836710059356763_k8

 real(k8),parameter,dimension(3) :: fex=(/45,-9,1/)/(mfbi*30.0_k8)

!===== ALLOCATABLE MAIN ARRAYS

 integer(k4),dimension(:,:),allocatable :: npc,lio
 integer(k4),dimension(:),allocatable :: li,lcsz,lctz,mxc
 integer(k4),dimension(:),allocatable :: lxim,letm,lzem,lpos
 integer(k4),dimension(:),allocatable :: lximb,letmb,lzemb,lhmb,mo

 real(k8),dimension(:,:),allocatable :: qo,qa,qb,de
 real(k8),dimension(:),allocatable :: txx,tyy,tzz,txy,tyz,tzx,hxx,hyy,hzz

 real(k8),dimension(:,:),allocatable :: xim,etm,zem
 real(k8),dimension(:),allocatable :: p,yaco
 real(k8),dimension(:),allocatable :: sbcc

 real(k8),dimension(:,:),allocatable :: rr,ss

 real(k8),dimension(:,:),allocatable :: xu,yu
 real(k8),dimension(:,:),allocatable :: xl,yl
 real(k8),dimension(:),allocatable :: sa,sb

 real(k8),dimension(:,:),allocatable :: ran,sit,ait
 real(k8),dimension(:),allocatable :: xit,yit,zit
 real(k8),dimension(:),allocatable :: asz,bsz,csz
 real(k8),dimension(:),allocatable :: times,vmpi
 real(4),dimension(:),allocatable :: varr

 real(k8),dimension(:,:,:),pointer :: drva,drvb,send,recv,cm

 real(k8),dimension(:,:,:),allocatable,target :: drva1,drva2,drva3
 real(k8),dimension(:,:,:),allocatable,target :: drvb1,drvb2,drvb3
 real(k8),dimension(:,:,:),allocatable,target :: send1,send2,send3
 real(k8),dimension(:,:,:),allocatable,target :: recv1,recv2,recv3
 real(k8),dimension(:,:,:),allocatable,target :: cm1,cm2,cm3

 real(k8),dimension(:,:),allocatable :: varm
 real(k8),dimension(:),allocatable :: varmin,varmax
 real(k4),dimension(:),allocatable :: vart,vara,varb,vmean

 character(13),dimension(:),allocatable :: ctecplt,cthead
 character(13),dimension(:),allocatable :: cfilet
 character(7),dimension(:),allocatable :: czonet

!===== CONSTANT-SIZED MAIN VARIABLES

 integer(k4),dimension(0:1,0:1,3) :: ndf
 integer(k4),dimension(3,3) :: ijk
 integer(k4),dimension(0:1,3) :: nbc,ncd,nsz
 integer(k4),dimension(0:4) :: no
 integer(k4),dimension(-2:2) :: ilag
 integer(k4),dimension(3) :: nbsize,nbcs,nbce,ncds,ncde,ms,me
 integer(8) :: lp,lq,ltomb
 integer(k4) :: lxio,leto,lzeo,lxi,let,lze,lmx,lim,nrec
 integer(k4) :: i,ii,is,ie,ip,iq,j,jj,js,je,jp,jq,jk,k,kk,kp,l,lh,ll
 integer(k4) :: m,ma,mb,mm,mp,mq,mbk,mps,mpe,n,ndt,nn,nk,ns,ne,np,nt,nz,ndati,nsigi,nout,nfile
 integer(k4) :: nts,nscrn,nsgnl,ndata,nviscous,nkrk,nsmf,nfskp,nrestart
 integer(k8) :: nlmx,llmb,llmo,lis,lie,ljs,lje

 real(k8),dimension(0:lmp,0:1,0:1) :: pbci,pbco
 real(k8),dimension(-2:2,0:2,0:1) :: albed,albef
 real(k8),dimension(mbci,mbci) :: cbca,cbcs
 real(k8),dimension(5,5) :: xt
 real(k8),dimension(0:5,0:2) :: abc
 real(k8),dimension(0:1,0:1) :: pbcot
 real(k8),dimension(0:lmp) :: sap
 real(k8),dimension(mbci) :: rbci,sbci
 real(k8),dimension(5) :: cha,dha
 real(k8),dimension(-2:2) :: alag,blag,tlag
 real(k8),dimension(3) :: ve,dm,rv,uoo,umf,dudtmf
 real(k8),dimension(0:2) :: fam,fbm,fcm
 real(k8) :: alphf,betf,fa,fb,fc
 real(k8) :: ra0,ra1,ra2,ra3,res,fctr,dfdt
 real(k8) :: reoo,tempoo,amach1,amach2,amach3,wtemp,cfl,tmax,timf,fltk,fltkbc,dto
 real(k8) :: rhooo,poo,aoo,amachoo,srefoo,srefp1dre
 real(k8) :: dt,dts,dte,dtk,dtko,dtsum,timo,tsam,wts,wte,wtime
 real(k8) :: vn,vs,hv2,ao,bo,co,ho,aoi,rhoi,progmf,sqrtrema,sqrtremai

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

 integer(k4),dimension(:,:),allocatable :: ista
 integer(k4),dimension(:),allocatable :: ireq
 integer(k4) :: ir,mpro,npro,myid,itag,info,icom,ierr

!===== INTEGER VARIABLES FOR RECORDING BY RPT

 integer(k4) :: ngrec,nwrec,nread,totVar

!===== VARIABLES FOR FORCING BY RPT

 real(k8)  ::  xfor,yfor,rfor,amfor,tsfor,tefor
 integer(k4)   ::  lfor
 integer(k4) ,allocatable,dimension(:)  ::  lcfor
 real(k8),allocatable,dimension(:,:)  ::  xafor,yafor,bfor

!===== VARIABLES FOR WALL OPERATIONS BY RPT
 integer(k4) :: lcwall
 integer(k4), dimension(:), allocatable ::lwall
 real(k8), dimension(:), allocatable ::area
 real(k8), dimension(:,:), allocatable ::wnor,wtan,tw
 logical :: wflag
 real(k8), dimension(:,:), allocatable,target :: xyz

!===== POST-PROCESSING VARIABLES BY RPT

 integer(k4) :: lsta
 logical :: tecplot,ispost
 real(k8), dimension(:,:), allocatable :: wplus
 real(k8), dimension(:), allocatable :: wvarr
 character(18),dimension(:),allocatable :: ofiles 

 real(k8), dimension(2,2) :: cl

!===== ADITIONAL INPUTO VARIABLES BY RPT
 integer(k4)  :: nto,iwrec
 integer(k4)  :: forcing,LES
 real(k8) :: tgustd,tguste
 real(k8) :: talphas,talphar,aoa
 real(k8) :: smago1,smago2
 integer(k4)  :: output,ogrid,osol,oblock

 integer(k4) :: bkx,bky,bkz
 integer(k4), dimension(:),allocatable:: lxibk,letbk,lzebk
!=====
contains
function int4(x) result(res)
integer :: x
integer(k4) :: res
res=int(x,k4)
end function int4


 end module mainvar3d

!*****
