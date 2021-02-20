 
module globvar
implicit none
include 'mpif.h'

integer, parameter::doub=8, sing=4 
integer, parameter::gamma=2
real(doub), parameter::sdz=0.027092 
real(doub), parameter::sds=0.003 
real(doub), parameter::muz=0 
real(doub), parameter::rhoz=0.948503
real(doub), parameter::coup=0.03
real(doub)::lambda
real(doub), parameter::rbase=0.01
real(doub), parameter::ent=0.0385
real(doub), parameter::beta=0.9540232420
real(doub)::qbase

integer, parameter::S=11, Z=200, A=350
real(doub), parameter::m=3.0, m2=2.0
real(doub), parameter::amax=0.0
real(doub)::amin
real(doub):: dela1, dela2

real(doub), dimension(A)::am
integer, dimension(A)::ambound
real(doub), dimension(S)::sm, smdiff
real(doub), dimension(S+1)::smthresh
real(doub), parameter::thrmin=-10, thrmax=10

real(doub):: lrsdz, delz, dels
real(doub), dimension(Z,Z)::pdfz, cdfz, dividz
real(doub),dimension(Z)::zm, logzm, temp
real(doub), dimension(S)::pdfs, divids
integer,parameter::xx=1
real(doub)::x1,x2, non,boun, x01, x02
integer:: sta
integer, parameter::time=3000 

real(doub), parameter::defp0=-0.1881927550, defp1=0.2455843389
integer:: srcmax
real(doub), dimension(S,Z):: def_fix
real(doub), dimension(Z):: def_fixn
real(doub), dimension(S,Z):: Vds
real(doub), dimension(Z):: vd_fut, Vd
real(doub), dimension(Z,A)::Vn_fut, zamvec, q, q2
real(doub),dimension(A)::fixnum
real(doub),dimension(Z)::zvec, zvec2
real(doub),dimension(S)::ucons
real(doub)::qpr, latval
real(doub), dimension(S)::q_fut
logical, dimension(A)::rela
real(doub), dimension(Z):: Vdsnone
real(doub), dimension(A):: defthr
real(doub), dimension(Z):: defz

integer, parameter:: DA=A
integer:: combint
real(doub), dimension(DA)::sthcomb, thcomb
integer, dimension(DA)::chass, asscomb
logical, dimension(DA)::chdef

integer, dimension(Z,A,DA)::assfin
real(doub), dimension(Z,A,DA)::thfin
integer, dimension(Z,A)::endint

integer:: t, i, j, j2, ij, k, int1, int2, int3
integer::is, iz, is2, iz2, ia, ia2, ia3, ia4
integer::kkk, ss, isind
logical:: defpr
real(doub)::num, num1, num2, num3, num4, Kons1, yy
integer, dimension(1:1)::dummy


real(doub),dimension(Z,A,A)::fixnumall

real(doub)::shockone
real(doub),dimension(A)::amser
real(doub),dimension(A,A)::amdiff
integer, parameter:: sim=300, TB=5000 
real(doub),dimension(TB, sim)::shockz, shockh, shocks
integer, dimension(TB)::zpath, apath, history
real(doub), dimension(TB):: spath, defpath, csim, ysim, nx, debts, qpath, assets
real(doub),dimension(TB-1):: sprpatha, rpath, spreadpath, spreadpatha
logical, dimension(TB-1)::samp1

real(doub),dimension(sim)::defaultpc, debtsvol
real(doub),dimension(sim)::meanspread, meandebt, meandebtserv

real(doub):: meanmeandef, meandebtsvol
real(doub):: meanmeanspr, meanmeandebt, meanmeansv
integer, parameter:: TBd=TB/4
real(doub),dimension(TBd)::debtssim
real(doub),dimension(4)::ds1, ds2, dsint

integer::aminloc

real, dimension(Z)::invpdfz, invpdfz2
real(doub)::Ez

integer,parameter:: nos=TB-1
real(doub), dimension(nos)::logy, nxy, logc, spreadt
real(doub),dimension(4, sim):: stdall1
real(doub),dimension(10,sim):: corrall1
real(doub), dimension(10):: corr1
real(doub), dimension(4)::stdn1

integer, dimension(TB)::inc, inc2
integer, parameter::maxinc=20


real(doub),parameter:: spreaddata=0.0815, debtdata=0.70, sprvoldata=0.0443

real(doub), parameter::one=1.0
integer, parameter::two=2
real(doub),parameter::eps=1e-15

integer,parameter:: idmaster=0
integer:: ierr,id,nproc
integer:: itop, iend, nii
integer::i0, i1, i2, i0pr
integer, parameter:: NZA=Z*A
integer, parameter:: NZAA=Z*A*A
integer, dimension(NZA):: izloop, ialoop
real(doub),dimension(A)::shorts
real(doub)::uexpt
integer::iins, iins2, iins3
integer, dimension(Z,A)::iiv
real(doub), dimension(NZA)::Vn, Vn2
real(doub), dimension(NZA)::qsum

real(doub):: errq, errv
real(doub), dimension(time)::errtot_v, errtot_q

real(doub),dimension(:),allocatable:: val1, val2, valagg, dar, dar2, daragg, daragg2
integer,dimension(:),allocatable:: valint, valintagg, darint, darintagg
integer, parameter::LL=11
real(sing), dimension(S,Z,A)::Vns
real(doub), dimension(LL)::lambinvvec
real(doub), dimension(LL)::aminvec
real(doub), dimension(LL,6)::funcreg

end module globvar
