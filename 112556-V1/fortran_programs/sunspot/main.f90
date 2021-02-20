!  Code for: 
!  Maturity, Indebtedness and Default Risk    
!  Authors: 
!  Satyajit Chatterjee and Burcu Eyigungor
!
!  Replicates Table 8 (sunspot equilibrium)
! 
!  DIRECT=1 Replicates Panel A of Table 8
!  lambda=0.05, solves for different sunspot probabilities
!
!  DIRECT=2 Replicates Panel B of Table 8
!  lambda=1.0, solves for different sunspot probabilities
!
!  In this code, shock z corresponds to shock y in the paper, and 
!  shock s corresponds to shock m in the paper.

program main

use globvar
use funcall

implicit none

integer:: DIRECT
integer, parameter::NN=4
real(doub),dimension(NN)::funcinp
real(doub):: res

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,id,ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)

DIRECT=1

call random_number(shockz)
call random_number(shockh)
call random_number(shocks)
call random_number(shockp)

! Creating transition matrix pdfz for the persistent shock (z)

lrsdz=sdz/sqrt(1-rhoz**2) !stdev of invariant distribution of log prod
delz=2*m/(Z-1)*lrsdz  !length of intervals between z in discretizations (equal distance)
zm(1)=muz-m*lrsdz

do i=2,Z
zm(i)=zm(i-1)+delz
end do

do i=1,Z
    do j=1,Z    
  call cdfnor(xx,x1,non,dble(zm(j)+delz/2),dble((1-rhoz)*muz+rhoz*zm(i)),dble(sdz),sta,boun)
  call cdfnor(xx,x2,non,dble(zm(j)-delz/2),dble((1-rhoz)*muz+rhoz*zm(i)),dble(sdz),sta,boun)
    pdfz(i,j)=x1-x2        
    end do
end do

dividz=spread(sum(pdfz,DIM=2),DIM=2,NCOPIES=Z)
pdfz=pdfz/dividz
logzm=zm
zm=exp(zm)

!cdfz: cdf of z.

cdfz(:,1)=pdfz(:,1)
do i=2,Z
cdfz(:,i)=cdfz(:,i-1)+pdfz(:,i)
end do
cdfz(:,Z)=1 

!smthresh: the end points for the intervals of s
dels=2*m2/S*sds
smthresh(1)=-m2*sds

do i=2,S+1
smthresh(i)=smthresh(i-1)+dels
end do

do i=1,S
sm(i)=(smthresh(i)+smthresh(i+1))*0.5
end do

do i=1,S
smdiff(i)=smthresh(i+1)-smthresh(i)
end do

!pdfs: the probability of drawing s

do j=1,S    
call cdfnor(xx,x1,non,dble(smthresh(j+1)),dble(0),dble(sds),sta,boun)
call cdfnor(xx,x2,non,dble(smthresh(j)),dble(0),dble(sds),sta,boun)
pdfs(j)=x1-x2        
end do
divids=spread(sum(pdfs),DIM=1,NCOPIES=S)
pdfs=pdfs/divids


! output shock z under default
do iz=1,Z
defz(iz)=(one-defp0-defp1*zm(iz))*zm(iz)
end do

dummy(:)=maxloc(defz)

if (dummy(1)<Z) then
defz(dummy(1):Z)=defz(dummy(1))
end if

! flow utility under default
do iz=1,Z
do is=1,S
num=sm(is)+defz(iz)
def_fix(is,iz)=num**(one-gamma)/(one-gamma)
end do
num=smthresh(1)+defz(iz)
def_fixn(iz)=num**(one-gamma)/(one-gamma)
end do

invpdfz=real(1)/real(Z)
do i=1,10000
do iz=1,Z
invpdfz2(iz)=dot_product(invpdfz,pdfz(:,iz))
end do
invpdfz=invpdfz2
end do
Ez=dot_product(zm,invpdfz)

sptvec(1)=0.00
sptvec(2)=0.005
sptvec(3)=0.01
sptvec(4)=0.02
sptvec(5)=0.05
sptvec(6)=0.1

Ez=dot_product(zm,invpdfz)

i1=0
do iz=1,Z
		do ia=1,A
			i1=i1+1			
			izloop(i1)=iz
			ialoop(i1)=ia
			iiv(iz,ia)=i1
end do
end do

nii=int(real(nza-1,8)/real(nproc,8))+1
itop=id*nii+1
iend=min((id+1)*nii,nza)

allocate(val1(nii), val2(nii), valagg(nii*nproc))
allocate(dar(nii*DA),dar2(nii*DA),daragg(nii*DA*nproc))
allocate(darint(nii*DA),darint2(nii*DA),darintagg(nii*DA*nproc))
allocate(valint(nii),valint2(nii), valintagg(nii*nproc))
allocate(dsv1(nii*S), dsagg(nii*S*nproc))

if (DIRECT==1) then

aminvec=-1.0
do kkk=1,LL
funcinp(1)=0.05
funcinp(2)=0.03
funcinp(3)=sptvec(kkk)
funcinp(4)=aminvec(kkk)
res=funcmin(funcinp)
end do

else if (DIRECT==2) then

aminvec=-1.3

do kkk=1,LL
funcinp(1)=1.0
funcinp(2)=0.0
funcinp(3)=sptvec(kkk)
funcinp(4)=aminvec(kkk)
res=funcmin(funcinp)
end do

end if

call mpi_finalize(ierr)

100 format (10000(1x, f12.8))
110 format (1000(1x,i4))
120 format (1000(1x, f7.3))

   
end program main
