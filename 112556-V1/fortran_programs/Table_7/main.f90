!  Code for: 
!  Maturity, Indebtedness and Default Risk    
!  Authors: 
!  Satyajit Chatterjee and Burcu Eyigungor
!
!  PROGRAM: Table 7
!  
!  The Program replicates Table 7
! 
! In this code, shock z corresponds to shock y in the paper, and 
! shock s corresponds to shock m in the paper.

program main

use globvar
use funcall

implicit none

integer, parameter::NN=6
real(doub),dimension(NN)::funcinp
real(doub):: res
integer::wild

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,id,ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)

call random_number(shockz)
call random_number(shockh)
call random_number(shocks)

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

!cdfz: cdf of z

cdfz(:,1)=pdfz(:,1)
do i=2,Z
cdfz(:,i)=cdfz(:,i-1)+pdfz(:,i)
end do
cdfz(:,Z)=1 

invpdfz=real(1)/real(Z)

do i=1,10000
do iz=1,Z
invpdfz2(iz)=dot_product(invpdfz,pdfz(:,iz))
end do
invpdfz=invpdfz2
end do

Ez=dot_product(zm,invpdfz)

! smthresh: the end points for the intervals of s
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

! pdfs: the probability of drawing s

do j=1,S    
call cdfnor(xx,x1,non,dble(smthresh(j+1)),dble(0),dble(sds),sta,boun)
call cdfnor(xx,x2,non,dble(smthresh(j)),dble(0),dble(sds),sta,boun)
pdfs(j)=x1-x2  
end do
divids=spread(sum(pdfs),DIM=1,NCOPIES=S)
pdfs=pdfs/divids


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
allocate(dar(nii*DA),daragg(nii*DA*nproc))
allocate(dar2(nii*S),daragg2(nii*S*nproc))
allocate(darint(nii*DA),darintagg(nii*DA*nproc))
allocate(valint(nii), valintagg(nii*nproc))

lambinvvec(1)=1.0
lambinvvec(2)=2.0
lambinvvec(3)=4.0
lambinvvec(4)=6.0
lambinvvec(5)=8.0
lambinvvec(6)=10.0
lambinvvec(7)=12.0
lambinvvec(8)=14.0
lambinvvec(9)=16.0
lambinvvec(10)=18.0
lambinvvec(11)=20.0

aminvec(1)=-1.35
aminvec(2)=-1.35
aminvec(3)=-1.3
aminvec(4)=-1.24
aminvec(5)=-1.2
aminvec(6)=-1.16
aminvec(7)=-1.13
aminvec(8)=-1.1
aminvec(9)=-1.08
aminvec(10)=-1.08
aminvec(11)=-1.08

do wild=1,LL
res=funcmin(wild)
end do

call mpi_finalize(ierr)

100 format (10000(1x, f12.5))
110 format (1000(1x,i4))
120 format (1000(1x, f7.3))


end program main
