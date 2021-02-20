module funcall
use globvar

implicit none

contains

real(doub) function funcmin (xx0)
implicit none
DOUBLE PRECISION dinvnr
EXTERNAL dinvnr

integer, intent(in)::xx0


if (xx0==1) then

beta=0.93696
defp0=-0.09774
defp1=0.15

else if (xx0==2) then

beta=0.95402
defp0=-0.18819
defp1=0.24558

else if (xx0==3) then

beta=0.96195
defp0=-0.28778
defp1=0.35
end if


! output shock y under default
do iz=1,Z
defz(iz)=(one-defp0-defp1*zm(iz))*zm(iz)
if (defz(iz)>zm(iz)) then
defz(iz)=zm(iz)
end if
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

! asset grids

dela1=(amax-amin)/real(A-1)
am(1)=amin
am(A)=amax
do i=2,A-1
am(i)=am(i-1)+dela1
end do

ambound=A

do iz=1,Z
do ia=1,A
zamvec(iz,ia)=zm(iz)+am(ia)*(lambda+(one-lambda)*coup)
end do
end do

! start with risk-free price
qbase=(lambda+(one-lambda)*coup)/(rbase+lambda)
q=qbase

errtot_q=0
errtot_v=0

do t=1,time

do iz=1,Z
Vd(iz)=dot_product(pdfs,Vds(:,iz))
end do

Vd_fut=beta*matmul(pdfz,Vd)

do ia=1,A
do iz=1,Z
i1=iiv(iz,ia)
zvec(iz)=Vn(i1)
end do
do iz=1,Z
Vn_fut(iz,ia)=beta*dot_product(pdfz(iz,:),zvec)
end do
end do

do iz=1,Z
num=ent*Vn_fut(iz,A)+(one-ent)*Vd_fut(iz) 
do is=1,S
Vds(is,iz)=def_fix(is,iz)+num 
end do
Vdsnone(iz)=def_fixn(iz)+num 
end do


do i1=itop,iend

shortS=thrmin
defthr=thrmin

iz=izloop(i1)
ia=ialoop(i1)

rela=.True.

do ia2=1,ambound(ia)
fixnum(ia2)=zamvec(iz,ia)-q(iz,ia2)*(am(ia2)-(one-lambda)*am(ia))
end do

! finding the highest asset choice that gives positive consumption

srcmax=1
do ia2=ambound(ia),1,-1
if ((fixnum(ia2)+smthresh(S+1))>0) then
srcmax=ia2
exit
endif
end do

! finding thresholds between adjacent asset grids where choice of asset changes

ia2=srcmax-1
do while (ia2>=1)

num1=fixnum(ia2)
num2=fixnum(ia2+1)
Kons1=-(Vn_fut(iz,ia2+1)-Vn_fut(iz,ia2)) ! Vn_fut is increasing in a' (ia2).

if (num2<num1) then ! if current consumption is higher for a higher debt level (or lower asset level) chosen, find the threshold s.
yy=(-(num1+num2)+sqrt((num1+num2)**two-real(4)*(num1*num2+(num1-num2)/min(Kons1,-eps))))/real(2)

shortS(ia2)=yy

else ! if not, then that asset level will never be chosen (rela(iz,ia,ia2)=0).
rela(ia2)=.false.
shortS(ia2)=thrmin
end if

ia2=ia2-1
end do

shortS(srcmax)=max(maxval(shorts(1:srcmax))+0.01,smthresh(S+1)+0.01)


! check which asset choices are dominated

ia2=srcmax-1 !Start with highest possible choice of asset

do while (ia2>1) 

if (rela(ia2-1)) then 

if (shortS(ia2-1)>shortS(ia2)) then ! if this holds then ia2 is dominated by ia2-1

rela(ia2)=.false.
shortS(ia2)=thrmin
defpr=.true.

if (ia2+1<=srcmax) then

do ia3=ia2+1,srcmax    ! we find out whether ia2-1 dominates other choices of a' (ia3) that are greater than ia2

if (defpr) then          ! as soon as we find the minimum ia3>ia2 that ia2-1 does not dominate, we end the loop
if (rela(ia3)) then 
! we find the threshold s between the choices of ia2-1 and ia3
num1=fixnum(ia2-1)
num2=fixnum(ia3)
Kons1=-(Vn_fut(iz,ia3)-Vn_fut(iz,ia2-1))

if (num2<num1) then
yy=(-(num1+num2)+sqrt((num1+num2)**two-real(4)*(num1*num2+(num1-num2)/min(Kons1,-eps))))/real(2)

if (yy>shortS(ia3).AND.ia3<srcmax) then  ! if this holds ia3 is dominated by ia2-1 also.
rela(ia3)=.false.
shortS(ia3)=thrmin

elseif (yy>shortS(ia3).AND.ia3==srcmax) then 
shortS(ia2-1)=yy            
shortS(srcmax)=yy+0.01

else  ! If ia2-1 does not dominate ia3 then we end the loop. 
shortS(ia2-1)=yy
defpr=.false.
end if


else ! In this case ia2-1 is dominated by ia3 (num2>num1, that is consumption is higher with the choice of higher asset level),
     ! we again end the loop.
rela(ia2-1)=.false.
shortS(ia2-1)=thrmin
defpr=.false.
end if

end if
end if
end do
end if
end if
end if
ia2=ia2-1

end do

! the asset choice thresholds that fall between smthresh(1) and smthresh(S+1)

sthcomb=0
i=0
chass=0
ia2=1

do while (ia2<=srcmax)
if (rela(ia2)) then
if (shortS(ia2)>=smthresh(1)) then
i=i+1
sthcomb(i)=shortS(ia2)
chass(i)=ia2
end if

if (shortS(ia2)>=smthresh(S+1)) then
sthcomb(i)=smthresh(S+1)
ia2=srcmax+1
end if
end if
ia2=ia2+1
end do

! combining asset thresholds with default threshold

combint=i
chdef=.False.
asscomb=0
thcomb=0

j=0
do i=1,combint
num=Vdsnone(iz)-Vn_fut(iz,chass(i))
if (num>=0) then
yy=thrmax
else
yy=((one-gamma)*num)**(one/(one-gamma))-fixnum(chass(i))
end if

if (i==1) then
num2=smthresh(1)
else
num2=sthcomb(i-1)
end if

if (yy<=num2) then
j=j+1
int1=combint-i
thcomb(j:j+int1)=sthcomb(i:combint)
chdef(j:j+int1)=.FALSE.
asscomb(j:j+int1)=chass(i:combint)
j=j+int1

exit

else if (yy>=sthcomb(i)) then
j=1
thcomb(j)=sthcomb(i)
chdef(j)=.TRUE.
asscomb(j)=0

else 
j=1
thcomb(j)=yy
chdef(j)=.TRUE.
asscomb(j)=0

j=j+1
int1=combint-i
thcomb(j:j+int1)=sthcomb(i:combint)
chdef(j:j+int1)=.FALSE.
asscomb(j:j+int1)=chass(i:combint)

j=j+int1

exit

end if

end do

combint=j

! integration over S

ucons=0
q_fut=0

isind=1
i=1

do while (i<=combint.AND.isind<=S) 
if (i==1) then
if (chdef(i)) then
isind=S
do j=2,S+1
if (thcomb(1)<smthresh(j)) then
isind=j-1
exit
end if
end do

ucons(isind)=ucons(isind)+((thcomb(1)-smthresh(isind))/smdiff(isind))*Vdsnone(iz)

if (isind>1) then
ucons(1:isind-1)=Vdsnone(iz)
end if

latval=thcomb(1)
i=i+1

else

if (thcomb(i)>=smthresh(isind+1)) then

num2=sm(isind)+fixnum(asscomb(i))
ucons(isind)=ucons(isind)+(num2**(1-gamma)/(1-gamma)+Vn_fut(iz,asscomb(i)))
q_fut(isind)=q_fut(isind)+(lambda+(one-lambda)*(coup+q(iz,asscomb(i))))

latval=smthresh(isind+1)
isind=isind+1

else

num=(thcomb(i)-smthresh(isind))/smdiff(isind)

num2=sm(isind)+fixnum(asscomb(i))
ucons(isind)=ucons(isind)+num*(num2**(1-gamma)/(1-gamma)+Vn_fut(iz,asscomb(i)))
q_fut(isind)=q_fut(isind)+num*(lambda+(one-lambda)*(coup+q(iz,asscomb(i))))

latval=thcomb(i)

i=i+1

end if

end if

else

if (thcomb(i)<smthresh(isind+1)) then

num=(thcomb(i)-latval)/smdiff(isind)
num2=sm(isind)+fixnum(asscomb(i))
ucons(isind)=ucons(isind)+num*(num2**(1-gamma)/(1-gamma)+Vn_fut(iz,asscomb(i)))
q_fut(isind)=q_fut(isind)+num*(lambda+(one-lambda)*(coup+q(iz,asscomb(i))))

latval=thcomb(i)
i=i+1

else  !if (thcomb(i)>smthresh(isind+1)) then

num=(smthresh(isind+1)-latval)/smdiff(isind)
num2=sm(isind)+fixnum(asscomb(i))
ucons(isind)=ucons(isind)+num*(num2**(1-gamma)/(1-gamma)+Vn_fut(iz,asscomb(i)))
q_fut(isind)=q_fut(isind)+num*(lambda+(one-lambda)*(coup+q(iz,asscomb(i))))

latval=smthresh(isind+1)
isind=isind+1

end if

end if

end do

qpr=dot_product(pdfs,q_fut)    
uexpt=dot_product(pdfs,ucons)

iins=i1-itop+1

val1(iins)=uexpt
val2(iins)=qpr
valint(iins)=combint

iins2=(iins-1)*DA
dar(iins2+1:iins2+DA)=thcomb
darint(iins2+1:iins2+DA)=asscomb

end do

! only price and value function will be called during solution of the model

call mpi_allgather(val1(1),nii,mpi_double_precision,&
valagg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
Vn2=valagg(1:nza)

call mpi_allgather(val2(1),nii,mpi_double_precision,&
valagg(1),nii,mpi_double_precision,mpi_comm_world,ierr)
qsum=valagg(1:nza)

do ia=1,A
do iz=1,Z
i1=iiv(iz,ia)
zvec(iz)=qsum(i1)
end do

do iz=1,Z
q2(iz,ia)=(one/(one+rbase))*dot_product(pdfz(iz,:),zvec)
end do

end do

errq=maxval(abs(q-q2))
errv=maxval(abs(Vn-Vn2))

num=0.5
q=num*q2+(one-num)*q
Vn=num*Vn2+(one-num)*Vn

errtot_q(t)=errq
errtot_v(t)=errv

if (t>250.AND.errq<0.00001 .AND. errv<0.00001) exit 

end do

! calling for other variables to simulate the model

call mpi_allgather(valint(1),nii,mpi_integer,&
valintagg(1),nii,mpi_integer,mpi_comm_world,ierr)

i1=0
do iz=1,Z
do ia=1,A
	i1=i1+1	
    endint(iz,ia)=valintagg(i1)
end do
end do

call mpi_allgather(dar(1),nii*DA,mpi_double_precision,&
daragg(1),nii*DA,mpi_double_precision,mpi_comm_world,ierr)

i1=0
do iz=1,Z
do ia=1,A
do ia2=1,DA
	i1=i1+1	
	thfin(iz,ia,ia2)=daragg(i1)
end do
end do
end do

call mpi_allgather(darint(1),nii*DA,mpi_integer,&
darintagg(1),nii*DA,mpi_integer,mpi_comm_world,ierr)

i1=0
do iz=1,Z
do ia=1,A
do ia2=1,DA
	i1=i1+1	
	assfin(iz,ia,ia2)=darintagg(i1)
end do
end do
end do

do iz=1,Z
do ia=1,A
do ia2=1,A
fixnumall(iz,ia,ia2)=zamvec(iz,ia)-q(iz,ia2)*(am(ia2)-(one-lambda)*am(ia))
end do
end do
end do

! debt service
do j=1,A
amser(j)=-am(j)*(lambda+(1-lambda)*coup)
end do

! debt buyback
do i=1,A
do j=1,A
amdiff(i,j)=am(j)-(one-lambda)*am(i)
end do
end do


! simulating the economy
zpath(TB)=Z/2 ! intiate zpath

do ss=1,sim 

! find the z-path using transition matrix
zpath(1)=zpath(TB)

do t=2,TB
temp=cdfz(zpath(t-1),:)
				do i=1,Z 
				if (shockz(t,ss).le.temp(i)) then
				zpath(t)=i
				exit
 				end if
				end do				
end do

! finding s-path
do i=1,TB
num=shocks(i,ss) 
spath(i)=dinvnr(shocks(i,ss),one-shocks(i,ss))*sds
if (spath(i)>smthresh(S+1).OR.spath(i)<smthresh(1))then
defpr=.true.
do while (defpr)
call random_number(shockone)
spath(i)=dinvnr(shockone,one-shockone)*sds
if (spath(i)<=smthresh(S+1).AND.spath(i)>=smthresh(1))then
defpr=.false.
end if
end do
end if
end do

apath(1)=A !initiating assets
history(1)=1 ! history: 1 if included in financial markets, 0 if excluded from financial markets.

do t=1,TB-1
        if (history(t)==1.OR.(history(t)==0.AND.shockh(t,ss)<ent)) then    ! if included in fin markets or redeemed into fin markets.
        
        ! find the asset level it will want to hold
        thcomb=thfin(zpath(t),apath(t),:)
        do i=1,endint(zpath(t),apath(t))
        if (spath(t)<=thcomb(i)) then 
        int1=i
        exit
        end if
        end do
	  
        if (spath(t)>thcomb(endint(zpath(t),apath(t)))) then
        int1=endint(zpath(t),apath(t))       
        end if

	  int3=assfin(zpath(t),apath(t),int1)

	 if (int3>0) then
		
  	    apath(t+1)=int3
        csim(t)=fixnumall(zpath(t),apath(t),apath(t+1))+spath(t)  ! consumption
        ysim(t)=zm(zpath(t))+spath(t)  ! output
        nx(t)=ysim(t)-csim(t)   ! exports
        qpath(t)=q(zpath(t),apath(t+1)) ! price of bonds
        defpath(t)=0            
        history(t+1)=1          
        history(t)=1
        
        if (amdiff(apath(t), apath(t+1))>0) then  ! if there is buyback of debt
        debts(t)=amser(apath(t))+qpath(t)*amdiff(apath(t), apath(t+1)) ! debt-service
        else
        debts(t)=amser(apath(t))
        end if     

        if (t==1) then
        inc(t)=1	
        inc2(t)=1
        else
        inc(t)=inc(t-1)+1
        inc2(t)=inc2(t-1)+1
        end if

	else
		
		inc(t)=0		   
		if (t==1) then        
        inc2(t)=1
        else
        inc2(t)=inc2(t-1)+1
        end if            
        ysim(t)=zm(zpath(t))+spath(t)
        csim(t)=ysim(t)
        apath(t+1)=A  ! starting with zero debt level
        qpath(t)=qbase
        defpath(t)=1           
        history(t+1)=0
        history(t)=1  
        debts(t)=0                    
        nx(t)=0    
    end if

 else 
 inc2(t)=0
 inc(t)=0
 ysim(t)=zm(zpath(t))+spath(t)
 csim(t)=ysim(t)
 defpath(t)=0
 history(t+1)=0            
 apath(t+1)=A        
 qpath(t)=qbase
 debts(t)=0
 nx(t)=0
 end if       
end do

rpath=(lambda+(one-lambda)*coup-lambda*qpath(1:TB-1))/qpath(1:TB-1) !interest rate
spreadpath=rpath-rbase 
spreadpatha=(1+rpath)**4-(1+rbase)**4 ! annual spreads
assets=am(apath)  ! outstanding obligations

dummy(:)=minloc(apath(1:TB))  ! finding minimum asset choice to see if it is binding
if (ss==1) then
aminloc=apath(dummy(1))
else if (aminloc>apath(dummy(1))) then
aminloc=apath(dummy(1))
end if

! creating the logical samp1 to exclude the first maxinc years after reentry from simulations
num=0
num2=0
do i=1,TB-1
if (inc(i)>maxinc) then 
num=num+1   ! size of sample
samp1(i)=.TRUE.
else
samp1(i)=.FALSE.
end if

if (inc2(i)>maxinc) then 
num2=num2+1 ! size of sample to calculate default probability
end if
end do

meanspread(ss)=sum(spreadpatha, MASK=samp1)/num
meandebt(ss)=sum(assets(2:TB)/ysim(1:TB-1), MASK=samp1)/num
num3=sum(defpath(1:TB-1),MASK=inc2(1:TB-1)>maxinc)/num2 ! quarterly default probability
defaultpc(ss)=one-(one-num3)**4  ! yearly default probability

!turning quarterly debt service data to yearly
t=0
j=0
do i=1,TB-1
j=j+1
ds1(j)=debts(i)
ds2(j)=ysim(i)
if (samp1(i)) then
dsint(j)=1
else
dsint(j)=0
end if
if (j==4) then
j=0
if (sum(dsint)>0) then
t=t+1
debtssim(t)=dot_product(dsint,ds1)/dot_product(dsint,ds2)
end if
end if
end do

meandebtserv(ss)=sum(debtssim(1:t))/real(t)
debtsvol(ss)=sqrt(sum((debtssim(1:t)-meandebtserv(ss))**2)/real(t))

logy=log(ysim(TB-nos:TB-1))
nxy = nx(TB-nos:TB-1)/ysim(TB-nos:TB-1)
logc = log(csim(TB-nos:TB-1))
spreadt=spreadpatha(TB-nos:TB-1) 

call moment(logy, stdall1(1,ss))
call moment(logc, stdall1(2,ss))
call moment(spreadt, stdall1(3,ss))
call moment(nxy, stdall1(4,ss))

corrall1(1,ss)=corr(logy, logc)
corrall1(2,ss)=corr(logy, spreadt)
corrall1(3,ss)=corr(logy, nxy)
corrall1(4,ss)=corr(logc, spreadt)
corrall1(5,ss)=corr(logc, nxy)
corrall1(6,ss)=corr(spreadt, nxy)

end do

corr1=sum(corrall1,DIM=2)/real(sim)
stdn1=sum(stdall1,DIM=2)/real(sim)

meanmeandef=sum(defaultpc)/real(sim)
meanmeanspr=sum(meanspread)/real(sim)
meanmeandebt=sum(meandebt)/real(sim)
meanmeansv=sum(meandebtserv)/real(sim)
meandebtsvol=sum(debtsvol)/real(sim)

funcreg(xx0,1)=defp1
funcreg(xx0,2)=defp0
funcreg(xx0,3)=beta
funcreg(xx0,4)=stdn1(3)
funcreg(xx0,5)=meanmeanspr
funcreg(xx0,6)=-meanmeandebt

!!!!!!!!!!!!!!!!!Results!!!!!!!!!!!!!!!!!!!!

open(20,file='Table 5.out') 

Write(20,'(10(3X,A))') 'd1', 'd0', 'beta', 'std. dev. spread', 'avg. spread', 'avg. b/y'
do i=1,xx0
Write(20,'(10(2X,f8.4))') (funcreg(i,j),j=1,6)
end do
close(20) 

funcmin=100

end function funcmin

subroutine moment(xx, std)
implicit none
integer::mm
real(doub),dimension(:), intent(in)::xx
real(doub), intent(out):: std
real(doub):: mean

mm=sum(history(TB-nos:TB-1),MASK=samp1(TB-nos:TB-1))
mean=sum(xx,MASK=samp1(TB-nos:TB-1))/mm
std=sum((xx-mean)**2,MASK=samp1(TB-nos:TB-1))/mm
std=sqrt(std)

return 
end subroutine moment

real(doub) function corr(xx1,xx2)
implicit none
integer::mm
real(doub)::mean1, mean2, std1, std2
real(doub), dimension(:), intent(in)::xx1, xx2

mm=sum(history(TB-nos:TB-1),MASK=samp1(TB-nos:TB-1))
mean1=sum(xx1,MASK=samp1(TB-nos:TB-1))/mm
mean2=sum(xx2,MASK=samp1(TB-nos:TB-1))/mm
std1=sqrt(sum((xx1-mean1)**2,MASK=samp1(TB-nos:TB-1))/mm)
std2=sqrt(sum((xx2-mean2)**2,MASK=samp1(TB-nos:TB-1))/mm)
corr=sum((xx1-mean1)*(xx2-mean2),MASK=samp1(TB-nos:TB-1))/mm
corr=corr/(std1*std2)

return
end function corr

end module funcall