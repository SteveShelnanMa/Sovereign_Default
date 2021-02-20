
% This program creates Figure 3 in Maturity, Indebtedness and Default Risk.
% 
% The following line loads the model output derived in the Fortran program
% in file xxx.

load input_figure3

% Deviation of log deseasonalized output from its linear trend
% Detrending/deseasonalizing is done using the output series between
% 1991Q1-2001Q4. The range: 1993Q1-2001Q4

y_dev=[-0.0114825065286013
    0.00191996320010190
    0.0128506379997937
    0.0197795310101867
    0.0368778236758232
    0.0449406444377978
    0.0384236033624674
    0.0510407369068844
    0.0336069659706091
    -0.0112465295180115
    -0.0255587983211960
    -0.0189751293788198
    0.00537009441107728
    0.0187543256414422
    0.0348071948770503
    0.0497100934167420
    0.0655113397863723
    0.0754476094526595
    0.0961495017920946
    0.103855801485729
    0.104370517685929
    0.120690280208986
    0.109153274015558
    0.0793871561683091
    0.0609529006360390
    0.0485696157747348
    0.0384264063942936
    0.0496140394739264
    0.0403661824601791
    0.0233410597087307
    0.0134915525765820
    0.00971206378687839
    0.000543111984008249
    0.00116508522358671
    -0.0561724150140819
    -0.121621314815158];


% Spread data from 1993Q1 to 2001Q4
dataspread=[9.311875
8.721841667
7.365241667
5.835208333
6.121908333
6.538837634
6.792622917
7.993606557
14.71641828
11.6785381
11.33520476
11.32854839
8.180645161
7.526349206
7.163560417
5.773353763
4.2315
3.736997917
3.002466667
4.574677419
4.38659071
4.492188889
7.610970833
7.480033333
7.685704372
6.860033333
8.090625
6.23309785
5.319047619
6.414124242
6.688503968
8.080638441
7.550907796
9.834818571
14.89385164
29.71842376];

dataspread=dataspread*0.01;
nn=size(y_dev,1);
rbase=0.01;
lambda=0.05;
coup=0.03;
A=350;
defpath=zeros(nn,1);
yy=exp(y_dev);

zind=zeros(nn,1);

for i=1:nn
    num=abs(zm-yy(i));
    zind(i)=find(num==min(num));
end

apath=zeros(nn+1,1);

apath(1)=106; 
zpath=zind;
spath=zeros(nn,1);
qbase=(lambda+(1-lambda)*coup)/(rbase+lambda);
qpath=zeros(nn,1);

for t=1:nn   
    
      zaind=(zpath(t)-1)*A+apath(t);
	  thcomb=thfin(zaind,:);

        for i=1:endint(zaind)
        if (spath(t)<=thcomb(i)) 
        int=i;
        break
        end 
        end 
	 
	  int3=assfin(zaind,int);

	  if (int3>0) 
		
  	  apath(t+1)=int3;
            qpath(t,1)= q(zpath(t),apath(t+1)); 
            defpath(t)=0;            
      else
            defpath(t)=1;
            qpath(t,1)= qbase;
            apath(t+1,1)=round(99);
      end 
end
        
        
rpath=zeros(nn,1);
spreadpatha=zeros(nn,1);

for i=1:nn-1
    
rpath(i)=(lambda+(1-lambda)*coup-lambda*qpath(i))/qpath(i);
spreadpatha(i)=(1+rpath(i)).^4-(1+rbase)^4;
end

xgraph=zeros(nn,1);
for i=1:nn
xgraph(i)=i;
end

 figure (1)
 plot(xgraph(1:nn-2),spreadpatha(1:nn-2),xgraph(1:nn-2),dataspread(1:nn-2))
 title('Simulated Spreads for Model and Data (1993Q1-2001Q3)')%,'fontsize',20)
 xlabel('Quarters')% ,'fontsize',20)
 ylabel('Spreads')%,'fontsize',20,'Rotation',0.0)
 legend('Model','Data')
