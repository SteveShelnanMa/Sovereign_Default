
%  Code for: 
%  Maturity, Indebtedness and Default Risk    
%  Authors:
%  Satyajit Chatterjee and Burcu Eyigungor
%  Replicates Figures 1-2-4-5

% The following line uploads the solution of the baseline program for 
% 50 grid points for z in order to make the non-concavitites more apparent
% in the figures. 

load inputs1245

lambdaL=0.05;
coupL=0.03;
A=350;

rbase=0.01;

zfix1=26;
for i=1:A
    budg(i)=-amL(i)*qL(zfix1,i);    
end

dispA1=90;
dispA2=150;

figure (1)
plot(amL(dispA1:dispA2),budg(dispA1:dispA2))
title('Figure 1: Revenue from Long-Term Debt')
xlabel('Debt Relative to Mean Output')
ylabel('(-q(y,bpr)bpr)')

gamma=2;
zfix1=20;
ca=90;

for i=1:A
    budga(i)=zm(zfix1)+lambdaL*amL(ca)+(1-lambdaL)*coupL*amL(ca)+((1-lambdaL)*amL(ca)-amL(i))*qL(zfix1,i);   
    if (budga(i)<0)
        dispA=i;
    end
    budgb(i)=budga(i)^(1-gamma)/(1-gamma)+Vn_futL(zfix1,i);      
end

dispA1=50;
dispA2=150;

figure (2)
plot(amL(dispA1:dispA2),budgb(dispA1:dispA2))
title('Figure 2: Variation over bpr in Total Lifetime Utility')
xlabel('Debt Relative to Mean Output')
ylabel('Lifetime Utility')

zfix1=26;
lambdaS=1;
coupS=0;

rS=zeros(A,1);
rL=zeros(A,1);
spreadS=zeros(A,1);
spreadL=zeros(A,1);
defproneS=zeros(A,1);
defproneL=zeros(A,1);

for i=1:A

rL(i)=(lambdaL+(1-lambdaL)*coupL-lambdaL*qL(zfix1,i))/qL(zfix1,i);
spreadL(i)=(1+rL(i))^4-(1+rbase)^4;
defproneS(i)=1-(1-defprobS(zfix1,i))^4;

rS(i)=(lambdaS+(1-lambdaS)*coupS-lambdaS*qS(zfix1,i))/qS(zfix1,i);
spreadS(i)=(1+rS(i))^4-(1+rbase)^4;
defproneL(i)=1-(1-defprobL(zfix1,i))^4;
 
end


i=A
while (spreadL(i)<0.3)
   i=i-1; 
end
dispA=i; 

figure (4)
plot(amL(dispA:A),spreadL(dispA:A), amL(dispA:A), defproneL(dispA:A))
title('Figure 4: Current Spreads and Next-Period Default Probability for Long-term Bonds')
xlabel('Debt Relative to Mean Output')
ylabel('Spreads and Default Probabilities')
legend('Spreads','Def Probability')

i=A;
while (spreadS(i)<0.3)
   i=i-1; 
end
dispA=i;    

figure (5)
plot(amS(dispA:A),spreadS(dispA:A), amS(dispA:A), defproneS(dispA:A))
title('Figure 5: Current Spreads and Next-Period Default Probability for Short-term Bonds')
xlabel('Debt Relative to Mean Output')
ylabel('Spreads and Default Probabilities')
legend('Spreads','Def Probability')





