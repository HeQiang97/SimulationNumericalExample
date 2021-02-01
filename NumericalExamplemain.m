close all;
clear;
clc;

%% The Simulation of Numerical Example in Section IV %%

%% System dynamics
A=[0 1 0;-1 0 2;0 0 -1]; B=[1; 1; 1]; C=[1 0 0];

%% communication topology
huaA=[0 1 0 1 0 0;1 0 0 0 -1 -1;0 0 0 1 -1 0;1 0 1 0 0 0;0 -1 -1 0 0 0;0 -1 0 0 0 0];
huaL=[2 -1 0 -1 0 0;-1 3 0 0 1 1;0 0 2 -1 1 0;-1 0 -1 2 0 0;0 1 1 0 2 0;0 1 0 0 0 1];
huaB=diag([1 0 0 1 0 0]);
huaH=[3 -1 0 -1 0 0;-1 3 0 0 1 1;0 0 2 -1 1 0;-1 0 -1 3 0 0;0 1 1 0 2 0;0 1 0 0 0 1];
a10=1; a20=0; a30=0; a40=1; a50=0; a60=0;
S=diag([1 1 1 1 -1 -1]);
s1=1;s2=1;s3=1;s4=1;s5=-1;s6=-1;

%% Ni£ºthe number of neighbors for the ith agent.
N1=2; N2=3; N3=2; N4=2; N5=2; N6=1;

%% select F and L such that Assumption 2 holds.
F=[-5, -5, -6]; eig(A+B*F)
L=[-5 -2 -1]'; eig(A+L*C)

%% The parameters of two triggering mechanisms (7) (11) and controller (14) are given
k=0.1;  n=3; N=6;
P=are(A, B*B', k*eye(n));
K=B'*P;
k1=0.05; k2=0.05; k3=0.05; k4=0.05; k5=0.05; k6=0.05;
beta=2; rrou1=2; rrou2=2; gama=2; 
rou1=(N1+a10)^2+(N1+a10)*N1+N*N1+(N-1)*N1;
rou2=(N2+a20)^2+(N2+a20)*N2+N*N2+(N-1)*N2;
rou3=(N3+a30)^2+(N3+a30)*N3+N*N3+(N-1)*N3;
rou4=(N4+a10)^2+(N4+a10)*N4+N*N4+(N-1)*N4;
rou5=(N5+a10)^2+(N5+a10)*N5+N*N5+(N-1)*N5;
rou6=(N6+a10)^2+(N6+a10)*N6+N*N6+(N-1)*N6;
detaphi=10; detatheta=10; mu1=5; v1=0.2; mu2=5; v2=0.2;
 
tmax=80; %% maximum simulation time (s)
T=0.01; %% Iteration step (s)
kmax=1+tmax/T; %% maximum number of iterations
acmax=30*100+1; %% Actuator failure time

[Num Num]=size(huaL);
[Ord_A Ord_A]=size(A);
[Ord_B1 Ord_B2]=size(B);
[Ord_C1 Ord_C2]=size(C);

%% All state variables are set to 0
%% x_{i}(t) in (1) and (2)
x1=zeros(Ord_A*kmax,1);
x2=zeros(Ord_A*kmax,1);
x3=zeros(Ord_A*kmax,1);
x4=zeros(Ord_A*kmax,1);
x5=zeros(Ord_A*kmax,1);
x6=zeros(Ord_A*kmax,1);
x0=zeros(Ord_A*kmax,1);

%% x_{0}(T_{k_{i}}^{i}) in (14)
x01last=zeros(Ord_A*kmax,1);
x02last=zeros(Ord_A*kmax,1);
x03last=zeros(Ord_A*kmax,1);
x04last=zeros(Ord_A*kmax,1);
x05last=zeros(Ord_A*kmax,1);
x06last=zeros(Ord_A*kmax,1);

%% output
y1=zeros(Ord_C1*kmax,1);
y2=zeros(Ord_C1*kmax,1);
y3=zeros(Ord_C1*kmax,1);
y4=zeros(Ord_C1*kmax,1);
y5=zeros(Ord_C1*kmax,1);
y6=zeros(Ord_C1*kmax,1);
y0=zeros(Ord_C1*kmax,1);

%% controller
u1=zeros(Ord_B2*kmax,1);
u2=zeros(Ord_B2*kmax,1);
u3=zeros(Ord_B2*kmax,1);
u4=zeros(Ord_B2*kmax,1);
u5=zeros(Ord_B2*kmax,1);
u6=zeros(Ord_B2*kmax,1);

%% real controller
uf1=zeros(Ord_B2*kmax,1);
uf2=zeros(Ord_B2*kmax,1);
uf3=zeros(Ord_B2*kmax,1);
uf4=zeros(Ord_B2*kmax,1);
uf5=zeros(Ord_B2*kmax,1);
uf6=zeros(Ord_B2*kmax,1);

%% multiplicative factors in (3)
Theta1=zeros(kmax,1);
Theta2=zeros(kmax,1);
Theta3=zeros(kmax,1);
Theta4=zeros(kmax,1);
Theta5=zeros(kmax,1);
Theta6=zeros(kmax,1);

%% additive factors in (3)
phi1=zeros(Ord_B2*kmax,1);
phi2=zeros(Ord_B2*kmax,1);
phi3=zeros(Ord_B2*kmax,1);
phi4=zeros(Ord_B2*kmax,1);
phi5=zeros(Ord_B2*kmax,1);
phi6=zeros(Ord_B2*kmax,1);

%% state observer \eta_{i}(t) in (6)
eita1=zeros(Ord_A*kmax,1);
eita2=zeros(Ord_A*kmax,1);
eita3=zeros(Ord_A*kmax,1);
eita4=zeros(Ord_A*kmax,1);
eita5=zeros(Ord_A*kmax,1);
eita6=zeros(Ord_A*kmax,1);

%% \eta_{i}(t_{\sigma_{i}}^{i}) in (9)
eita1last=zeros(Ord_A*kmax,1);
eita2last=zeros(Ord_A*kmax,1);
eita3last=zeros(Ord_A*kmax,1);
eita4last=zeros(Ord_A*kmax,1);
eita5last=zeros(Ord_A*kmax,1);
eita6last=zeros(Ord_A*kmax,1);

%% \eta_{i}(t_{\sigma_{i}}^{i}) in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
eita1last0=zeros(Ord_A*kmax,1);
eita2last0=zeros(Ord_A*kmax,1);
eita3last0=zeros(Ord_A*kmax,1);
eita4last0=zeros(Ord_A*kmax,1);
eita5last0=zeros(Ord_A*kmax,1);
eita6last0=zeros(Ord_A*kmax,1);

%% \hat{\zeta}_{i}(T_{k_{ij}}^{ij}) in (14)
zetaj1last=zeros(Ord_A*kmax,1);
zetaj2last=zeros(Ord_A*kmax,1);
zetaj3last=zeros(Ord_A*kmax,1);
zetaj4last=zeros(Ord_A*kmax,1);
zetaj5last=zeros(Ord_A*kmax,1);
zetaj6last=zeros(Ord_A*kmax,1);

%% c_{i}(t) in (16)
c1=zeros(kmax,1);
c2=zeros(kmax,1);
c3=zeros(kmax,1);
c4=zeros(kmax,1);
c5=zeros(kmax,1);
c6=zeros(kmax,1);

%% c_{i}(T_{k_{i}}^{i}) in (14)
c1last=zeros(kmax,1);
c2last=zeros(kmax,1);
c3last=zeros(kmax,1);
c4last=zeros(kmax,1);
c5last=zeros(kmax,1);
c6last=zeros(kmax,1);

%% e_{i}(t) in (9)
e1=zeros(Ord_A*kmax,1);
e2=zeros(Ord_A*kmax,1);
e3=zeros(Ord_A*kmax,1);
e4=zeros(Ord_A*kmax,1);
e5=zeros(Ord_A*kmax,1);
e6=zeros(Ord_A*kmax,1);

%% \epsilon_{i}(t) in (13)
ee1=zeros(Ord_A*kmax,1);
ee2=zeros(Ord_A*kmax,1);
ee3=zeros(Ord_A*kmax,1);
ee4=zeros(Ord_A*kmax,1);
ee5=zeros(Ord_A*kmax,1);
ee6=zeros(Ord_A*kmax,1);

%% The first two terms of formula (8)
err11=zeros(kmax,1);
err12=zeros(kmax,1);
err13=zeros(kmax,1);
err14=zeros(kmax,1);
err15=zeros(kmax,1);
err16=zeros(kmax,1);

%% The first two terms of formula (12)
err21=zeros(kmax,1);
err22=zeros(kmax,1);
err23=zeros(kmax,1);
err24=zeros(kmax,1);
err25=zeros(kmax,1);
err26=zeros(kmax,1);

%% \mu_{1 i} e^{-v_{1 i} t} in (8)
thre11=zeros(kmax,1);
thre12=zeros(kmax,1);
thre13=zeros(kmax,1);
thre14=zeros(kmax,1);
thre15=zeros(kmax,1);
thre16=zeros(kmax,1);

%% \mu_{2 i} e^{-v_{2 i} t} in (12)
thre21=zeros(kmax,1);
thre22=zeros(kmax,1);
thre23=zeros(kmax,1);
thre24=zeros(kmax,1);
thre25=zeros(kmax,1);
thre26=zeros(kmax,1);

%% f_{1 i}(t) in (8)
ff11=zeros(kmax,1);
ff12=zeros(kmax,1);
ff13=zeros(kmax,1);
ff14=zeros(kmax,1);
ff15=zeros(kmax,1);
ff16=zeros(kmax,1);

%% f_{2 i}(t) in (12)
ff21=zeros(kmax,1);
ff22=zeros(kmax,1);
ff23=zeros(kmax,1);
ff24=zeros(kmax,1);
ff25=zeros(kmax,1);
ff26=zeros(kmax,1);

%% trigger flag in ETM-a, if trigger triIS1i=1, else triIS1i=-1
triIS11=zeros(kmax,1);
triIS12=zeros(kmax,1);
triIS13=zeros(kmax,1);
triIS14=zeros(kmax,1);
triIS15=zeros(kmax,1);
triIS16=zeros(kmax,1);

%% trigger flag in ETM-b, if trigger triIS2i=1, else triIS2i=-1
triIS21=zeros(kmax,1);
triIS22=zeros(kmax,1);
triIS23=zeros(kmax,1);
triIS24=zeros(kmax,1);
triIS25=zeros(kmax,1);
triIS26=zeros(kmax,1);

%% trigger number of ETM-a
tricount11=0;
tricount12=0;
tricount13=0;
tricount14=0;
tricount15=0;
tricount16=0;

%% trigger number of ETM-b
tricount21=0;
tricount22=0;
tricount23=0;
tricount24=0;
tricount25=0;
tricount26=0;

for k=1:kmax  
  if k==1
  
     %% The actuator is normal    
      Theta1(k)=1; Theta2(k)=1; Theta3(k)=1; Theta4(k)=1; Theta5(k)=1; Theta6(k)=1;
      phi1(k)=0; phi2(k)=0; phi3(k)=0; phi4(k)=0; phi5(k)=0; phi6(k)=0;

      load data %% load the initial values of the variable x_{i}(t), \eta_{i}(t), and c_{i}(t)
 
       %% Initialization of variables at the initial time
       %% x_{i}(t) in (1) and (2)
        x1(3*k-2:3*k)=aa1;
        x2(3*k-2:3*k)=aa2;
        x3(3*k-2:3*k)=aa3;
        x4(3*k-2:3*k)=aa4;
        x5(3*k-2:3*k)=aa5;
        x6(3*k-2:3*k)=aa6;
        x0(3*k-2:3*k)=aa7;
        
       %% x_{0}(T_{k_{i}}^{i}) in (14)
        x01last(3*k-2:3*k)=x0(3*k-2:3*k);
        x02last(3*k-2:3*k)=x0(3*k-2:3*k);
        x03last(3*k-2:3*k)=x0(3*k-2:3*k);
        x04last(3*k-2:3*k)=x0(3*k-2:3*k);
        x05last(3*k-2:3*k)=x0(3*k-2:3*k);
        x06last(3*k-2:3*k)=x0(3*k-2:3*k);
        
       %% state observer \eta_{i}(t) in (6)
        eita1(3*k-2:3*k)=aa8;
        eita2(3*k-2:3*k)=aa9;
        eita3(3*k-2:3*k)=aa10;
        eita4(3*k-2:3*k)=aa11;
        eita5(3*k-2:3*k)=aa12;
        eita6(3*k-2:3*k)=aa13;
        
       %% \eta_{i}(t_{\sigma_{i}}^{i}) in (9)
        eita1last(3*k-2:3*k)=eita1(3*k-2:3*k);
        eita2last(3*k-2:3*k)=eita2(3*k-2:3*k);
        eita3last(3*k-2:3*k)=eita3(3*k-2:3*k);
        eita4last(3*k-2:3*k)=eita4(3*k-2:3*k);
        eita5last(3*k-2:3*k)=eita5(3*k-2:3*k);
        eita6last(3*k-2:3*k)=eita6(3*k-2:3*k);
        
       %% \eta_{i}(t_{\sigma_{i}}^{i}) in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        eita1last0(3*k-2:3*k)=eita1(3*k-2:3*k);
        eita2last0(3*k-2:3*k)=eita2(3*k-2:3*k);
        eita3last0(3*k-2:3*k)=eita3(3*k-2:3*k);
        eita4last0(3*k-2:3*k)=eita4(3*k-2:3*k);
        eita5last0(3*k-2:3*k)=eita5(3*k-2:3*k);
        eita6last0(3*k-2:3*k)=eita6(3*k-2:3*k);

       %% c_{i}(t) in (16)
        c1(k)=aa14;
        c2(k)=aa15;
        c3(k)=aa16;
        c4(k)=aa17;
        c5(k)=aa18;
        c6(k)=aa19;
      
       %% c_{i}(T_{k_{i}}^{i}) in (14)
        c1last(k)=c1(k);
        c2last(k)=c2(k);
        c3last(k)=c3(k);
        c4last(k)=c4(k);
        c5last(k)=c5(k);
        c6last(k)=c6(k);
        
       %% output at the initial time
        y1(k)=C*x1(3*k-2:3*k);
        y2(k)=C*x2(3*k-2:3*k);
        y3(k)=C*x3(3*k-2:3*k);
        y4(k)=C*x4(3*k-2:3*k);
        y5(k)=C*x5(3*k-2:3*k);
        y6(k)=C*x6(3*k-2:3*k);
        y0(k)=C*x0(3*k-2:3*k);
        
       %% triggering instant t_{\sigma_{i}}^{i} in ETM-a
        triT11=1;
        triT12=1;
        triT13=1;
        triT14=1;
        triT15=1;
        triT16=1;
        
       %% triggering instant t_{\sigma_{i}}^{i}, and \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        triT110=1;
        triT120=1;
        triT130=1;
        triT140=1;
        triT150=1;
        triT160=1;
        
       %% triggering instant T_{k_{i}}^{i} in ETM-b
        triT21=1;
        triT22=1;
        triT23=1;
        triT24=1;
        triT25=1;
        triT26=1;
        
       %% e^{A(t-t_{\sigma_{i}}^{i})} in (10)
        GTstep_S11=expm(A*(k-triT11)*T);
        GTstep_S12=expm(A*(k-triT12)*T);
        GTstep_S13=expm(A*(k-triT13)*T);
        GTstep_S14=expm(A*(k-triT14)*T);
        GTstep_S15=expm(A*(k-triT15)*T);
        GTstep_S16=expm(A*(k-triT16)*T);
        
       %% e^{A(T_{k_{i}}^{i}-t_{\sigma_{i}}^{i})} in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        GTstep_S21=expm(A*(triT21-triT110)*T);
        GTstep_S22=expm(A*(triT22-triT120)*T);
        GTstep_S23=expm(A*(triT23-triT130)*T);
        GTstep_S24=expm(A*(triT24-triT140)*T);
        GTstep_S25=expm(A*(triT25-triT150)*T);
        GTstep_S26=expm(A*(triT26-triT160)*T);

       %% \hat{\zeta}_{i}(T_{k_{ij}}^{ij}) in (14)
        zetaj1last(3*k-2:3*k)=abs(huaA(1,1))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(1,2))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(1,3))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(1,4))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(1,5))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(1,6))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(1,1)*(GTstep_S21*eita1last0(3*k-2:3*k)-S(1,1)*x01last(3*k-2:3*k));      
        zetaj2last(3*k-2:3*k)=abs(huaA(2,1))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(2,2))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(2,3))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(2,4))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(2,5))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(2,6))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(2,2)*(GTstep_S22*eita2last0(3*k-2:3*k)-S(2,2)*x02last(3*k-2:3*k));      
        zetaj3last(3*k-2:3*k)=abs(huaA(3,1))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(3,2))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(3,3))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(3,4))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(3,5))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(3,6))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(3,3)*(GTstep_S23*eita3last0(3*k-2:3*k)-S(3,3)*x03last(3*k-2:3*k));      
        zetaj4last(3*k-2:3*k)=abs(huaA(4,1))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(4,2))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(4,3))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(4,4))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(4,5))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(4,6))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(4,4)*(GTstep_S24*eita4last0(3*k-2:3*k)-S(4,4)*x04last(3*k-2:3*k));      
        zetaj5last(3*k-2:3*k)=abs(huaA(5,1))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(5,2))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(5,3))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(5,4))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(5,5))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(5,6))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(5,5)*(GTstep_S25*eita5last0(3*k-2:3*k)-S(5,5)*x05last(3*k-2:3*k));      
        zetaj6last(3*k-2:3*k)=abs(huaA(6,1))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(6,2))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(6,3))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(6,4))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(6,5))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(6,6))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(6,6)*(GTstep_S26*eita6last0(3*k-2:3*k)-S(6,6)*x06last(3*k-2:3*k));      
         
       %% u_{i}(t) in (14)
        u1(k)=-c1last(k)*K*zetaj1last(3*k-2:3*k)-c1last(k)*hh(K*zetaj1last(3*k-2:3*k));
        u2(k)=-c2last(k)*K*zetaj2last(3*k-2:3*k)-c2last(k)*hh(K*zetaj2last(3*k-2:3*k));
        u3(k)=-c3last(k)*K*zetaj3last(3*k-2:3*k)-c3last(k)*hh(K*zetaj3last(3*k-2:3*k));
        u4(k)=-c4last(k)*K*zetaj4last(3*k-2:3*k)-c4last(k)*hh(K*zetaj4last(3*k-2:3*k));
        u5(k)=-c5last(k)*K*zetaj5last(3*k-2:3*k)-c5last(k)*hh(K*zetaj5last(3*k-2:3*k));
        u6(k)=-c6last(k)*K*zetaj6last(3*k-2:3*k)-c6last(k)*hh(K*zetaj6last(3*k-2:3*k));
        
       %% real controller in (1)
        uf1(k)=Theta1(k)*u1(k)+phi1(k);
        uf2(k)=Theta2(k)*u2(k)+phi2(k);
        uf3(k)=Theta3(k)*u3(k)+phi3(k);
        uf4(k)=Theta4(k)*u4(k)+phi4(k);
        uf5(k)=Theta5(k)*u5(k)+phi5(k);
        uf6(k)=Theta6(k)*u6(k)+phi6(k);

       %% The iteration value of the variables at the next time
       %% x_{i}(t) in (1) and (2)
        x1(3*(k+1)-2:3*(k+1))=x1(3*k-2:3*k)+T*(A*x1(3*k-2:3*k)+B*uf1(k));
        x2(3*(k+1)-2:3*(k+1))=x2(3*k-2:3*k)+T*(A*x2(3*k-2:3*k)+B*uf2(k));
        x3(3*(k+1)-2:3*(k+1))=x3(3*k-2:3*k)+T*(A*x3(3*k-2:3*k)+B*uf3(k));
        x4(3*(k+1)-2:3*(k+1))=x4(3*k-2:3*k)+T*(A*x4(3*k-2:3*k)+B*uf4(k));
        x5(3*(k+1)-2:3*(k+1))=x5(3*k-2:3*k)+T*(A*x5(3*k-2:3*k)+B*uf5(k));
        x6(3*(k+1)-2:3*(k+1))=x6(3*k-2:3*k)+T*(A*x6(3*k-2:3*k)+B*uf6(k));
        x0(3*(k+1)-2:3*(k+1))=x0(3*k-2:3*k)+T*(A*x0(3*k-2:3*k));   
        
       %% state observer \eta_{i}(t) in (6) 
        eita1(3*(k+1)-2:3*(k+1))=eita1(3*k-2:3*k)+T*(A*eita1(3*k-2:3*k)+B*uf1(k)+L*(C*eita1(3*k-2:3*k)-y1(k)));
        eita2(3*(k+1)-2:3*(k+1))=eita2(3*k-2:3*k)+T*(A*eita2(3*k-2:3*k)+B*uf2(k)+L*(C*eita2(3*k-2:3*k)-y2(k)));
        eita3(3*(k+1)-2:3*(k+1))=eita3(3*k-2:3*k)+T*(A*eita3(3*k-2:3*k)+B*uf3(k)+L*(C*eita3(3*k-2:3*k)-y3(k)));
        eita4(3*(k+1)-2:3*(k+1))=eita4(3*k-2:3*k)+T*(A*eita4(3*k-2:3*k)+B*uf4(k)+L*(C*eita4(3*k-2:3*k)-y4(k)));
        eita5(3*(k+1)-2:3*(k+1))=eita5(3*k-2:3*k)+T*(A*eita5(3*k-2:3*k)+B*uf5(k)+L*(C*eita5(3*k-2:3*k)-y5(k)));
        eita6(3*(k+1)-2:3*(k+1))=eita6(3*k-2:3*k)+T*(A*eita6(3*k-2:3*k)+B*uf6(k)+L*(C*eita6(3*k-2:3*k)-y6(k)));
        
       %% c_{i}(t) in (16)
        c1(k+1)=c1(k)+T*k1*(zetaj1last(3*k-2:3*k)'*P*B*B'*P*zetaj1last(3*k-2:3*k)+norm(B'*P*zetaj1last(3*k-2:3*k)));
        c2(k+1)=c2(k)+T*k2*(zetaj2last(3*k-2:3*k)'*P*B*B'*P*zetaj2last(3*k-2:3*k)+norm(B'*P*zetaj2last(3*k-2:3*k)));
        c3(k+1)=c3(k)+T*k3*(zetaj3last(3*k-2:3*k)'*P*B*B'*P*zetaj3last(3*k-2:3*k)+norm(B'*P*zetaj3last(3*k-2:3*k)));
        c4(k+1)=c4(k)+T*k4*(zetaj4last(3*k-2:3*k)'*P*B*B'*P*zetaj4last(3*k-2:3*k)+norm(B'*P*zetaj4last(3*k-2:3*k)));
        c5(k+1)=c5(k)+T*k5*(zetaj5last(3*k-2:3*k)'*P*B*B'*P*zetaj5last(3*k-2:3*k)+norm(B'*P*zetaj5last(3*k-2:3*k)));
        c6(k+1)=c6(k)+T*k6*(zetaj6last(3*k-2:3*k)'*P*B*B'*P*zetaj6last(3*k-2:3*k)+norm(B'*P*zetaj6last(3*k-2:3*k)));
        
       %% e_{i}(t) in (9)
        e1(3*k-2:3*k)=GTstep_S11*eita1last(3*k-2:3*k)-eita1(3*k-2:3*k);
        e2(3*k-2:3*k)=GTstep_S12*eita2last(3*k-2:3*k)-eita2(3*k-2:3*k);
        e3(3*k-2:3*k)=GTstep_S13*eita3last(3*k-2:3*k)-eita3(3*k-2:3*k);
        e4(3*k-2:3*k)=GTstep_S14*eita4last(3*k-2:3*k)-eita4(3*k-2:3*k);
        e5(3*k-2:3*k)=GTstep_S15*eita5last(3*k-2:3*k)-eita5(3*k-2:3*k);  
        e6(3*k-2:3*k)=GTstep_S16*eita6last(3*k-2:3*k)-eita6(3*k-2:3*k);  
        
       %% \epsilon_{i}(t) in (13)
        ee1(3*k-2:3*k)=GTstep_S21*eita1last0(3*k-2:3*k)-GTstep_S11*eita1last(3*k-2:3*k);
        ee2(3*k-2:3*k)=GTstep_S22*eita2last0(3*k-2:3*k)-GTstep_S12*eita2last(3*k-2:3*k);
        ee3(3*k-2:3*k)=GTstep_S23*eita3last0(3*k-2:3*k)-GTstep_S13*eita3last(3*k-2:3*k);
        ee4(3*k-2:3*k)=GTstep_S24*eita4last0(3*k-2:3*k)-GTstep_S14*eita4last(3*k-2:3*k);
        ee5(3*k-2:3*k)=GTstep_S25*eita5last0(3*k-2:3*k)-GTstep_S15*eita5last(3*k-2:3*k);
        ee6(3*k-2:3*k)=GTstep_S26*eita6last0(3*k-2:3*k)-GTstep_S16*eita6last(3*k-2:3*k);

       %% f_{1 i}(t) in (8)
        err11(k)=(1+beta+rrou1*c1last(k)) *(1+1/gama)*rou1*e1(3*k-2:3*k)'*P*B*B'*P*e1(3*k-2:3*k)+...
            2*(2*N1+a10)*(detaphi+detatheta*c1last(k))*norm(e1(3*k-2:3*k)'*P*B);
        thre11(k)=mu1*exp(-v1*(k-1)*T);
        ff11(k)=err11(k)-thre11(k);
        
        err12(k)=(1+beta+rrou1*c2last(k)) *(1+1/gama)*rou2*e2(3*k-2:3*k)'*P*B*B'*P*e2(3*k-2:3*k)+...
            2*(2*N2+a20)*(detaphi+detatheta*c2last(k))*norm(e2(3*k-2:3*k)'*P*B);
        thre12(k)=mu1*exp(-v1*(k-1)*T);
        ff12(k)=err12(k)-thre12(k);
          
        err13(k)=(1+beta+rrou1*c3last(k)) *(1+1/gama)*rou3*e3(3*k-2:3*k)'*P*B*B'*P*e3(3*k-2:3*k)+...
            2*(2*N3+a30)*(detaphi+detatheta*c3last(k))*norm(e3(3*k-2:3*k)'*P*B);
        thre13(k)=mu1*exp(-v1*(k-1)*T);
        ff13(k)=err13(k)-thre13(k);
        
        err14(k)=(1+beta+rrou1*c4last(k)) *(1+1/gama)*rou4*e4(3*k-2:3*k)'*P*B*B'*P*e4(3*k-2:3*k)+...
            2*(2*N4+a40)*(detaphi+detatheta*c4last(k))*norm(e4(3*k-2:3*k)'*P*B);
        thre14(k)=mu1*exp(-v1*(k-1)*T);
        ff14(k)=err14(k)-thre14(k);
        
        err15(k)=(1+beta+rrou1*c5last(k)) *(1+1/gama)*rou5*e5(3*k-2:3*k)'*P*B*B'*P*e5(3*k-2:3*k)+...
            2*(2*N5+a50)*(detaphi+detatheta*c5last(k))*norm(e5(3*k-2:3*k)'*P*B);
        thre15(k)=mu1*exp(-v1*(k-1)*T);
        ff15(k)=err15(k)-thre15(k);
        
        err16(k)=(1+beta+rrou1*c6last(k)) *(1+1/gama)*rou6*e6(3*k-2:3*k)'*P*B*B'*P*e6(3*k-2:3*k)+...
            2*(2*N6+a60)*(detaphi+detatheta*c6last(k))*norm(e6(3*k-2:3*k)'*P*B);
        thre16(k)=mu1*exp(-v1*(k-1)*T);
        ff16(k)=err16(k)-thre16(k);
       
        
      %% f_{2 i}(t) in (12)
       err21(k)= (1+beta+rrou2*c1last(k)) *(1+gama)*ee1(3*k-2:3*k)'*P*B*B'*P*ee1(3*k-2:3*k)+...
           2*(detaphi+detatheta*c1last(k))*norm(ee1(3*k-2:3*k)'*P*B);
       thre21(k)=mu2*exp(-v2*(k-1)*T);
       ff21(k)=err21(k)-thre21(k);
       
       err22(k)= (1+beta+rrou2*c2last(k)) *(1+gama)*ee2(3*k-2:3*k)'*P*B*B'*P*ee2(3*k-2:3*k)+...
           2*(detaphi+detatheta*c2last(k))*norm(ee2(3*k-2:3*k)'*P*B);
       thre22(k)=mu2*exp(-v2*(k-1)*T);
       ff22(k)=err22(k)-thre22(k);
       
       err23(k)= (1+beta+rrou2*c3last(k)) *(1+gama)*ee3(3*k-2:3*k)'*P*B*B'*P*ee3(3*k-2:3*k)+...
           2*(detaphi+detatheta*c3last(k))*norm(ee3(3*k-2:3*k)'*P*B);
       thre23(k)=mu2*exp(-v2*(k-1)*T);
       ff23(k)=err23(k)-thre23(k);
       
       err24(k)= (1+beta+rrou2*c4last(k)) *(1+gama)*ee4(3*k-2:3*k)'*P*B*B'*P*ee4(3*k-2:3*k)+...
           2*(detaphi+detatheta*c4last(k))*norm(ee4(3*k-2:3*k)'*P*B);
       thre24(k)=mu2*exp(-v2*(k-1)*T);
       ff24(k)=err24(k)-thre24(k);
       
       err25(k)= (1+beta+rrou2*c5last(k)) *(1+gama)*ee5(3*k-2:3*k)'*P*B*B'*P*ee5(3*k-2:3*k)+...
           2*(detaphi+detatheta*c5last(k))*norm(ee5(3*k-2:3*k)'*P*B);
       thre25(k)=mu2*exp(-v2*(k-1)*T);
       ff25(k)=err25(k)-thre25(k);
       
       err26(k)= (1+beta+rrou2*c6last(k)) *(1+gama)*ee6(3*k-2:3*k)'*P*B*B'*P*ee6(3*k-2:3*k)+...
           2*(detaphi+detatheta*c6last(k))*norm(ee6(3*k-2:3*k)'*P*B);
       thre26(k)=mu2*exp(-v2*(k-1)*T);
       ff26(k)=err26(k)-thre26(k);
 
       %% trigger flag in ETM-a, if trigger triIS1i=1, else triIS1i=-1
        triIS11(k)=1;
        triIS12(k)=1;
        triIS13(k)=1;
        triIS14(k)=1;
        triIS15(k)=1;
        triIS16(k)=1;
        
       %% trigger flag in ETM-b, if trigger triIS2i=1, else triIS2i=-1
        triIS21(k)=1;
        triIS22(k)=1;
        triIS23(k)=1;
        triIS24(k)=1;
        triIS25(k)=1;
        triIS26(k)=1;
        
       %% trigger number of ETM-a
        tricount11=tricount11+1;
        tricount12=tricount12+1;
        tricount13=tricount13+1;
        tricount14=tricount14+1;
        tricount15=tricount15+1;
        tricount16=tricount16+1;

       %% trigger number of ETM-b
        tricount21=tricount21+1;
        tricount22=tricount22+1;
        tricount23=tricount23+1;
        tricount24=tricount24+1;
        tricount25=tricount25+1;
        tricount26=tricount26+1;
       
  else if 1<k&&k<acmax
          
     %% The actuator is normal         
      Theta1(k)=1; Theta2(k)=1; Theta3(k)=1; Theta4(k)=1; Theta5(k)=1; Theta6(k)=1; 
      phi1(k)=0; phi2(k)=0; phi3(k)=0; phi4(k)=0; phi5(k)=0; phi6(k)=0;

       %% e^{A(t-t_{\sigma_{i}}^{i})} in (10)
        GTstep_S11=expm(A*(k-triT11)*T);
        GTstep_S12=expm(A*(k-triT12)*T);
        GTstep_S13=expm(A*(k-triT13)*T);
        GTstep_S14=expm(A*(k-triT14)*T);
        GTstep_S15=expm(A*(k-triT15)*T);
        GTstep_S16=expm(A*(k-triT16)*T);
        
       %% e^{A(T_{k_{i}}^{i}-t_{\sigma_{i}}^{i})} in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        GTstep_S21=expm(A*(triT21-triT110)*T);
        GTstep_S22=expm(A*(triT22-triT120)*T);
        GTstep_S23=expm(A*(triT23-triT130)*T);
        GTstep_S24=expm(A*(triT24-triT140)*T);
        GTstep_S25=expm(A*(triT25-triT150)*T);
        GTstep_S26=expm(A*(triT26-triT160)*T);
         
       %% judge whether the ETM-a is triggered
        e1(3*k-2:3*k)=GTstep_S11*eita1last(3*(k-1)-2:3*(k-1))-eita1(3*k-2:3*k); 
        err11(k)=(1+beta+rrou1*c1last(k-1)) *(1+1/gama)*rou1*e1(3*k-2:3*k)'*P*B*B'*P*e1(3*k-2:3*k)+...
            2*(2*N1+a10)*(detaphi+detatheta*c1last(k-1))*norm(e1(3*k-2:3*k)'*P*B);
        thre11(k)=mu1*exp(-v1*(k-1)*T); 
        ff11(k)=err11(k)-thre11(k);
        
         if ff11(k)>0
                 eita1last(3*k-2:3*k)=eita1(3*k-2:3*k); %% if tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (9) is updated
                 triT11=k;  %% triggering instant t_{\sigma_{i}}^{i} in ETM-a
                 triIS11(k)=1; %% trigger flag in ETM-a, if trigger triIS1i=1, else triIS1i=-1
                 tricount11=tricount11+1; %% trigger number of ETM-a
                
            else
               eita1last(3*k-2:3*k)=eita1last(3*(k-1)-2:3*(k-1)); %% if not tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (9) remains unchanged
               triIS11(k)=-1; %% trigger flag in ETM-a, if trigger triIS1i=1, else triIS1i=-1
         end 
         
      e2(3*k-2:3*k)=GTstep_S12*eita2last(3*(k-1)-2:3*(k-1))-eita2(3*k-2:3*k);   
      err12(k)=(1+beta+rrou1*c2last(k-1)) *(1+1/gama)*rou2*e2(3*k-2:3*k)'*P*B*B'*P*e2(3*k-2:3*k)+...
            2*(2*N2+a20)*(detaphi+detatheta*c2last(k-1))*norm(e2(3*k-2:3*k)'*P*B);
      thre12(k)=mu1*exp(-v1*(k-1)*T);   
      ff12(k)=err12(k)-thre12(k);
         
     if ff12(k)>0
                 eita2last(3*k-2:3*k)=eita2(3*k-2:3*k);
                 triT12=k;
                 triIS12(k)=1;
                 tricount12=tricount12+1;
                
         else
               eita2last(3*k-2:3*k)=eita2last(3*(k-1)-2:3*(k-1));
               triIS12(k)=-1; 
      end     
         
      e3(3*k-2:3*k)=GTstep_S13*eita3last(3*(k-1)-2:3*(k-1))-eita3(3*k-2:3*k);   
      err13(k)=(1+beta+rrou1*c3last(k-1))*(1+1/gama)*rou3*e3(3*k-2:3*k)'*P*B*B'*P*e3(3*k-2:3*k)+...
            2*(2*N3+a30)*(detaphi+detatheta*c3last(k-1))*norm(e3(3*k-2:3*k)'*P*B);
      thre13(k)=mu1*exp(-v1*(k-1)*T);   
      ff13(k)=err13(k)-thre13(k);
      
       if ff13(k)>0
                 eita3last(3*k-2:3*k)=eita3(3*k-2:3*k);
                 triT13=k;
                 triIS13(k)=1;
                 tricount13=tricount13+1;
                
         else
               eita3last(3*k-2:3*k)=eita3last(3*(k-1)-2:3*(k-1));
               triIS13(k)=-1; 
       end     
            
       e4(3*k-2:3*k)=GTstep_S14*eita4last(3*(k-1)-2:3*(k-1))-eita4(3*k-2:3*k);  
       err14(k)=(1+beta+rrou1*c4last(k-1)) *(1+1/gama)*rou4*e4(3*k-2:3*k)'*P*B*B'*P*e4(3*k-2:3*k)+...
            2*(2*N4+a40)*(detaphi+detatheta*c4last(k-1))*norm(e4(3*k-2:3*k)'*P*B);
       thre14(k)=mu1*exp(-v1*(k-1)*T);
       ff14(k)=err14(k)-thre14(k);
         
        if ff14(k)>0
                 eita4last(3*k-2:3*k)=eita4(3*k-2:3*k);
                 triT14=k;
                 triIS14(k)=1;
                 tricount14=tricount14+1;
                
         else
               eita4last(3*k-2:3*k)=eita4last(3*(k-1)-2:3*(k-1));
               triIS14(k)=-1; 
       end       
         
        e5(3*k-2:3*k)=GTstep_S15*eita5last(3*(k-1)-2:3*(k-1))-eita5(3*k-2:3*k);    
        err15(k)=(1+beta+rrou1*c5last(k-1)) *(1+1/gama)*rou5*e5(3*k-2:3*k)'*P*B*B'*P*e5(3*k-2:3*k)+...
            2*(2*N5+a50)*(detaphi+detatheta*c5last(k-1))*norm(e5(3*k-2:3*k)'*P*B);
        thre15(k)=mu1*exp(-v1*(k-1)*T);
        ff15(k)=err15(k)-thre15(k);
        
         if ff15(k)>0
                 eita5last(3*k-2:3*k)=eita5(3*k-2:3*k);
                 triT15=k;
                 triIS15(k)=1;
                 tricount15=tricount15+1;
                
         else
               eita5last(3*k-2:3*k)=eita5last(3*(k-1)-2:3*(k-1));
               triIS15(k)=-1; 
       end       

       e6(3*k-2:3*k)=GTstep_S16*eita6last(3*(k-1)-2:3*(k-1))-eita6(3*k-2:3*k); 
       err16(k)=(1+beta+rrou1*c6last(k-1)) *(1+1/gama)*rou6*e6(3*k-2:3*k)'*P*B*B'*P*e6(3*k-2:3*k)+...
            2*(2*N6+a60)*(detaphi+detatheta*c6last(k-1))*norm(e6(3*k-2:3*k)'*P*B);
       thre16(k)=mu1*exp(-v1*(k-1)*T);
       ff16(k)=err16(k)-thre16(k);
       
        if ff16(k)>0
                 eita6last(3*k-2:3*k)=eita6(3*k-2:3*k);
                 triT16=k;
                 triIS16(k)=1;
                 tricount16=tricount16+1;
                
         else
               eita6last(3*k-2:3*k)=eita6last(3*(k-1)-2:3*(k-1));
               triIS16(k)=-1; 
        end      
        
      %% judge whether the ETM-b is triggered 
       ee1(3*k-2:3*k)=GTstep_S21*eita1last0(3*(k-1)-2:3*(k-1))-GTstep_S11*eita1last(3*(k-1)-2:3*(k-1));
       err21(k)= (1+beta+rrou2*c1last(k-1)) *(1+gama)*ee1(3*k-2:3*k)'*P*B*B'*P*ee1(3*k-2:3*k)+...
           2*(detaphi+detatheta*c1last(k-1))*norm(ee1(3*k-2:3*k)'*P*B);
       thre21(k)=mu2*exp(-v2*(k-1)*T); 
       ff21(k)=err21(k)-thre21(k);
          
       if ff21(k)>0
                 x01last(3*k-2:3*k)=x0(3*k-2:3*k);  %% if tigger, x_{0}(T_{k_{i}}^{i}) in (14) is updated
                 eita1last0(3*k-2:3*k)=eita1last(3*k-2:3*k); %% if tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (13) is updated, \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
                 c1last(k)=c1(k); %%  if tigger, c_{i}(T_{k_{i}}^{i}) in (14) is updated
                 triT21=k; %% triggering instant t_{\sigma_{i}}^{i} in ETM-b
                 triT110=triT11; %% if tigger, triggering instant t_{\sigma_{i}}^{i} (\sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}) is updated.
                 triIS21(k)=1; %% trigger flag in ETM-b, if trigger triIS2i=1, else triIS2i=-1
                 tricount21=tricount21+1; %% trigger number of ETM-b  
                
       else
               x01last(3*k-2:3*k)=x01last(3*(k-1)-2:3*(k-1)); %% if not tigger, x_{0}(T_{k_{i}}^{i}) in (14) remains unchanged
               eita1last0(3*k-2:3*k)=eita1last0(3*(k-1)-2:3*(k-1)); %% if not tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (13) remains unchanged, \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
               c1last(k)=c1last(k-1); %% if not tigger, c_{i}(T_{k_{i}}^{i}) in (14) remains unchanged.
               triIS21(k)=-1; %% trigger flag in ETM-b, if trigger triIS2i=1, else triIS2i=-1
        end        
        
       ee2(3*k-2:3*k)=GTstep_S22*eita2last0(3*(k-1)-2:3*(k-1))-GTstep_S12*eita2last(3*(k-1)-2:3*(k-1)); 
       err22(k)= (1+beta+rrou2*c2last(k-1)) *(1+gama)*ee2(3*k-2:3*k)'*P*B*B'*P*ee2(3*k-2:3*k)+...
           2*(detaphi+detatheta*c2last(k-1))*norm(ee2(3*k-2:3*k)'*P*B);
       thre22(k)=mu2*exp(-v2*(k-1)*T); 
       ff22(k)=err22(k)-thre22(k);
       
      if ff22(k)>0
                 x02last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita2last0(3*k-2:3*k)=eita2last(3*k-2:3*k);
                 c2last(k)=c2(k);
                 triT22=k;
                 triT120=triT12;
                 triIS22(k)=1;
                 tricount22=tricount22+1;
                
      else
               x02last(3*k-2:3*k)=x02last(3*(k-1)-2:3*(k-1));
               eita2last0(3*k-2:3*k)=eita2last0(3*(k-1)-2:3*(k-1));
               c2last(k)=c2last(k-1);
               triIS22(k)=-1; 
      end       
       
     ee3(3*k-2:3*k)=GTstep_S23*eita3last0(3*(k-1)-2:3*(k-1))-GTstep_S13*eita3last(3*(k-1)-2:3*(k-1));   
     err23(k)= (1+beta+rrou2*c3last(k-1)) *(1+gama)*ee3(3*k-2:3*k)'*P*B*B'*P*ee3(3*k-2:3*k)+...
           2*(detaphi+detatheta*c3last(k-1))*norm(ee3(3*k-2:3*k)'*P*B);
     thre23(k)=mu2*exp(-v2*(k-1)*T); 
     ff23(k)=err23(k)-thre23(k);
      
     if ff23(k)>0
                 x03last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita3last0(3*k-2:3*k)=eita3last(3*k-2:3*k);
                 c3last(k)=c3(k);
                 triT23=k;
                 triT130=triT13;
                 triIS23(k)=1;
                 tricount23=tricount23+1;
                
     else
             x03last(3*k-2:3*k)=x03last(3*(k-1)-2:3*(k-1));
             eita3last0(3*k-2:3*k)=eita3last0(3*(k-1)-2:3*(k-1));
               c3last(k)=c3last(k-1);
               triIS23(k)=-1; 
     end        
       
      ee4(3*k-2:3*k)=GTstep_S24*eita4last0(3*(k-1)-2:3*(k-1))-GTstep_S14*eita4last(3*(k-1)-2:3*(k-1));
      err24(k)= (1+beta+rrou2*c4last(k-1)) *(1+gama)*ee4(3*k-2:3*k)'*P*B*B'*P*ee4(3*k-2:3*k)+...
           2*(detaphi+detatheta*c4last(k-1))*norm(ee4(3*k-2:3*k)'*P*B);
      thre24(k)=mu2*exp(-v2*(k-1)*T);
      ff24(k)=err24(k)-thre24(k);
     
     if ff24(k)>0
                 x04last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita4last0(3*k-2:3*k)=eita4last(3*k-2:3*k);
                 c4last(k)=c4(k);
                 triT24=k;
                 triT140=triT14;
                 triIS24(k)=1;
                 tricount24=tricount24+1;
                
     else
              x04last(3*k-2:3*k)=x04last(3*(k-1)-2:3*(k-1));
              eita4last0(3*k-2:3*k)=eita4last0(3*(k-1)-2:3*(k-1));
               c4last(k)=c4last(k-1);
               triIS24(k)=-1; 
     end     
     
     ee5(3*k-2:3*k)=GTstep_S25*eita5last0(3*(k-1)-2:3*(k-1))-GTstep_S15*eita5last(3*(k-1)-2:3*(k-1));
     err25(k)= (1+beta+rrou2*c5last(k-1)) *(1+gama)*ee5(3*k-2:3*k)'*P*B*B'*P*ee5(3*k-2:3*k)+...
           2*(detaphi+detatheta*c5last(k-1))*norm(ee5(3*k-2:3*k)'*P*B);
     thre25(k)=mu2*exp(-v2*(k-1)*T);
     ff25(k)=err25(k)-thre25(k);
     
     if ff25(k)>0
                 x05last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita5last0(3*k-2:3*k)=eita5last(3*k-2:3*k);
                 c5last(k)=c5(k);
                 triT25=k;
                 triT150=triT15;
                 triIS25(k)=1;
                 tricount25=tricount25+1;
                
     else
              x05last(3*k-2:3*k)=x05last(3*(k-1)-2:3*(k-1));
               eita5last0(3*k-2:3*k)=eita5last0(3*(k-1)-2:3*(k-1));
               c5last(k)=c5last(k-1);
               triIS25(k)=-1; 
     end     
     
     ee6(3*k-2:3*k)=GTstep_S26*eita6last0(3*(k-1)-2:3*(k-1))-GTstep_S16*eita6last(3*(k-1)-2:3*(k-1));
     err26(k)= (1+beta+rrou2*c6last(k-1)) *(1+gama)*ee6(3*k-2:3*k)'*P*B*B'*P*ee6(3*k-2:3*k)+...
           2*(detaphi+detatheta*c6last(k-1))*norm(ee6(3*k-2:3*k)'*P*B);
     thre26(k)=mu2*exp(-v2*(k-1)*T);
     ff26(k)=err26(k)-thre26(k);
     
      if ff26(k)>0
                 x06last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita6last0(3*k-2:3*k)=eita6last(3*k-2:3*k);
                 c6last(k)=c6(k);
                 triT26=k;
                 triT160=triT16;
                 triIS26(k)=1;
                 tricount26=tricount26+1;
                
      else
              x06last(3*k-2:3*k)=x06last(3*(k-1)-2:3*(k-1));
              eita6last0(3*k-2:3*k)=eita6last0(3*(k-1)-2:3*(k-1));
               c6last(k)=c6last(k-1);
               triIS26(k)=-1; 
      end      
           
       %% output
        y1(k)=C*x1(3*k-2:3*k);
        y2(k)=C*x2(3*k-2:3*k);
        y3(k)=C*x3(3*k-2:3*k);
        y4(k)=C*x4(3*k-2:3*k);
        y5(k)=C*x5(3*k-2:3*k);
        y6(k)=C*x6(3*k-2:3*k);
        y0(k)=C*x0(3*k-2:3*k);
       
       %% updata e^{A(t-t_{\sigma_{i}}^{i})} in (10)
        GTstep_S11=expm(A*(k-triT11)*T);
        GTstep_S12=expm(A*(k-triT12)*T);
        GTstep_S13=expm(A*(k-triT13)*T);
        GTstep_S14=expm(A*(k-triT14)*T);
        GTstep_S15=expm(A*(k-triT15)*T);
        GTstep_S16=expm(A*(k-triT16)*T);
        
       %% updata e^{A(T_{k_{i}}^{i}-t_{\sigma_{i}}^{i})} in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        GTstep_S21=expm(A*(triT21-triT110)*T);
        GTstep_S22=expm(A*(triT22-triT120)*T);
        GTstep_S23=expm(A*(triT23-triT130)*T);
        GTstep_S24=expm(A*(triT24-triT140)*T);
        GTstep_S25=expm(A*(triT25-triT150)*T);
        GTstep_S26=expm(A*(triT26-triT160)*T);
        
       %% updata \hat{\zeta}_{i}(T_{k_{ij}}^{ij}) in (14)
        zetaj1last(3*k-2:3*k)=abs(huaA(1,1))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(1,2))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(1,3))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(1,4))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(1,5))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(1,6))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(1,1)*(GTstep_S21*eita1last0(3*k-2:3*k)-S(1,1)*x01last(3*k-2:3*k));      
        zetaj2last(3*k-2:3*k)=abs(huaA(2,1))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(2,2))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(2,3))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(2,4))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(2,5))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(2,6))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(2,2)*(GTstep_S22*eita2last0(3*k-2:3*k)-S(2,2)*x02last(3*k-2:3*k));      
        zetaj3last(3*k-2:3*k)=abs(huaA(3,1))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(3,2))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(3,3))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(3,4))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(3,5))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(3,6))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(3,3)*(GTstep_S23*eita3last0(3*k-2:3*k)-S(3,3)*x03last(3*k-2:3*k));      
        zetaj4last(3*k-2:3*k)=abs(huaA(4,1))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(4,2))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(4,3))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(4,4))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(4,5))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(4,6))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(4,4)*(GTstep_S24*eita4last0(3*k-2:3*k)-S(4,4)*x04last(3*k-2:3*k));      
        zetaj5last(3*k-2:3*k)=abs(huaA(5,1))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(5,2))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(5,3))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(5,4))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(5,5))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(5,6))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(5,5)*(GTstep_S25*eita5last0(3*k-2:3*k)-S(5,5)*x05last(3*k-2:3*k));      
        zetaj6last(3*k-2:3*k)=abs(huaA(6,1))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(6,2))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(6,3))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(6,4))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(6,5))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(6,6))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(6,6)*(GTstep_S26*eita6last0(3*k-2:3*k)-S(6,6)*x06last(3*k-2:3*k));      
         
       %% updata u_{i}(t) in (14)
        u1(k)=-c1last(k)*K*zetaj1last(3*k-2:3*k)-c1last(k)*hh(K*zetaj1last(3*k-2:3*k));
        u2(k)=-c2last(k)*K*zetaj2last(3*k-2:3*k)-c2last(k)*hh(K*zetaj2last(3*k-2:3*k));
        u3(k)=-c3last(k)*K*zetaj3last(3*k-2:3*k)-c3last(k)*hh(K*zetaj3last(3*k-2:3*k));
        u4(k)=-c4last(k)*K*zetaj4last(3*k-2:3*k)-c4last(k)*hh(K*zetaj4last(3*k-2:3*k));
        u5(k)=-c5last(k)*K*zetaj5last(3*k-2:3*k)-c5last(k)*hh(K*zetaj5last(3*k-2:3*k));
        u6(k)=-c6last(k)*K*zetaj6last(3*k-2:3*k)-c6last(k)*hh(K*zetaj6last(3*k-2:3*k));
        
       %% updata real controller in (1)
        uf1(k)=Theta1(k)*u1(k)+phi1(k);
        uf2(k)=Theta2(k)*u2(k)+phi2(k);
        uf3(k)=Theta3(k)*u3(k)+phi3(k);
        uf4(k)=Theta4(k)*u4(k)+phi4(k);
        uf5(k)=Theta5(k)*u5(k)+phi5(k);
        uf6(k)=Theta6(k)*u6(k)+phi6(k);

       %% The iteration value of the variables at the next time
       %% x_{i}(t) in (1) and (2)
        x1(3*(k+1)-2:3*(k+1))=x1(3*k-2:3*k)+T*(A*x1(3*k-2:3*k)+B*uf1(k));
        x2(3*(k+1)-2:3*(k+1))=x2(3*k-2:3*k)+T*(A*x2(3*k-2:3*k)+B*uf2(k));
        x3(3*(k+1)-2:3*(k+1))=x3(3*k-2:3*k)+T*(A*x3(3*k-2:3*k)+B*uf3(k));
        x4(3*(k+1)-2:3*(k+1))=x4(3*k-2:3*k)+T*(A*x4(3*k-2:3*k)+B*uf4(k));
        x5(3*(k+1)-2:3*(k+1))=x5(3*k-2:3*k)+T*(A*x5(3*k-2:3*k)+B*uf5(k));
        x6(3*(k+1)-2:3*(k+1))=x6(3*k-2:3*k)+T*(A*x6(3*k-2:3*k)+B*uf6(k));
        x0(3*(k+1)-2:3*(k+1))=x0(3*k-2:3*k)+T*(A*x0(3*k-2:3*k));   
        
       %% state observer \eta_{i}(t) in (6) 
        eita1(3*(k+1)-2:3*(k+1))=eita1(3*k-2:3*k)+T*(A*eita1(3*k-2:3*k)+B*uf1(k)+L*(C*eita1(3*k-2:3*k)-y1(k)));
        eita2(3*(k+1)-2:3*(k+1))=eita2(3*k-2:3*k)+T*(A*eita2(3*k-2:3*k)+B*uf2(k)+L*(C*eita2(3*k-2:3*k)-y2(k)));
        eita3(3*(k+1)-2:3*(k+1))=eita3(3*k-2:3*k)+T*(A*eita3(3*k-2:3*k)+B*uf3(k)+L*(C*eita3(3*k-2:3*k)-y3(k)));
        eita4(3*(k+1)-2:3*(k+1))=eita4(3*k-2:3*k)+T*(A*eita4(3*k-2:3*k)+B*uf4(k)+L*(C*eita4(3*k-2:3*k)-y4(k)));
        eita5(3*(k+1)-2:3*(k+1))=eita5(3*k-2:3*k)+T*(A*eita5(3*k-2:3*k)+B*uf5(k)+L*(C*eita5(3*k-2:3*k)-y5(k)));
        eita6(3*(k+1)-2:3*(k+1))=eita6(3*k-2:3*k)+T*(A*eita6(3*k-2:3*k)+B*uf6(k)+L*(C*eita6(3*k-2:3*k)-y6(k)));
        
       %% c_{i}(t) in (16)
        c1(k+1)=c1(k)+T*k1*(zetaj1last(3*k-2:3*k)'*P*B*B'*P*zetaj1last(3*k-2:3*k)+norm(B'*P*zetaj1last(3*k-2:3*k)));
        c2(k+1)=c2(k)+T*k2*(zetaj2last(3*k-2:3*k)'*P*B*B'*P*zetaj2last(3*k-2:3*k)+norm(B'*P*zetaj2last(3*k-2:3*k)));
        c3(k+1)=c3(k)+T*k3*(zetaj3last(3*k-2:3*k)'*P*B*B'*P*zetaj3last(3*k-2:3*k)+norm(B'*P*zetaj3last(3*k-2:3*k)));
        c4(k+1)=c4(k)+T*k4*(zetaj4last(3*k-2:3*k)'*P*B*B'*P*zetaj4last(3*k-2:3*k)+norm(B'*P*zetaj4last(3*k-2:3*k)));
        c5(k+1)=c5(k)+T*k6*(zetaj5last(3*k-2:3*k)'*P*B*B'*P*zetaj5last(3*k-2:3*k)+norm(B'*P*zetaj5last(3*k-2:3*k)));
        c6(k+1)=c6(k)+T*k6*(zetaj6last(3*k-2:3*k)'*P*B*B'*P*zetaj6last(3*k-2:3*k)+norm(B'*P*zetaj6last(3*k-2:3*k))); 
  
  else if acmax<=k && k< kmax    

       %% actuator failure
        Theta1(k)=1; Theta2(k)=0.6; Theta3(k)=1; Theta4(k)=1; Theta5(k)=1; Theta6(k)=1;
        phi1(k)=0; phi2(k)=0; phi3(k)=0.2+0.1*sin(0.1*(k-1)*T)+50*exp(-0.1*(k-1)*T); phi4(k)=0; phi5(k)=0; phi6(k)=0;

       %% e^{A(t-t_{\sigma_{i}}^{i})} in (10)
        GTstep_S11=expm(A*(k-triT11)*T);
        GTstep_S12=expm(A*(k-triT12)*T);
        GTstep_S13=expm(A*(k-triT13)*T);
        GTstep_S14=expm(A*(k-triT14)*T);
        GTstep_S15=expm(A*(k-triT15)*T);
        GTstep_S16=expm(A*(k-triT16)*T);
        
       %% e^{A(T_{k_{i}}^{i}-t_{\sigma_{i}}^{i})} in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        GTstep_S21=expm(A*(triT21-triT110)*T);
        GTstep_S22=expm(A*(triT22-triT120)*T);
        GTstep_S23=expm(A*(triT23-triT130)*T);
        GTstep_S24=expm(A*(triT24-triT140)*T);
        GTstep_S25=expm(A*(triT25-triT150)*T);
        GTstep_S26=expm(A*(triT26-triT160)*T);
         
       %% judge whether the ETM-a is triggered
        e1(3*k-2:3*k)=GTstep_S11*eita1last(3*(k-1)-2:3*(k-1))-eita1(3*k-2:3*k); 
        err11(k)=(1+beta+rrou1*c1last(k-1)) *(1+1/gama)*rou1*e1(3*k-2:3*k)'*P*B*B'*P*e1(3*k-2:3*k)+...
            2*(2*N1+a10)*(detaphi+detatheta*c1last(k-1))*norm(e1(3*k-2:3*k)'*P*B);
        thre11(k)=mu1*exp(-v1*(k-1)*T); 
        ff11(k)=err11(k)-thre11(k);
        
         if ff11(k)>0
                 eita1last(3*k-2:3*k)=eita1(3*k-2:3*k); %% if tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (9) is updated
                 triT11=k;  %% triggering instant t_{\sigma_{i}}^{i} in ETM-a
                 triIS11(k)=1; %% trigger flag in ETM-a, if trigger triIS1i=1, else triIS1i=-1
                 tricount11=tricount11+1; %% trigger number of ETM-a
                
            else
               eita1last(3*k-2:3*k)=eita1last(3*(k-1)-2:3*(k-1)); %% if not tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (9) remains unchanged
               triIS11(k)=-1; %% trigger flag in ETM-a, if trigger triIS1i=1, else triIS1i=-1
         end 
         
      e2(3*k-2:3*k)=GTstep_S12*eita2last(3*(k-1)-2:3*(k-1))-eita2(3*k-2:3*k);   
      err12(k)=(1+beta+rrou1*c2last(k-1)) *(1+1/gama)*rou2*e2(3*k-2:3*k)'*P*B*B'*P*e2(3*k-2:3*k)+...
            2*(2*N2+a20)*(detaphi+detatheta*c2last(k-1))*norm(e2(3*k-2:3*k)'*P*B);
      thre12(k)=mu1*exp(-v1*(k-1)*T);   
      ff12(k)=err12(k)-thre12(k);
         
     if ff12(k)>0
                 eita2last(3*k-2:3*k)=eita2(3*k-2:3*k);
                 triT12=k;
                 triIS12(k)=1;
                 tricount12=tricount12+1;
                
         else
               eita2last(3*k-2:3*k)=eita2last(3*(k-1)-2:3*(k-1));
               triIS12(k)=-1; 
      end     
         
      e3(3*k-2:3*k)=GTstep_S13*eita3last(3*(k-1)-2:3*(k-1))-eita3(3*k-2:3*k);   
      err13(k)=(1+beta+rrou1*c3last(k-1))*(1+1/gama)*rou3*e3(3*k-2:3*k)'*P*B*B'*P*e3(3*k-2:3*k)+...
            2*(2*N3+a30)*(detaphi+detatheta*c3last(k-1))*norm(e3(3*k-2:3*k)'*P*B);
      thre13(k)=mu1*exp(-v1*(k-1)*T);   
      ff13(k)=err13(k)-thre13(k);
      
       if ff13(k)>0
                 eita3last(3*k-2:3*k)=eita3(3*k-2:3*k);
                 triT13=k;
                 triIS13(k)=1;
                 tricount13=tricount13+1;
                
         else
               eita3last(3*k-2:3*k)=eita3last(3*(k-1)-2:3*(k-1));
               triIS13(k)=-1; 
       end     
            
       e4(3*k-2:3*k)=GTstep_S14*eita4last(3*(k-1)-2:3*(k-1))-eita4(3*k-2:3*k);  
       err14(k)=(1+beta+rrou1*c4last(k-1)) *(1+1/gama)*rou4*e4(3*k-2:3*k)'*P*B*B'*P*e4(3*k-2:3*k)+...
            2*(2*N4+a40)*(detaphi+detatheta*c4last(k-1))*norm(e4(3*k-2:3*k)'*P*B);
       thre14(k)=mu1*exp(-v1*(k-1)*T);
       ff14(k)=err14(k)-thre14(k);
         
        if ff14(k)>0
                 eita4last(3*k-2:3*k)=eita4(3*k-2:3*k);
                 triT14=k;
                 triIS14(k)=1;
                 tricount14=tricount14+1;
                
         else
               eita4last(3*k-2:3*k)=eita4last(3*(k-1)-2:3*(k-1));
               triIS14(k)=-1; 
       end       
         
        e5(3*k-2:3*k)=GTstep_S15*eita5last(3*(k-1)-2:3*(k-1))-eita5(3*k-2:3*k);    
        err15(k)=(1+beta+rrou1*c5last(k-1)) *(1+1/gama)*rou5*e5(3*k-2:3*k)'*P*B*B'*P*e5(3*k-2:3*k)+...
            2*(2*N5+a50)*(detaphi+detatheta*c5last(k-1))*norm(e5(3*k-2:3*k)'*P*B);
        thre15(k)=mu1*exp(-v1*(k-1)*T);
        ff15(k)=err15(k)-thre15(k);
        
         if ff15(k)>0
                 eita5last(3*k-2:3*k)=eita5(3*k-2:3*k);
                 triT15=k;
                 triIS15(k)=1;
                 tricount15=tricount15+1;
                
         else
               eita5last(3*k-2:3*k)=eita5last(3*(k-1)-2:3*(k-1));
               triIS15(k)=-1; 
       end       

       e6(3*k-2:3*k)=GTstep_S16*eita6last(3*(k-1)-2:3*(k-1))-eita6(3*k-2:3*k); 
       err16(k)=(1+beta+rrou1*c6last(k-1)) *(1+1/gama)*rou6*e6(3*k-2:3*k)'*P*B*B'*P*e6(3*k-2:3*k)+...
            2*(2*N6+a60)*(detaphi+detatheta*c6last(k-1))*norm(e6(3*k-2:3*k)'*P*B);
       thre16(k)=mu1*exp(-v1*(k-1)*T);
       ff16(k)=err16(k)-thre16(k);
       
        if ff16(k)>0
                 eita6last(3*k-2:3*k)=eita6(3*k-2:3*k);
                 triT16=k;
                 triIS16(k)=1;
                 tricount16=tricount16+1;
                
         else
               eita6last(3*k-2:3*k)=eita6last(3*(k-1)-2:3*(k-1));
               triIS16(k)=-1; 
        end      
        
      %% judge whether the ETM-b is triggered 
       ee1(3*k-2:3*k)=GTstep_S21*eita1last0(3*(k-1)-2:3*(k-1))-GTstep_S11*eita1last(3*(k-1)-2:3*(k-1));
       err21(k)= (1+beta+rrou2*c1last(k-1)) *(1+gama)*ee1(3*k-2:3*k)'*P*B*B'*P*ee1(3*k-2:3*k)+...
           2*(detaphi+detatheta*c1last(k-1))*norm(ee1(3*k-2:3*k)'*P*B);
       thre21(k)=mu2*exp(-v2*(k-1)*T); 
       ff21(k)=err21(k)-thre21(k);
          
       if ff21(k)>0
                 x01last(3*k-2:3*k)=x0(3*k-2:3*k);  %% if tigger, x_{0}(T_{k_{i}}^{i}) in (14) is updated
                 eita1last0(3*k-2:3*k)=eita1last(3*k-2:3*k); %% if tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (13) is updated, \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
                 c1last(k)=c1(k); %% if tigger, c_{i}(T_{k_{i}}^{i}) in (14) is updated
                 triT21=k; %% triggering instant t_{\sigma_{i}}^{i} in ETM-b
                 triT110=triT11; %% if tigger, triggering instant t_{\sigma_{i}}^{i} (\sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}) is updated.
                 triIS21(k)=1; %% trigger flag in ETM-b, if trigger triIS2i=1, else triIS2i=-1
                 tricount21=tricount21+1; %% trigger number of ETM-b  
                
       else
               x01last(3*k-2:3*k)=x01last(3*(k-1)-2:3*(k-1)); %% if not tigger, x_{0}(T_{k_{i}}^{i}) in (14) remains unchanged
               eita1last0(3*k-2:3*k)=eita1last0(3*(k-1)-2:3*(k-1)); %% if not tigger, \eta_{i}(t_{\sigma_{i}}^{i}) in (13) remains unchanged, \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
               c1last(k)=c1last(k-1); %% if not tigger, c_{i}(T_{k_{i}}^{i}) in (14) remains unchanged.
               triIS21(k)=-1; %% trigger flag in ETM-b, if trigger triIS2i=1, else triIS2i=-1
        end        
        
       ee2(3*k-2:3*k)=GTstep_S22*eita2last0(3*(k-1)-2:3*(k-1))-GTstep_S12*eita2last(3*(k-1)-2:3*(k-1)); 
       err22(k)= (1+beta+rrou2*c2last(k-1)) *(1+gama)*ee2(3*k-2:3*k)'*P*B*B'*P*ee2(3*k-2:3*k)+...
           2*(detaphi+detatheta*c2last(k-1))*norm(ee2(3*k-2:3*k)'*P*B);
       thre22(k)=mu2*exp(-v2*(k-1)*T); 
       ff22(k)=err22(k)-thre22(k);
       
      if ff22(k)>0
                 x02last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita2last0(3*k-2:3*k)=eita2last(3*k-2:3*k);
                 c2last(k)=c2(k);
                 triT22=k;
                 triT120=triT12;
                 triIS22(k)=1;
                 tricount22=tricount22+1;
                
      else
               x02last(3*k-2:3*k)=x02last(3*(k-1)-2:3*(k-1));
               eita2last0(3*k-2:3*k)=eita2last0(3*(k-1)-2:3*(k-1));
               c2last(k)=c2last(k-1);
               triIS22(k)=-1; 
      end       
       
     ee3(3*k-2:3*k)=GTstep_S23*eita3last0(3*(k-1)-2:3*(k-1))-GTstep_S13*eita3last(3*(k-1)-2:3*(k-1));   
     err23(k)= (1+beta+rrou2*c3last(k-1)) *(1+gama)*ee3(3*k-2:3*k)'*P*B*B'*P*ee3(3*k-2:3*k)+...
           2*(detaphi+detatheta*c3last(k-1))*norm(ee3(3*k-2:3*k)'*P*B);
     thre23(k)=mu2*exp(-v2*(k-1)*T); 
     ff23(k)=err23(k)-thre23(k);
      
     if ff23(k)>0
                 x03last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita3last0(3*k-2:3*k)=eita3last(3*k-2:3*k);
                 c3last(k)=c3(k);
                 triT23=k;
                 triT130=triT13;
                 triIS23(k)=1;
                 tricount23=tricount23+1;
                
     else
             x03last(3*k-2:3*k)=x03last(3*(k-1)-2:3*(k-1));
             eita3last0(3*k-2:3*k)=eita3last0(3*(k-1)-2:3*(k-1));
               c3last(k)=c3last(k-1);
               triIS23(k)=-1; 
     end        
       
      ee4(3*k-2:3*k)=GTstep_S24*eita4last0(3*(k-1)-2:3*(k-1))-GTstep_S14*eita4last(3*(k-1)-2:3*(k-1));
      err24(k)= (1+beta+rrou2*c4last(k-1)) *(1+gama)*ee4(3*k-2:3*k)'*P*B*B'*P*ee4(3*k-2:3*k)+...
           2*(detaphi+detatheta*c4last(k-1))*norm(ee4(3*k-2:3*k)'*P*B);
      thre24(k)=mu2*exp(-v2*(k-1)*T);
      ff24(k)=err24(k)-thre24(k);
     
     if ff24(k)>0
                 x04last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita4last0(3*k-2:3*k)=eita4last(3*k-2:3*k);
                 c4last(k)=c4(k);
                 triT24=k;
                 triT140=triT14;
                 triIS24(k)=1;
                 tricount24=tricount24+1;
                
     else
              x04last(3*k-2:3*k)=x04last(3*(k-1)-2:3*(k-1));
              eita4last0(3*k-2:3*k)=eita4last0(3*(k-1)-2:3*(k-1));
               c4last(k)=c4last(k-1);
               triIS24(k)=-1; 
     end     
     
     ee5(3*k-2:3*k)=GTstep_S25*eita5last0(3*(k-1)-2:3*(k-1))-GTstep_S15*eita5last(3*(k-1)-2:3*(k-1));
     err25(k)= (1+beta+rrou2*c5last(k-1)) *(1+gama)*ee5(3*k-2:3*k)'*P*B*B'*P*ee5(3*k-2:3*k)+...
           2*(detaphi+detatheta*c5last(k-1))*norm(ee5(3*k-2:3*k)'*P*B);
     thre25(k)=mu2*exp(-v2*(k-1)*T);
     ff25(k)=err25(k)-thre25(k);
     
     if ff25(k)>0
                 x05last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita5last0(3*k-2:3*k)=eita5last(3*k-2:3*k);
                 c5last(k)=c5(k);
                 triT25=k;
                 triT150=triT15;
                 triIS25(k)=1;
                 tricount25=tricount25+1;
                
     else
              x05last(3*k-2:3*k)=x05last(3*(k-1)-2:3*(k-1));
               eita5last0(3*k-2:3*k)=eita5last0(3*(k-1)-2:3*(k-1));
               c5last(k)=c5last(k-1);
               triIS25(k)=-1; 
     end     
     
     ee6(3*k-2:3*k)=GTstep_S26*eita6last0(3*(k-1)-2:3*(k-1))-GTstep_S16*eita6last(3*(k-1)-2:3*(k-1));
     err26(k)= (1+beta+rrou2*c6last(k-1)) *(1+gama)*ee6(3*k-2:3*k)'*P*B*B'*P*ee6(3*k-2:3*k)+...
           2*(detaphi+detatheta*c6last(k-1))*norm(ee6(3*k-2:3*k)'*P*B);
     thre26(k)=mu2*exp(-v2*(k-1)*T);
     ff26(k)=err26(k)-thre26(k);
     
      if ff26(k)>0
                 x06last(3*k-2:3*k)=x0(3*k-2:3*k);
                 eita6last0(3*k-2:3*k)=eita6last(3*k-2:3*k);
                 c6last(k)=c6(k);
                 triT26=k;
                 triT160=triT16;
                 triIS26(k)=1;
                 tricount26=tricount26+1;
                
      else
              x06last(3*k-2:3*k)=x06last(3*(k-1)-2:3*(k-1));
              eita6last0(3*k-2:3*k)=eita6last0(3*(k-1)-2:3*(k-1));
               c6last(k)=c6last(k-1);
               triIS26(k)=-1; 
      end      
  
       %% output
        y1(k)=C*x1(3*k-2:3*k);
        y2(k)=C*x2(3*k-2:3*k);
        y3(k)=C*x3(3*k-2:3*k);
        y4(k)=C*x4(3*k-2:3*k);
        y5(k)=C*x5(3*k-2:3*k);
        y6(k)=C*x6(3*k-2:3*k);
        y0(k)=C*x0(3*k-2:3*k);
       
       %% updata e^{A(t-t_{\sigma_{i}}^{i})} in (10)
        GTstep_S11=expm(A*(k-triT11)*T);
        GTstep_S12=expm(A*(k-triT12)*T);
        GTstep_S13=expm(A*(k-triT13)*T);
        GTstep_S14=expm(A*(k-triT14)*T);
        GTstep_S15=expm(A*(k-triT15)*T);
        GTstep_S16=expm(A*(k-triT16)*T);
        
       %% updata e^{A(T_{k_{i}}^{i}-t_{\sigma_{i}}^{i})} in (13), \sigma_{i} \triangleq \arg \min _{l \in \mathbb{Z}: T_{k_{i}}^{i} \geq t_{l}^{i}}\{T_{k_{i}}^{i}-t_{l}^{i}\}.
        GTstep_S21=expm(A*(triT21-triT110)*T);
        GTstep_S22=expm(A*(triT22-triT120)*T);
        GTstep_S23=expm(A*(triT23-triT130)*T);
        GTstep_S24=expm(A*(triT24-triT140)*T);
        GTstep_S25=expm(A*(triT25-triT150)*T);
        GTstep_S26=expm(A*(triT26-triT160)*T);
        
       %% updata \hat{\zeta}_{i}(T_{k_{ij}}^{ij}) in (14)
        zetaj1last(3*k-2:3*k)=abs(huaA(1,1))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(1,2))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(1,3))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(1,4))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(1,5))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(1,6))*(GTstep_S21*eita1last0(3*k-2:3*k)-sign(huaA(1,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(1,1)*(GTstep_S21*eita1last0(3*k-2:3*k)-S(1,1)*x01last(3*k-2:3*k));      
        zetaj2last(3*k-2:3*k)=abs(huaA(2,1))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(2,2))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(2,3))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(2,4))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(2,5))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(2,6))*(GTstep_S22*eita2last0(3*k-2:3*k)-sign(huaA(2,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(2,2)*(GTstep_S22*eita2last0(3*k-2:3*k)-S(2,2)*x02last(3*k-2:3*k));      
        zetaj3last(3*k-2:3*k)=abs(huaA(3,1))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(3,2))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(3,3))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(3,4))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(3,5))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(3,6))*(GTstep_S23*eita3last0(3*k-2:3*k)-sign(huaA(3,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(3,3)*(GTstep_S23*eita3last0(3*k-2:3*k)-S(3,3)*x03last(3*k-2:3*k));      
        zetaj4last(3*k-2:3*k)=abs(huaA(4,1))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(4,2))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(4,3))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(4,4))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(4,5))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(4,6))*(GTstep_S24*eita4last0(3*k-2:3*k)-sign(huaA(4,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(4,4)*(GTstep_S24*eita4last0(3*k-2:3*k)-S(4,4)*x04last(3*k-2:3*k));      
        zetaj5last(3*k-2:3*k)=abs(huaA(5,1))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(5,2))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(5,3))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(5,4))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(5,5))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(5,6))*(GTstep_S25*eita5last0(3*k-2:3*k)-sign(huaA(5,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(5,5)*(GTstep_S25*eita5last0(3*k-2:3*k)-S(5,5)*x05last(3*k-2:3*k));      
        zetaj6last(3*k-2:3*k)=abs(huaA(6,1))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,1))*GTstep_S21*eita1last0(3*k-2:3*k))+abs(huaA(6,2))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,2))*GTstep_S22*eita2last0(3*k-2:3*k))+abs(huaA(6,3))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,3))*GTstep_S23*eita3last0(3*k-2:3*k))+abs(huaA(6,4))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,4))*GTstep_S24*eita4last0(3*k-2:3*k))+abs(huaA(6,5))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,5))*GTstep_S25*eita5last0(3*k-2:3*k))+abs(huaA(6,6))*(GTstep_S26*eita6last0(3*k-2:3*k)-sign(huaA(6,6))*GTstep_S26*eita6last0(3*k-2:3*k))+huaB(6,6)*(GTstep_S26*eita6last0(3*k-2:3*k)-S(6,6)*x06last(3*k-2:3*k));      
         
       %% updata u_{i}(t) in (14)
        u1(k)=-c1last(k)*K*zetaj1last(3*k-2:3*k)-c1last(k)*hh(K*zetaj1last(3*k-2:3*k));
        u2(k)=-c2last(k)*K*zetaj2last(3*k-2:3*k)-c2last(k)*hh(K*zetaj2last(3*k-2:3*k));
        u3(k)=-c3last(k)*K*zetaj3last(3*k-2:3*k)-c3last(k)*hh(K*zetaj3last(3*k-2:3*k));
        u4(k)=-c4last(k)*K*zetaj4last(3*k-2:3*k)-c4last(k)*hh(K*zetaj4last(3*k-2:3*k));
        u5(k)=-c5last(k)*K*zetaj5last(3*k-2:3*k)-c5last(k)*hh(K*zetaj5last(3*k-2:3*k));
        u6(k)=-c6last(k)*K*zetaj6last(3*k-2:3*k)-c6last(k)*hh(K*zetaj6last(3*k-2:3*k));
        
       %% updata real controller in (1)
        uf1(k)=Theta1(k)*u1(k)+phi1(k);
        uf2(k)=Theta2(k)*u2(k)+phi2(k);
        uf3(k)=Theta3(k)*u3(k)+phi3(k);
        uf4(k)=Theta4(k)*u4(k)+phi4(k);
        uf5(k)=Theta5(k)*u5(k)+phi5(k);
        uf6(k)=Theta6(k)*u6(k)+phi6(k);

       %% The iteration value of the variables at the next time
       %% x_{i}(t) in (1) and (2)
        x1(3*(k+1)-2:3*(k+1))=x1(3*k-2:3*k)+T*(A*x1(3*k-2:3*k)+B*uf1(k));
        x2(3*(k+1)-2:3*(k+1))=x2(3*k-2:3*k)+T*(A*x2(3*k-2:3*k)+B*uf2(k));
        x3(3*(k+1)-2:3*(k+1))=x3(3*k-2:3*k)+T*(A*x3(3*k-2:3*k)+B*uf3(k));
        x4(3*(k+1)-2:3*(k+1))=x4(3*k-2:3*k)+T*(A*x4(3*k-2:3*k)+B*uf4(k));
        x5(3*(k+1)-2:3*(k+1))=x5(3*k-2:3*k)+T*(A*x5(3*k-2:3*k)+B*uf5(k));
        x6(3*(k+1)-2:3*(k+1))=x6(3*k-2:3*k)+T*(A*x6(3*k-2:3*k)+B*uf6(k));
        x0(3*(k+1)-2:3*(k+1))=x0(3*k-2:3*k)+T*(A*x0(3*k-2:3*k));   
        
       %% state observer \eta_{i}(t) in (6) 
        eita1(3*(k+1)-2:3*(k+1))=eita1(3*k-2:3*k)+T*(A*eita1(3*k-2:3*k)+B*uf1(k)+L*(C*eita1(3*k-2:3*k)-y1(k)));
        eita2(3*(k+1)-2:3*(k+1))=eita2(3*k-2:3*k)+T*(A*eita2(3*k-2:3*k)+B*uf2(k)+L*(C*eita2(3*k-2:3*k)-y2(k)));
        eita3(3*(k+1)-2:3*(k+1))=eita3(3*k-2:3*k)+T*(A*eita3(3*k-2:3*k)+B*uf3(k)+L*(C*eita3(3*k-2:3*k)-y3(k)));
        eita4(3*(k+1)-2:3*(k+1))=eita4(3*k-2:3*k)+T*(A*eita4(3*k-2:3*k)+B*uf4(k)+L*(C*eita4(3*k-2:3*k)-y4(k)));
        eita5(3*(k+1)-2:3*(k+1))=eita5(3*k-2:3*k)+T*(A*eita5(3*k-2:3*k)+B*uf5(k)+L*(C*eita5(3*k-2:3*k)-y5(k)));
        eita6(3*(k+1)-2:3*(k+1))=eita6(3*k-2:3*k)+T*(A*eita6(3*k-2:3*k)+B*uf6(k)+L*(C*eita6(3*k-2:3*k)-y6(k)));
        
       %% c_{i}(t) in (16)
        c1(k+1)=c1(k)+T*k1*(zetaj1last(3*k-2:3*k)'*P*B*B'*P*zetaj1last(3*k-2:3*k)+norm(B'*P*zetaj1last(3*k-2:3*k)));
        c2(k+1)=c2(k)+T*k2*(zetaj2last(3*k-2:3*k)'*P*B*B'*P*zetaj2last(3*k-2:3*k)+norm(B'*P*zetaj2last(3*k-2:3*k)));
        c3(k+1)=c3(k)+T*k3*(zetaj3last(3*k-2:3*k)'*P*B*B'*P*zetaj3last(3*k-2:3*k)+norm(B'*P*zetaj3last(3*k-2:3*k)));
        c4(k+1)=c4(k)+T*k4*(zetaj4last(3*k-2:3*k)'*P*B*B'*P*zetaj4last(3*k-2:3*k)+norm(B'*P*zetaj4last(3*k-2:3*k)));
        c5(k+1)=c5(k)+T*k6*(zetaj5last(3*k-2:3*k)'*P*B*B'*P*zetaj5last(3*k-2:3*k)+norm(B'*P*zetaj5last(3*k-2:3*k)));
        c6(k+1)=c6(k)+T*k6*(zetaj6last(3*k-2:3*k)'*P*B*B'*P*zetaj6last(3*k-2:3*k)+norm(B'*P*zetaj6last(3*k-2:3*k))); 
       
      else
  
       %%  output
        y1(k)=C*x1(3*k-2:3*k);
        y2(k)=C*x2(3*k-2:3*k);
        y3(k)=C*x3(3*k-2:3*k);
        y4(k)=C*x4(3*k-2:3*k);
        y5(k)=C*x5(3*k-2:3*k);
        y6(k)=C*x6(3*k-2:3*k);
        y0(k)=C*x0(3*k-2:3*k);
    
    end  
   end  
  end
end

%% The observer errors $\delta_{i}(t)$ of six agents
figure(1)
subplot(3,1,1)
plot(0:T:T*(kmax-1),eita1(1:3:length(eita1))-x1(1:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita2(1:3:length(eita1))-x2(1:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita3(1:3:length(eita1))-x3(1:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita4(1:3:length(eita1))-x4(1:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita5(1:3:length(eita1))-x5(1:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita6(1:3:length(eita1))-x6(1:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
axis([0 30 -Inf Inf ]);
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$\delta_{i1}(t)$','Interpreter','latex','FontSize',24); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',22);
set(h,'Orientation','horizon')
subplot(3,1,2)
plot(0:T:T*(kmax-1),eita1(2:3:length(eita1))-x1(2:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita2(2:3:length(eita1))-x2(2:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita3(2:3:length(eita1))-x3(2:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita4(2:3:length(eita1))-x4(2:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita5(2:3:length(eita1))-x5(2:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita6(2:3:length(eita1))-x6(2:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
axis([0 30 -Inf Inf ]);
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$\delta_{i2}(t)$','Interpreter','latex','FontSize',24); 
subplot(3,1,3)
plot(0:T:T*(kmax-1),eita1(3:3:length(eita1))-x1(3:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita2(3:3:length(eita1))-x2(3:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita3(3:3:length(eita1))-x3(3:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita4(3:3:length(eita1))-x4(3:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita5(3:3:length(eita1))-x5(3:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita6(3:3:length(eita1))-x6(3:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
axis([0 30 -Inf Inf ]);
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$\delta_{i3}(t)$','Interpreter','latex','FontSize',24);


%% The adaptive coupling weights $c_{i}(t)$ of six agents
figure(2)
plot(0:T:T*(kmax-1),c1,'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),c2,'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),c3,'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),c4,'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),c5,'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),c6,'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$c_{i}(t)$','Interpreter','latex','FontSize',24); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',22);
set(h,'Orientation','horizon')

%% The triggering instants of ETM-a for six agents
figure(3); 
scatter(0:T:T*(kmax-1),1*triIS11,15,'ko');
hold on;
scatter(0:T:T*(kmax-1),2*triIS12,15,'r*');
hold on;
scatter(0:T:T*(kmax-1),3*triIS13,15,'bd');
hold on;
scatter(0:T:T*(kmax-1),4*triIS14,15,'k>');
hold on;
scatter(0:T:T*(kmax-1),5*triIS15,15,'rs');
hold on;
scatter(0:T:T*(kmax-1),6*triIS16,15,'bh');
hold on;
axis([0,5,0,7]);
set(gca,'FontSize',22,'FontName','Arial');
box on;
set(gca,'yTick',[1:1:6]);
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel({'trigger instants of six agents'},'Interpreter','latex','FontSize',24);
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',22);
set(h,'Orientation','horizon')


%% The triggering instants of ETM-b for six agents
figure(4);
scatter(0:T:T*(kmax-1),1*triIS21,15,'ko');
hold on;
scatter(0:T:T*(kmax-1),2*triIS22,15,'r*');
hold on;
scatter(0:T:T*(kmax-1),3*triIS23,15,'bd');
hold on;
scatter(0:T:T*(kmax-1),4*triIS24,15,'k>');
hold on;
scatter(0:T:T*(kmax-1),5*triIS25,15,'rs');
hold on;
scatter(0:T:T*(kmax-1),6*triIS26,15,'bh');
hold on;
axis([0,5,0,7]);
set(gca,'FontSize',22,'FontName','Arial');
box on;
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
set(gca,'yTick',[1:1:6]);
ylabel({'trigger instants of six agents'},'Interpreter','latex','FontSize',24);
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',22);
set(h,'Orientation','horizon')

%% The state trajectories x_{i1}(t) of all agents
figure(5)
plot(0:T:T*(kmax-1),x1(1:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x2(1:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x3(1:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x4(1:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x5(1:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x6(1:3:length(eita1)),'-.y','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x0(1:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$x_{i1}(t)$','Interpreter','latex','FontSize',24); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6','Leader');
set(h,'Fontsize',22);

%% The state trajectories x_{i2}(t) of all agents
figure(6)
plot(0:T:T*(kmax-1),x1(2:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x2(2:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x3(2:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x4(2:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x5(2:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x6(2:3:length(eita1)),'-.y','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x0(2:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$x_{i2}(t)$','Interpreter','latex','FontSize',24); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6','Leader');
set(h,'Fontsize',22);

%% The state trajectories x_{i3}(t) of all agents
figure(7)
plot(0:T:T*(kmax-1),x1(3:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x2(3:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x3(3:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x4(3:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x5(3:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x6(3:3:length(eita1)),'-.y','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),x0(3:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$x_{i3}(t)$','Interpreter','latex','FontSize',24); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6','Leader');
set(h,'Fontsize',22);

%% The consensus errors \xi_{i}(t) of six agents
figure(8)
subplot(3,1,1)
plot(0:T:T*(kmax-1),eita1(1:3:length(eita1))-s1*x0(1:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita2(1:3:length(eita1))-s2*x0(1:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita3(1:3:length(eita1))-s3*x0(1:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita4(1:3:length(eita1))-s4*x0(1:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita5(1:3:length(eita1))-s5*x0(1:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita6(1:3:length(eita1))-s6*x0(1:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$\xi_{i1}(t)$','Interpreter','latex','FontSize',24); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',22);
set(h,'Orientation','horizon')
subplot(3,1,2)
plot(0:T:T*(kmax-1),eita1(2:3:length(eita1))-s1*x0(2:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita2(2:3:length(eita1))-s2*x0(2:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita3(2:3:length(eita1))-s3*x0(2:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita4(2:3:length(eita1))-s4*x0(2:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita5(2:3:length(eita1))-s5*x0(2:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita6(2:3:length(eita1))-s6*x0(2:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$\xi_{i2}(t)$','Interpreter','latex','FontSize',24); 
subplot(3,1,3)
plot(0:T:T*(kmax-1),eita1(3:3:length(eita1))-s1*x0(3:3:length(eita1)),'-.k','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita2(3:3:length(eita1))-s2*x0(3:3:length(eita1)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita3(3:3:length(eita1))-s3*x0(3:3:length(eita1)),'-.b','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita4(3:3:length(eita1))-s4*x0(3:3:length(eita1)),'-.g','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita5(3:3:length(eita1))-s5*x0(3:3:length(eita1)),'-.m','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold on;
plot(0:T:T*(kmax-1),eita6(3:3:length(eita1))-s6*x0(3:3:length(eita1)),'-.c','LineWidth',3,'MarkerFaceColor','r','MarkerSize',2)
hold off;
set(gca,'FontSize',22,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',24);
ylabel('$\xi_{i3}(t)$','Interpreter','latex','FontSize',24); 


