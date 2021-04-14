close all;
clear;
clc;
%% The Simulation of Numerical Example in Section V %%

%% System dynamics
A1=[0 1;0 -1];B1=[0 1]';C1=[1 0];
A2=[0 1;0 -1];B2=[0 1]';C2=[1 0];
A3=[0 1 0;0 0 1;0 -1 -2];B3=[0;0;1];C3=[1 0 0];
A4=[0 1 0;0 0 1;0 -1 -2];B4=[0;0;1];C4=[1 0 0];
A5=[0 1 0 0;0 1 0 1;1 0 0 1;1 0 -1 -2];B5=[0;0;0;1];C5=[1 0 0 0];
A6=[0 1 0 0;0 1 0 1;1 0 0 1;1 0 -1 -2];B6=[0;0;0;1];C6=[1 0 0 0];
A0=[0 1;-1 0];
C0=[1 0];

%% The solutions (\Phi_{i},\Psi_{i}) of Equation (28) 
xPhi1=[1 0;0 1];  %% \Phi_{1}=xPhi1
xPhi2=[1 0;0 1];  
xPhi3=[1 0;0 1;-1 0];  
xPhi4=[1 0;0 1;-1 0];
xPhi5=[1 0;0 1;-1 0;-1 -1];
xPhi6=[1 0;0 1;-1 0;-1 -1];

dPhi1=[-1 1];  %% \Psi_{1}=dPhi1
dPhi2=[-1 1];
dPhi3=[-2 0];
dPhi4=[-2 0];
dPhi5=[-3 -3];
dPhi6=[-3 -3];

%% Select K1i such that Assumption 3 holds.
K11=[-27 -27]; K21=dPhi1-K11*xPhi1; 
K12=[-27 -27]; K22=dPhi2-K12*xPhi2;
K13=[-27 -26 -7]; K23=dPhi3-K13*xPhi3;
K14=[-27 -26 -7]; K24=dPhi4-K14*xPhi4;
K15=[-27 -26 -6 -10]; K25=dPhi5-K15*xPhi5;
K16=[-27 -26 -6 -10]; K26=dPhi6-K16*xPhi6;

%% Select Fi such that Assumption 4 holds.
F1=[-8;-19]; 
F2=[-8;-19];
F3=[-7;-12;4];
F4=[-7;-12;4];
F5=[-6;-10;0;-10];
F6=[-6;-10;0;-10];

%% Communication topology
huaA=[0 -1 -1 -1 -1 1;
       -1 0 1 0 0 0;
       -1 1 0 0 0 0;
       -1 0 0 0 1 0;
       -1 0 0 1 0 -1;
       1 0 0 0 -1 0];  
huaD=diag([4 2 2 2 3 2]);
huaL=huaD-huaA; 
huaB=diag([1 0 0 0 0 1]);
huaH=huaB+huaL;
D=diag([1 -1 -1 -1 -1 1]);

%% The solution P of Riccati inequality (25)
n=length(A0);
setlmis([]);          
X=lmivar(1,[n 1]);  %% X=P^{-1}
lmiterm([-1 1 1 X],1,1);   %% LMI #1: X>0
lmiterm([2 1 1 X],A0,1,'s');  %% LMI #2: A0*X+X*A0'
lmiterm([2 1 1 0],-eye(n));   %% LMI #2: -eye(n)
lmiterm([2 2 1 X],1,1);   %% LMI #2: X
lmiterm([2 2 2 0],-eye(n));  %% LMI #2: -eye(2) 
lmis=getlmis;
[tminn,xfeas]=feasp(lmis);
tminn
X=dec2mat(lmis,xfeas,X);
P=inv(X);
K=-P; %% the parameter K in compensator (19)
Tao=P^2; %% the parameter \Gamma=Tao in compensator (19)

%% The parameters \delta, \mu_{ij} and \nu_{ij} of two triggering functions (23) and (24) are given 
deta=2; %% \delta=deta
mu=0.3; %% \mu_{ij}=\mu_{i0}=mu
v=0.2; %% \nu_{ij}=\nu_{i0}=v

%% The parameters \iota_{ij} and \iota_{i0} of compensator (19) are given
l11=1;l12=1;l13=1;l14=1;l15=1;l16=1; %% \iota_{ij}
l21=1;l22=1;l23=1;l24=1;l25=1;l26=1;
l31=1;l32=1;l33=1;l34=1;l35=1;l36=1;
l41=1;l42=1;l43=1;l44=1;l45=1;l46=1;
l51=1;l52=1;l53=1;l54=1;l55=1;l56=1;
l61=1;l62=1;l63=1;l64=1;l65=1;l66=1;

l1=1; l2=1; l3=1; l4=1; l5=1; l6=1; %% \iota_{i0}

tmax=25; %% maximum simulation time (s)
T=0.01; %% Iteration step (s)
kmax=1+tmax/T; %% maximum number of iterations

[Num Num]=size(huaL); 
[Ord_A1 Ord_A1]=size(A1);
[Ord_A2 Ord_A2]=size(A2);
[Ord_A3 Ord_A3]=size(A3);
[Ord_A4 Ord_A4]=size(A4);
[Ord_A5 Ord_A5]=size(A5);
[Ord_A6 Ord_A6]=size(A6);
[Ord_S Ord_S]=size(A0);

%% All state variables are set to 0
%% x_{i}(t) in (1) and (2)
x1=zeros(Ord_A1*kmax,1);
x2=zeros(Ord_A2*kmax,1);
x3=zeros(Ord_A3*kmax,1);
x4=zeros(Ord_A4*kmax,1);
x5=zeros(Ord_A5*kmax,1);
x6=zeros(Ord_A6*kmax,1);
x0=zeros(Ord_S*kmax,1);

%% y_{i}(t) in (1) and (2)
y1=zeros(kmax,1);
y2=zeros(kmax,1);
y3=zeros(kmax,1);
y4=zeros(kmax,1);
y5=zeros(kmax,1);
y6=zeros(kmax,1);
y0=zeros(kmax,1);

%% u_{i}(t) in (1) 
u1=zeros(kmax,1);
u2=zeros(kmax,1);
u3=zeros(kmax,1);
u4=zeros(kmax,1);
u5=zeros(kmax,1);
u6=zeros(kmax,1);

%% Compensator \vartheta_{i}(t) in (19a)
kesi1=zeros(Ord_S*kmax,1);
kesi2=zeros(Ord_S*kmax,1);
kesi3=zeros(Ord_S*kmax,1);
kesi4=zeros(Ord_S*kmax,1);
kesi5=zeros(Ord_S*kmax,1);
kesi6=zeros(Ord_S*kmax,1);

%% Observer \psi_{i}(t) in (27b)
eita1=zeros(Ord_A1*kmax,1);
eita2=zeros(Ord_A2*kmax,1);
eita3=zeros(Ord_A3*kmax,1);
eita4=zeros(Ord_A4*kmax,1);
eita5=zeros(Ord_A5*kmax,1);
eita6=zeros(Ord_A6*kmax,1);

%% \vartheta_{i}(t_{\sigma}^{ij})
kesi11last=zeros(Ord_S*kmax,1);
kesi12last=zeros(Ord_S*kmax,1);
kesi13last=zeros(Ord_S*kmax,1);
kesi14last=zeros(Ord_S*kmax,1);
kesi15last=zeros(Ord_S*kmax,1);
kesi16last=zeros(Ord_S*kmax,1);

kesi21last=zeros(Ord_S*kmax,1);
kesi22last=zeros(Ord_S*kmax,1);
kesi23last=zeros(Ord_S*kmax,1);
kesi24last=zeros(Ord_S*kmax,1);
kesi25last=zeros(Ord_S*kmax,1);
kesi26last=zeros(Ord_S*kmax,1);

kesi31last=zeros(Ord_S*kmax,1);
kesi32last=zeros(Ord_S*kmax,1);
kesi33last=zeros(Ord_S*kmax,1);
kesi34last=zeros(Ord_S*kmax,1);
kesi35last=zeros(Ord_S*kmax,1);
kesi36last=zeros(Ord_S*kmax,1);

kesi41last=zeros(Ord_S*kmax,1);
kesi42last=zeros(Ord_S*kmax,1);
kesi43last=zeros(Ord_S*kmax,1);
kesi44last=zeros(Ord_S*kmax,1);
kesi45last=zeros(Ord_S*kmax,1);
kesi46last=zeros(Ord_S*kmax,1);

kesi51last=zeros(Ord_S*kmax,1);
kesi52last=zeros(Ord_S*kmax,1);
kesi53last=zeros(Ord_S*kmax,1);
kesi54last=zeros(Ord_S*kmax,1);
kesi55last=zeros(Ord_S*kmax,1);
kesi56last=zeros(Ord_S*kmax,1);

kesi61last=zeros(Ord_S*kmax,1);
kesi62last=zeros(Ord_S*kmax,1);
kesi63last=zeros(Ord_S*kmax,1);
kesi64last=zeros(Ord_S*kmax,1);
kesi65last=zeros(Ord_S*kmax,1);
kesi66last=zeros(Ord_S*kmax,1);

%% \vartheta_{i}(t_{\sigma}^{i0})
kesi10last=zeros(Ord_S*kmax,1);
kesi20last=zeros(Ord_S*kmax,1);
kesi30last=zeros(Ord_S*kmax,1);
kesi40last=zeros(Ord_S*kmax,1);
kesi50last=zeros(Ord_S*kmax,1);
kesi60last=zeros(Ord_S*kmax,1);

%% f_{ij}(t) in (23)
f11=zeros(kmax,1);
f12=zeros(kmax,1);
f13=zeros(kmax,1);
f14=zeros(kmax,1);
f15=zeros(kmax,1);
f16=zeros(kmax,1);

f21=zeros(kmax,1);
f22=zeros(kmax,1);
f23=zeros(kmax,1);
f24=zeros(kmax,1);
f25=zeros(kmax,1);
f26=zeros(kmax,1);

f31=zeros(kmax,1);
f32=zeros(kmax,1);
f33=zeros(kmax,1);
f34=zeros(kmax,1);
f35=zeros(kmax,1);
f36=zeros(kmax,1);

f41=zeros(kmax,1);
f42=zeros(kmax,1);
f43=zeros(kmax,1);
f44=zeros(kmax,1);
f45=zeros(kmax,1);
f46=zeros(kmax,1);

f51=zeros(kmax,1);
f52=zeros(kmax,1);
f53=zeros(kmax,1);
f54=zeros(kmax,1);
f55=zeros(kmax,1);
f56=zeros(kmax,1);

f61=zeros(kmax,1);
f62=zeros(kmax,1);
f63=zeros(kmax,1);
f64=zeros(kmax,1);
f65=zeros(kmax,1);
f66=zeros(kmax,1);

%% f_{i0}(t) in (24)
f10=zeros(kmax,1);
f20=zeros(kmax,1);
f30=zeros(kmax,1);
f40=zeros(kmax,1);
f50=zeros(kmax,1);
f60=zeros(kmax,1);

%% e_{ij}(t) in (5a)
e11=zeros(Ord_S*kmax,1);
e12=zeros(Ord_S*kmax,1);
e13=zeros(Ord_S*kmax,1);
e14=zeros(Ord_S*kmax,1);
e15=zeros(Ord_S*kmax,1);
e16=zeros(Ord_S*kmax,1);

e21=zeros(Ord_S*kmax,1);
e22=zeros(Ord_S*kmax,1);
e23=zeros(Ord_S*kmax,1);
e24=zeros(Ord_S*kmax,1);
e25=zeros(Ord_S*kmax,1);
e26=zeros(Ord_S*kmax,1);

e31=zeros(Ord_S*kmax,1);
e32=zeros(Ord_S*kmax,1);
e33=zeros(Ord_S*kmax,1);
e34=zeros(Ord_S*kmax,1);
e35=zeros(Ord_S*kmax,1);
e36=zeros(Ord_S*kmax,1);

e41=zeros(Ord_S*kmax,1);
e42=zeros(Ord_S*kmax,1);
e43=zeros(Ord_S*kmax,1);
e44=zeros(Ord_S*kmax,1);
e45=zeros(Ord_S*kmax,1);
e46=zeros(Ord_S*kmax,1);

e51=zeros(Ord_S*kmax,1);
e52=zeros(Ord_S*kmax,1);
e53=zeros(Ord_S*kmax,1);
e54=zeros(Ord_S*kmax,1);
e55=zeros(Ord_S*kmax,1);
e56=zeros(Ord_S*kmax,1);

e61=zeros(Ord_S*kmax,1);
e62=zeros(Ord_S*kmax,1);
e63=zeros(Ord_S*kmax,1);
e64=zeros(Ord_S*kmax,1);
e65=zeros(Ord_S*kmax,1);
e66=zeros(Ord_S*kmax,1);

%% e_{i0}(t) in (5b)
e10=zeros(Ord_S*kmax,1);
e20=zeros(Ord_S*kmax,1);
e30=zeros(Ord_S*kmax,1);
e40=zeros(Ord_S*kmax,1);
e50=zeros(Ord_S*kmax,1);
e60=zeros(Ord_S*kmax,1);

%% c_{ij}(t) in (19b)
c11=zeros(kmax,1);c12=zeros(kmax,1);c13=zeros(kmax,1);c14=zeros(kmax,1);c15=zeros(kmax,1);c16=zeros(kmax,1);
c21=zeros(kmax,1);c22=zeros(kmax,1);c23=zeros(kmax,1);c24=zeros(kmax,1);c25=zeros(kmax,1);c26=zeros(kmax,1);
c31=zeros(kmax,1);c32=zeros(kmax,1);c33=zeros(kmax,1);c34=zeros(kmax,1);c35=zeros(kmax,1);c36=zeros(kmax,1);
c41=zeros(kmax,1);c42=zeros(kmax,1);c43=zeros(kmax,1);c44=zeros(kmax,1);c45=zeros(kmax,1);c46=zeros(kmax,1);
c51=zeros(kmax,1);c52=zeros(kmax,1);c53=zeros(kmax,1);c54=zeros(kmax,1);c55=zeros(kmax,1);c56=zeros(kmax,1);
c61=zeros(kmax,1);c62=zeros(kmax,1);c63=zeros(kmax,1);c64=zeros(kmax,1);c65=zeros(kmax,1);c66=zeros(kmax,1);

%% c_{i0}(t) in (19c)
c1=zeros(kmax,1);c2=zeros(kmax,1);c3=zeros(kmax,1);c4=zeros(kmax,1);c5=zeros(kmax,1);c6=zeros(kmax,1); 

%% trigger flag of ETM-a for edge (i,j), if trigger triISij=1, else triISij=-1
triIS11=zeros(kmax,1);
triIS12=zeros(kmax,1);
triIS13=zeros(kmax,1);
triIS14=zeros(kmax,1);
triIS15=zeros(kmax,1);
triIS16=zeros(kmax,1);

triIS21=zeros(kmax,1);
triIS22=zeros(kmax,1);
triIS23=zeros(kmax,1);
triIS24=zeros(kmax,1);
triIS25=zeros(kmax,1);
triIS26=zeros(kmax,1);

triIS31=zeros(kmax,1);
triIS32=zeros(kmax,1);
triIS33=zeros(kmax,1);
triIS34=zeros(kmax,1);
triIS35=zeros(kmax,1);
triIS36=zeros(kmax,1);

triIS41=zeros(kmax,1);
triIS42=zeros(kmax,1);
triIS43=zeros(kmax,1);
triIS44=zeros(kmax,1);
triIS45=zeros(kmax,1);
triIS46=zeros(kmax,1);

triIS51=zeros(kmax,1);
triIS52=zeros(kmax,1);
triIS53=zeros(kmax,1);
triIS54=zeros(kmax,1);
triIS55=zeros(kmax,1);
triIS56=zeros(kmax,1);

triIS61=zeros(kmax,1);
triIS62=zeros(kmax,1);
triIS63=zeros(kmax,1);
triIS64=zeros(kmax,1);
triIS65=zeros(kmax,1);
triIS66=zeros(kmax,1);

%% trigger flag of ETM-b for edge (i,0), if trigger triISi0=1, else triISi0=-1
triIS10=zeros(kmax,1);
triIS20=zeros(kmax,1);
triIS30=zeros(kmax,1);
triIS40=zeros(kmax,1);
triIS50=zeros(kmax,1);
triIS60=zeros(kmax,1);

%% trigger number of ETM-a for edge (i,j) 
tricount11=0;
tricount12=0;
tricount13=0;
tricount14=0;
tricount15=0;
tricount16=0;

tricount21=0;
tricount22=0;
tricount23=0;
tricount24=0;
tricount25=0;
tricount26=0;

tricount31=0;
tricount32=0;
tricount33=0;
tricount34=0;
tricount35=0;
tricount36=0;

tricount41=0;
tricount42=0;
tricount43=0;
tricount44=0;
tricount45=0;
tricount46=0;

tricount51=0;
tricount52=0;
tricount53=0;
tricount54=0;
tricount55=0;
tricount56=0;

tricount61=0;
tricount62=0;
tricount63=0;
tricount64=0;
tricount65=0;
tricount66=0;

%% trigger number of ETM-b for edge (i,0) 
tricount10=0;
tricount20=0;
tricount30=0;
tricount40=0;
tricount50=0;
tricount60=0;

for k=1:kmax  
     if k==1       
        load date %% load the initial values of the variable x_{i}(t), \psi_{i}(t), \vartheta_{i}(t), c_{ij}(t) and c_{i0}(t)
       
       %% Initialization of variables at the initial time
       %% x_{i}(t) in (1) and (2)
        x1(2*k-1:2*k)=aa1;
        x2(2*k-1:2*k)=aa2;
        x3(3*k-2:3*k)=aa3;
        x4(3*k-2:3*k)=aa4;
        x5(4*k-3:4*k)=aa5;
        x6(4*k-3:4*k)=aa6;
        
        x0(2*k-1:2*k)=aa7;
        
       %% state observer \psi_{i}(t) in (27b)
        eita1(2*k-1:2*k)=aa8;
        eita2(2*k-1:2*k)=aa9;
        eita3(3*k-2:3*k)=aa10;
        eita4(3*k-2:3*k)=aa11;
        eita5(4*k-3:4*k)=aa12;
        eita6(4*k-3:4*k)=aa13;
        
       %% compensator \vartheta_{i}(t) in (19a) 
        kesi1(2*k-1:2*k)=aa14;
        kesi2(2*k-1:2*k)=aa15;
        kesi3(2*k-1:2*k)=aa16;
        kesi4(2*k-1:2*k)=aa17;
        kesi5(2*k-1:2*k)=aa18;
        kesi6(2*k-1:2*k)=aa19;
        
       %% c_{ij}(t) in (19b)
        c11(k)=aa20; c12(k)=aa21; c13(k)=aa22; c14(k)=aa23; c15(k)=aa24; c16(k)=aa25;
        c21(k)=aa26; c22(k)=aa27; c23(k)=aa28; c24(k)=aa29; c25(k)=aa30; c26(k)=aa31;
        c31(k)=aa32; c32(k)=aa33; c33(k)=aa34; c34(k)=aa35; c35(k)=aa36; c36(k)=aa37;
        c41(k)=aa38; c42(k)=aa39; c43(k)=aa40; c44(k)=aa41; c45(k)=aa42; c46(k)=aa43;
        c51(k)=aa44; c52(k)=aa45; c53(k)=aa46; c54(k)=aa47; c55(k)=aa48; c56(k)=aa49;
        c61(k)=aa50; c62(k)=aa51; c63(k)=aa52; c64(k)=aa53; c65(k)=aa54; c66(k)=aa55;
        
       %% c_{i0}(t) in (19c)
        c1(k)=aa56; c2(k)=aa57; c3(k)=aa58; c4(k)=aa59; c5(k)=aa60; c6(k)=aa61;
        
       %% u_{i}(t) in (27a) 
        u1(k)=K11*eita1(2*k-1:2*k)+K21*kesi1(2*k-1:2*k);
        u2(k)=K12*eita2(2*k-1:2*k)+K22*kesi2(2*k-1:2*k);
        u3(k)=K13*eita3(3*k-2:3*k)+K23*kesi3(2*k-1:2*k);
        u4(k)=K14*eita4(3*k-2:3*k)+K24*kesi4(2*k-1:2*k);
        u5(k)=K15*eita5(4*k-3:4*k)+K25*kesi5(2*k-1:2*k);
        u6(k)=K16*eita6(4*k-3:4*k)+K26*kesi6(2*k-1:2*k);
        
       %% \vartheta_{i}(t_{\sigma}^{ij})
        kesi11last(2*k-1:2*k)=kesi1(2*k-1:2*k);
        kesi12last(2*k-1:2*k)=kesi1(2*k-1:2*k);
        kesi13last(2*k-1:2*k)=kesi1(2*k-1:2*k);
        kesi14last(2*k-1:2*k)=kesi1(2*k-1:2*k);
        kesi15last(2*k-1:2*k)=kesi1(2*k-1:2*k);   
        kesi16last(2*k-1:2*k)=kesi1(2*k-1:2*k);       
        
        kesi21last(2*k-1:2*k)=kesi2(2*k-1:2*k);
        kesi22last(2*k-1:2*k)=kesi2(2*k-1:2*k);
        kesi23last(2*k-1:2*k)=kesi2(2*k-1:2*k);
        kesi24last(2*k-1:2*k)=kesi2(2*k-1:2*k);
        kesi25last(2*k-1:2*k)=kesi2(2*k-1:2*k);   
        kesi26last(2*k-1:2*k)=kesi2(2*k-1:2*k);    
        
        kesi31last(2*k-1:2*k)=kesi3(2*k-1:2*k);
        kesi32last(2*k-1:2*k)=kesi3(2*k-1:2*k);
        kesi33last(2*k-1:2*k)=kesi3(2*k-1:2*k);
        kesi34last(2*k-1:2*k)=kesi3(2*k-1:2*k);
        kesi35last(2*k-1:2*k)=kesi3(2*k-1:2*k);   
        kesi36last(2*k-1:2*k)=kesi3(2*k-1:2*k); 
        
        kesi41last(2*k-1:2*k)=kesi4(2*k-1:2*k);
        kesi42last(2*k-1:2*k)=kesi4(2*k-1:2*k);
        kesi43last(2*k-1:2*k)=kesi4(2*k-1:2*k);
        kesi44last(2*k-1:2*k)=kesi4(2*k-1:2*k);
        kesi45last(2*k-1:2*k)=kesi4(2*k-1:2*k);   
        kesi46last(2*k-1:2*k)=kesi4(2*k-1:2*k); 
        
        kesi51last(2*k-1:2*k)=kesi5(2*k-1:2*k);
        kesi52last(2*k-1:2*k)=kesi5(2*k-1:2*k);
        kesi53last(2*k-1:2*k)=kesi5(2*k-1:2*k);
        kesi54last(2*k-1:2*k)=kesi5(2*k-1:2*k);
        kesi55last(2*k-1:2*k)=kesi5(2*k-1:2*k);   
        kesi56last(2*k-1:2*k)=kesi5(2*k-1:2*k); 
        
        kesi61last(2*k-1:2*k)=kesi6(2*k-1:2*k);
        kesi62last(2*k-1:2*k)=kesi6(2*k-1:2*k);
        kesi63last(2*k-1:2*k)=kesi6(2*k-1:2*k);
        kesi64last(2*k-1:2*k)=kesi6(2*k-1:2*k);
        kesi65last(2*k-1:2*k)=kesi6 (2*k-1:2*k);   
        kesi66last(2*k-1:2*k)=kesi6(2*k-1:2*k); 
        
       %% \vartheta_{i}(t_{\sigma}^{i0})
        kesi10last(2*k-1:2*k)=kesi1(2*k-1:2*k);
        kesi20last(2*k-1:2*k)=kesi2(2*k-1:2*k);
        kesi30last(2*k-1:2*k)=kesi3(2*k-1:2*k);
        kesi40last(2*k-1:2*k)=kesi4(2*k-1:2*k);
        kesi50last(2*k-1:2*k)=kesi5(2*k-1:2*k);   
        kesi60last(2*k-1:2*k)=kesi6(2*k-1:2*k); 
        
       %% y_{i}(t) in (1) 
        y1(k)=C1*x1(2*k-1:2*k);
        y2(k)=C2*x2(2*k-1:2*k);
        y3(k)=C3*x3(3*k-2:3*k);
        y4(k)=C4*x4(3*k-2:3*k);
        y5(k)=C5*x5(4*k-3:4*k);
        y6(k)=C6*x6(4*k-3:4*k);
        y0(k)=C0*x0(2*k-1:2*k);
          
       %% The iteration value of the variables at the next time
       %% x_{i}(t) in (1) and (2)
        x1(2*(k+1)-1:2*(k+1))=x1(2*k-1:2*k)+T*(A1*x1(2*k-1:2*k)+B1*u1(k));
        x2(2*(k+1)-1:2*(k+1))=x2(2*k-1:2*k)+T*(A2*x2(2*k-1:2*k)+B2*u2(k));
        x3(3*(k+1)-2:3*(k+1))=x3(3*k-2:3*k)+T*(A3*x3(3*k-2:3*k)+B3*u3(k));
        x4(3*(k+1)-2:3*(k+1))=x4(3*k-2:3*k)+T*(A4*x4(3*k-2:3*k)+B4*u4(k));
        x5(4*(k+1)-3:4*(k+1))=x5(4*k-3:4*k)+T*(A5*x5(4*k-3:4*k)+B5*u5(k));
        x6(4*(k+1)-3:4*(k+1))=x6(4*k-3:4*k)+T*(A6*x6(4*k-3:4*k)+B6*u6(k));
        x0(2*(k+1)-1:2*(k+1))=x0(2*k-1:2*k)+T*(A0*x0(2*k-1:2*k));       
        
       %% triggering instant t_{\sigma}^{ij} of ETM-a for edge (i,j) 
        triT11=1;
        triT12=1;
        triT13=1;
        triT14=1;
        triT15=1;
        triT16=1;
        
        triT21=1;
        triT22=1;
        triT23=1;
        triT24=1;
        triT25=1;
        triT26=1;
        
        triT31=1;
        triT32=1;
        triT33=1;
        triT34=1;
        triT35=1;
        triT36=1;
        
        triT41=1;
        triT42=1;
        triT43=1;
        triT44=1;
        triT45=1;
        triT46=1;
        
        triT51=1;
        triT52=1;
        triT53=1;
        triT54=1;
        triT55=1;
        triT56=1;
        
        triT61=1;
        triT62=1;
        triT63=1;
        triT64=1;
        triT65=1;
        triT66=1;
        
       %% triggering instant t_{\sigma}^{i0} of ETM-b for edge (i,0) 
        triT10=1;
        triT20=1;
        triT30=1;
        triT40=1;
        triT50=1;
        triT60=1;
        
       %% e^{S(t-t_{\sigma}^{ij})}
        GTstep_S11=expm(A0*(k-triT11)*T);
        GTstep_S12=expm(A0*(k-triT12)*T);
        GTstep_S13=expm(A0*(k-triT13)*T);
        GTstep_S14=expm(A0*(k-triT14)*T);
        GTstep_S15=expm(A0*(k-triT15)*T);
        GTstep_S16=expm(A0*(k-triT16)*T);
        
        GTstep_S21=expm(A0*(k-triT21)*T);
        GTstep_S22=expm(A0*(k-triT22)*T);
        GTstep_S23=expm(A0*(k-triT23)*T);
        GTstep_S24=expm(A0*(k-triT24)*T);
        GTstep_S25=expm(A0*(k-triT25)*T);
        GTstep_S26=expm(A0*(k-triT26)*T);
        
        GTstep_S31=expm(A0*(k-triT31)*T);
        GTstep_S32=expm(A0*(k-triT32)*T);
        GTstep_S33=expm(A0*(k-triT33)*T);
        GTstep_S34=expm(A0*(k-triT34)*T);
        GTstep_S35=expm(A0*(k-triT35)*T);
        GTstep_S36=expm(A0*(k-triT36)*T);
        
        GTstep_S41=expm(A0*(k-triT41)*T);
        GTstep_S42=expm(A0*(k-triT42)*T);
        GTstep_S43=expm(A0*(k-triT43)*T);
        GTstep_S44=expm(A0*(k-triT44)*T);
        GTstep_S45=expm(A0*(k-triT45)*T);
        GTstep_S46=expm(A0*(k-triT46)*T);
        
        GTstep_S51=expm(A0*(k-triT51)*T);
        GTstep_S52=expm(A0*(k-triT52)*T);
        GTstep_S53=expm(A0*(k-triT53)*T);
        GTstep_S54=expm(A0*(k-triT54)*T);
        GTstep_S55=expm(A0*(k-triT55)*T);
        GTstep_S56=expm(A0*(k-triT56)*T);
        
        GTstep_S61=expm(A0*(k-triT61)*T);
        GTstep_S62=expm(A0*(k-triT62)*T);
        GTstep_S63=expm(A0*(k-triT63)*T);
        GTstep_S64=expm(A0*(k-triT64)*T);
        GTstep_S65=expm(A0*(k-triT65)*T);
        GTstep_S66=expm(A0*(k-triT66)*T);
        
       %% e^{S(t-t_{\sigma}^{i0})}
        GTstep_S10=expm(A0*(k-triT10)*T);
        GTstep_S20=expm(A0*(k-triT20)*T);
        GTstep_S30=expm(A0*(k-triT30)*T);
        GTstep_S40=expm(A0*(k-triT40)*T);
        GTstep_S50=expm(A0*(k-triT50)*T);
        GTstep_S60=expm(A0*(k-triT60)*T);
        
       %% compensator \vartheta_{i}(t) in (19a)
        kesi1(2*(k+1)-1:2*(k+1))=kesi1(2*k-1:2*k)+T*(A0*kesi1(2*k-1:2*k))+T*K*(c11(k)*abs(huaA(1,1))*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k))+c12(k)*abs(huaA(1,2))*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k))+c13(k)*abs(huaA(1,3))*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k))+c14(k)*abs(huaA(1,4))*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k))+c15(k)*abs(huaA(1,5))*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k))+c16(k)*abs(huaA(1,6))*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k))+c1(k)*huaB(1,1)*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k)));
        kesi2(2*(k+1)-1:2*(k+1))=kesi2(2*k-1:2*k)+T*(A0*kesi2(2*k-1:2*k))+T*K*(c21(k)*abs(huaA(2,1))*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k))+c22(k)*abs(huaA(2,2))*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k))+c23(k)*abs(huaA(2,3))*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k))+c24(k)*abs(huaA(2,4))*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k))+c25(k)*abs(huaA(2,5))*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k))+c26(k)*abs(huaA(2,6))*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k))+c2(k)*huaB(2,2)*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k)));
        kesi3(2*(k+1)-1:2*(k+1))=kesi3(2*k-1:2*k)+T*(A0*kesi3(2*k-1:2*k))+T*K*(c31(k)*abs(huaA(3,1))*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k))+c32(k)*abs(huaA(3,2))*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k))+c33(k)*abs(huaA(3,3))*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k))+c34(k)*abs(huaA(3,4))*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k))+c35(k)*abs(huaA(3,5))*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k))+c36(k)*abs(huaA(3,6))*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k))+c3(k)*huaB(3,3)*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k)));
        kesi4(2*(k+1)-1:2*(k+1))=kesi4(2*k-1:2*k)+T*(A0*kesi4(2*k-1:2*k))+T*K*(c41(k)*abs(huaA(4,1))*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k))+c42(k)*abs(huaA(4,2))*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k))+c43(k)*abs(huaA(4,3))*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k))+c44(k)*abs(huaA(4,4))*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k))+c45(k)*abs(huaA(4,5))*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k))+c46(k)*abs(huaA(4,6))*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k))+c4(k)*huaB(4,4)*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k)));
        kesi5(2*(k+1)-1:2*(k+1))=kesi5(2*k-1:2*k)+T*(A0*kesi5(2*k-1:2*k))+T*K*(c51(k)*abs(huaA(5,1))*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k))+c52(k)*abs(huaA(5,2))*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k))+c53(k)*abs(huaA(5,3))*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k))+c54(k)*abs(huaA(5,4))*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k))+c55(k)*abs(huaA(5,5))*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k))+c56(k)*abs(huaA(5,6))*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k))+c5(k)*huaB(5,5)*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k)));
        kesi6(2*(k+1)-1:2*(k+1))=kesi6(2*k-1:2*k)+T*(A0*kesi6(2*k-1:2*k))+T*K*(c61(k)*abs(huaA(6,1))*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k))+c62(k)*abs(huaA(6,2))*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k))+c63(k)*abs(huaA(6,3))*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k))+c64(k)*abs(huaA(6,4))*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k))+c65(k)*abs(huaA(6,5))*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k))+c66(k)*abs(huaA(6,6))*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k))+c6(k)*huaB(6,6)*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k)));          
 
       %% state observer \psi_{i}(t) in (27b)
        eita1(2*(k+1)-1:2*(k+1))=eita1(2*k-1:2*k)+T*(A1*eita1(2*k-1:2*k))+T*(B1*u1(k))+T*F1*(C1*eita1(2*k-1:2*k)-y1(k));
        eita2(2*(k+1)-1:2*(k+1))=eita2(2*k-1:2*k)+T*(A2*eita2(2*k-1:2*k))+T*(B2*u2(k))+T*F2*(C2*eita2(2*k-1:2*k)-y2(k));
        eita3(3*(k+1)-2:3*(k+1))=eita3(3*k-2:3*k)+T*(A3*eita3(3*k-2:3*k))+T*(B3*u3(k))+T*F3*(C3*eita3(3*k-2:3*k)-y3(k));
        eita4(3*(k+1)-2:3*(k+1))=eita4(3*k-2:3*k)+T*(A4*eita4(3*k-2:3*k))+T*(B4*u4(k))+T*F4*(C4*eita4(3*k-2:3*k)-y4(k));
        eita5(4*(k+1)-3:4*(k+1))=eita5(4*k-3:4*k)+T*(A5*eita5(4*k-3:4*k))+T*(B5*u5(k))+T*F5*(C5*eita5(4*k-3:4*k)-y5(k));
        eita6(4*(k+1)-3:4*(k+1))=eita6(4*k-3:4*k)+T*(A6*eita6(4*k-3:4*k))+T*(B6*u6(k))+T*F6*(C6*eita6(4*k-3:4*k)-y6(k));
        
       %% c_{ij}(t) in (19b)
        c11(k+1)= c11(k)+T*l11*(abs(huaA(1,1))*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k))'*Tao*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k)));
        c12(k+1)= c12(k)+T*l12*(abs(huaA(1,2))*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k))'*Tao*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k)));
        c13(k+1)= c13(k)+T*l13*(abs(huaA(1,3))*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k))'*Tao*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k)));
        c14(k+1)= c14(k)+T*l14*(abs(huaA(1,4))*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k))'*Tao*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k)));
        c15(k+1)= c15(k)+T*l15*(abs(huaA(1,5))*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k))'*Tao*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k)));
        c16(k+1)= c16(k)+T*l16*(abs(huaA(1,6))*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k))'*Tao*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k)));
        
        c21(k+1)= c21(k)+T*l21*(abs(huaA(2,1))*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k))'*Tao*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k)));
        c22(k+1)= c22(k)+T*l22*(abs(huaA(2,2))*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k))'*Tao*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k)));
        c23(k+1)= c23(k)+T*l23*(abs(huaA(2,3))*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k))'*Tao*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k)));
        c24(k+1)= c24(k)+T*l24*(abs(huaA(2,4))*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k))'*Tao*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k)));
        c25(k+1)= c25(k)+T*l25*(abs(huaA(2,5))*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k))'*Tao*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k)));
        c26(k+1)= c26(k)+T*l26*(abs(huaA(2,6))*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k))'*Tao*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k)));
       
        c31(k+1)= c31(k)+T*l31*(abs(huaA(3,1))*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k))'*Tao*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k)));
        c32(k+1)= c32(k)+T*l32*(abs(huaA(3,2))*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k))'*Tao*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k)));
        c33(k+1)= c33(k)+T*l33*(abs(huaA(3,3))*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k))'*Tao*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k)));
        c34(k+1)= c34(k)+T*l34*(abs(huaA(3,4))*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k))'*Tao*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k)));
        c35(k+1)= c35(k)+T*l35*(abs(huaA(3,5))*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k))'*Tao*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k)));
        c36(k+1)= c36(k)+T*l36*(abs(huaA(3,6))*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k))'*Tao*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k)));
        
        c41(k+1)= c41(k)+T*l41*(abs(huaA(4,1))*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k))'*Tao*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k)));
        c42(k+1)= c42(k)+T*l42*(abs(huaA(4,2))*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k))'*Tao*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k)));
        c43(k+1)= c43(k)+T*l43*(abs(huaA(4,3))*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k))'*Tao*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k)));
        c44(k+1)= c44(k)+T*l44*(abs(huaA(4,4))*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k))'*Tao*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k)));
        c45(k+1)= c45(k)+T*l45*(abs(huaA(4,5))*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k))'*Tao*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k)));
        c46(k+1)= c46(k)+T*l46*(abs(huaA(4,6))*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k))'*Tao*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k)));
        
        c51(k+1)= c51(k)+T*l51*(abs(huaA(5,1))*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k))'*Tao*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k)));
        c52(k+1)= c52(k)+T*l52*(abs(huaA(5,2))*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k))'*Tao*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k)));
        c53(k+1)= c53(k)+T*l53*(abs(huaA(5,3))*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k))'*Tao*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k)));
        c54(k+1)= c54(k)+T*l54*(abs(huaA(5,4))*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k))'*Tao*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k)));
        c55(k+1)= c55(k)+T*l55*(abs(huaA(5,5))*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k))'*Tao*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k)));
        c56(k+1)= c56(k)+T*l56*(abs(huaA(5,6))*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k))'*Tao*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k)));
        
        c61(k+1)= c61(k)+T*l61*(abs(huaA(6,1))*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k))'*Tao*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k)));
        c62(k+1)= c62(k)+T*l62*(abs(huaA(6,2))*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k))'*Tao*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k)));
        c63(k+1)= c63(k)+T*l63*(abs(huaA(6,3))*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k))'*Tao*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k)));
        c64(k+1)= c64(k)+T*l64*(abs(huaA(6,4))*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k))'*Tao*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k)));
        c65(k+1)= c65(k)+T*l65*(abs(huaA(6,5))*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k))'*Tao*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k)));
        c66(k+1)= c66(k)+T*l66*(abs(huaA(6,6))*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k))'*Tao*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k)));
        
       %% c_{i0}(t) in (19c)
        c1(k+1)= c1(k)+T*l1*(huaB(1,1)*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k))'*Tao*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k)));
        c2(k+1)= c2(k)+T*l2*(huaB(2,2)*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k))'*Tao*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k)));
        c3(k+1)= c3(k)+T*l3*(huaB(3,3)*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k))'*Tao*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k)));
        c4(k+1)= c4(k)+T*l4*(huaB(4,4)*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k))'*Tao*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k)));
        c5(k+1)= c5(k)+T*l5*(huaB(5,5)*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k))'*Tao*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k)));
        c6(k+1)= c6(k)+T*l6*(huaB(6,6)*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k))'*Tao*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k)));
        
       %% e_{ij}(t) in (5a)
        e11(2*k-1:2*k)=GTstep_S11*kesi11last(2*k-1:2*k)-kesi1(2*k-1:2*k);
        e12(2*k-1:2*k)=GTstep_S12*kesi12last(2*k-1:2*k)-kesi1(2*k-1:2*k);
        e13(2*k-1:2*k)=GTstep_S13*kesi13last(2*k-1:2*k)-kesi1(2*k-1:2*k);
        e14(2*k-1:2*k)=GTstep_S14*kesi14last(2*k-1:2*k)-kesi1(2*k-1:2*k);
        e15(2*k-1:2*k)=GTstep_S15*kesi15last(2*k-1:2*k)-kesi1(2*k-1:2*k);  
        e16(2*k-1:2*k)=GTstep_S16*kesi16last(2*k-1:2*k)-kesi1(2*k-1:2*k);  
        
        e21(2*k-1:2*k)=GTstep_S21*kesi21last(2*k-1:2*k)-kesi2(2*k-1:2*k);
        e22(2*k-1:2*k)=GTstep_S22*kesi22last(2*k-1:2*k)-kesi2(2*k-1:2*k);
        e23(2*k-1:2*k)=GTstep_S23*kesi23last(2*k-1:2*k)-kesi2(2*k-1:2*k);
        e24(2*k-1:2*k)=GTstep_S24*kesi24last(2*k-1:2*k)-kesi2(2*k-1:2*k);
        e25(2*k-1:2*k)=GTstep_S25*kesi25last(2*k-1:2*k)-kesi2(2*k-1:2*k);  
        e26(2*k-1:2*k)=GTstep_S26*kesi26last(2*k-1:2*k)-kesi2(2*k-1:2*k); 
        
        e31(2*k-1:2*k)=GTstep_S31*kesi31last(2*k-1:2*k)-kesi3(2*k-1:2*k);
        e32(2*k-1:2*k)=GTstep_S32*kesi32last(2*k-1:2*k)-kesi3(2*k-1:2*k);
        e33(2*k-1:2*k)=GTstep_S33*kesi33last(2*k-1:2*k)-kesi3(2*k-1:2*k);
        e34(2*k-1:2*k)=GTstep_S34*kesi34last(2*k-1:2*k)-kesi3(2*k-1:2*k);
        e35(2*k-1:2*k)=GTstep_S35*kesi35last(2*k-1:2*k)-kesi3(2*k-1:2*k);  
        e36(2*k-1:2*k)=GTstep_S36*kesi36last(2*k-1:2*k)-kesi3(2*k-1:2*k); 
        
        e41(2*k-1:2*k)=GTstep_S41*kesi41last(2*k-1:2*k)-kesi4(2*k-1:2*k);
        e42(2*k-1:2*k)=GTstep_S42*kesi42last(2*k-1:2*k)-kesi4(2*k-1:2*k);
        e43(2*k-1:2*k)=GTstep_S43*kesi43last(2*k-1:2*k)-kesi4(2*k-1:2*k);
        e44(2*k-1:2*k)=GTstep_S44*kesi44last(2*k-1:2*k)-kesi4(2*k-1:2*k);
        e45(2*k-1:2*k)=GTstep_S45*kesi45last(2*k-1:2*k)-kesi4(2*k-1:2*k);  
        e46(2*k-1:2*k)=GTstep_S46*kesi46last(2*k-1:2*k)-kesi4(2*k-1:2*k); 
        
        e51(2*k-1:2*k)=GTstep_S51*kesi51last(2*k-1:2*k)-kesi5(2*k-1:2*k);
        e52(2*k-1:2*k)=GTstep_S52*kesi52last(2*k-1:2*k)-kesi5(2*k-1:2*k);
        e53(2*k-1:2*k)=GTstep_S53*kesi53last(2*k-1:2*k)-kesi5(2*k-1:2*k);
        e54(2*k-1:2*k)=GTstep_S54*kesi54last(2*k-1:2*k)-kesi5(2*k-1:2*k);
        e55(2*k-1:2*k)=GTstep_S55*kesi55last(2*k-1:2*k)-kesi5(2*k-1:2*k);  
        e56(2*k-1:2*k)=GTstep_S56*kesi56last(2*k-1:2*k)-kesi5(2*k-1:2*k); 
        
        e61(2*k-1:2*k)=GTstep_S61*kesi61last(2*k-1:2*k)-kesi6(2*k-1:2*k);
        e62(2*k-1:2*k)=GTstep_S62*kesi62last(2*k-1:2*k)-kesi6(2*k-1:2*k);
        e63(2*k-1:2*k)=GTstep_S63*kesi63last(2*k-1:2*k)-kesi6(2*k-1:2*k);
        e64(2*k-1:2*k)=GTstep_S64*kesi64last(2*k-1:2*k)-kesi6(2*k-1:2*k);
        e65(2*k-1:2*k)=GTstep_S65*kesi65last(2*k-1:2*k)-kesi6(2*k-1:2*k);  
        e66(2*k-1:2*k)=GTstep_S66*kesi66last(2*k-1:2*k)-kesi6(2*k-1:2*k); 
        
       %% e_{i0}(t) in (5b)
        e10(2*k-1:2*k)=GTstep_S10*kesi10last(2*k-1:2*k)-kesi1(2*k-1:2*k);
        e20(2*k-1:2*k)=GTstep_S20*kesi20last(2*k-1:2*k)-kesi2(2*k-1:2*k);
        e30(2*k-1:2*k)=GTstep_S30*kesi30last(2*k-1:2*k)-kesi3(2*k-1:2*k);
        e40(2*k-1:2*k)=GTstep_S40*kesi40last(2*k-1:2*k)-kesi4(2*k-1:2*k);
        e50(2*k-1:2*k)=GTstep_S50*kesi50last(2*k-1:2*k)-kesi5(2*k-1:2*k);  
        e60(2*k-1:2*k)=GTstep_S60*kesi60last(2*k-1:2*k)-kesi6(2*k-1:2*k); 
        
       %% f_{ij}(t) in (23)
        f111=(1+2*deta*c11(k))*e11(2*k-1:2*k)'*Tao*e11(2*k-1:2*k);
        f112=(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k))'*Tao*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k));
        f113=mu*exp(-v*(k-1)*T);
        f11(k)=f111-1/4*f112-f113;

        f121=(1+2*deta*c12(k))*e12(2*k-1:2*k)'*Tao*e12(2*k-1:2*k);
        f122=(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k))'*Tao*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k));
        f123=mu*exp(-v*(k-1)*T);
        f12(k)=f121-1/4*f122-f123;

        f131=(1+2*deta*c13(k))*e13(2*k-1:2*k)'*Tao*e13(2*k-1:2*k);
        f132=(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k))'*Tao*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k));
        f133=mu*exp(-v*(k-1)*T);
        f13(k)=f131-1/4*f132-f133;
        
        f141=(1+2*deta*c14(k))*e14(2*k-1:2*k)'*Tao*e14(2*k-1:2*k);
        f142=(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k))'*Tao*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k));
        f143=mu*exp(-v*(k-1)*T);
        f14(k)=f141-1/4*f142-f143;
        
        f151=(1+2*deta*c15(k))*e15(2*k-1:2*k)'*Tao*e15(2*k-1:2*k);
        f152=(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k))'*Tao*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k));
        f153=mu*exp(-v*(k-1)*T);
        f15(k)=f151-1/4*f152-f153;
        
        f161=(1+2*deta*c16(k))*e16(2*k-1:2*k)'*Tao*e16(2*k-1:2*k);
        f162=(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k))'*Tao*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k));
        f163=mu*exp(-v*(k-1)*T);
        f16(k)=f161-1/4*f162-f163;
        
        f211=(1+2*deta*c21(k))*e21(2*k-1:2*k)'*Tao*e21(2*k-1:2*k);
        f212=(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k))'*Tao*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k));
        f213=mu*exp(-v*(k-1)*T);
        f21(k)=f211-1/4*f212-f213;
        
        f221=(1+2*deta*c22(k))*e22(2*k-1:2*k)'*Tao*e22(2*k-1:2*k);
        f222=(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k))'*Tao*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k));
        f223=mu*exp(-v*(k-1)*T);
        f22(k)=f221-1/4*f222-f223;
        
        f231=(1+2*deta*c23(k))*e23(2*k-1:2*k)'*Tao*e23(2*k-1:2*k);
        f232=(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k))'*Tao*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k));
        f233=mu*exp(-v*(k-1)*T);
        f23(k)=f231-1/4*f232-f233;
        
        f241=(1+2*deta*c24(k))*e24(2*k-1:2*k)'*Tao*e24(2*k-1:2*k);
        f242=(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k))'*Tao*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k));
        f243=mu*exp(-v*(k-1)*T);
        f24(k)=f241-1/4*f242-f243;
        
        f251=(1+2*deta*c25(k))*e25(2*k-1:2*k)'*Tao*e25(2*k-1:2*k);
        f252=(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k))'*Tao*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k));
        f253=mu*exp(-v*(k-1)*T);
        f25(k)=f251-1/4*f252-f253;
        
        f261=(1+2*deta*c26(k))*e26(2*k-1:2*k)'*Tao*e26(2*k-1:2*k);
        f262=(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k))'*Tao*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k));
        f263=mu*exp(-v*(k-1)*T);
        f26(k)=f261-1/4*f262-f263;
        
        f311=(1+2*deta*c31(k))*e31(2*k-1:2*k)'*Tao*e31(2*k-1:2*k);
        f312=(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k))'*Tao*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k));
        f313=mu*exp(-v*(k-1)*T);
        f31(k)=f311-1/4*f312-f313;
        
        f321=(1+2*deta*c32(k))*e32(2*k-1:2*k)'*Tao*e32(2*k-1:2*k);
        f322=(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k))'*Tao*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k));
        f323=mu*exp(-v*(k-1)*T);
        f32(k)=f321-1/4*f322-f323;
        
        f331=(1+2*deta*c33(k))*e33(2*k-1:2*k)'*Tao*e33(2*k-1:2*k);
        f332=(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k))'*Tao*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k));
        f333=mu*exp(-v*(k-1)*T);
        f33(k)=f331-1/4*f332-f333;
        
        f341=(1+2*deta*c34(k))*e34(2*k-1:2*k)'*Tao*e34(2*k-1:2*k);
        f342=(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k))'*Tao*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k));
        f343=mu*exp(-v*(k-1)*T);
        f34(k)=f341-1/4*f342-f343;
        
        f351=(1+2*deta*c35(k))*e35(2*k-1:2*k)'*Tao*e35(2*k-1:2*k);
        f352=(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k))'*Tao*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k));
        f353=mu*exp(-v*(k-1)*T);
        f35(k)=f351-1/4*f352-f353;
        
        f361=(1+2*deta*c36(k))*e36(2*k-1:2*k)'*Tao*e36(2*k-1:2*k);
        f362=(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k))'*Tao*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k));
        f363=mu*exp(-v*(k-1)*T);
        f36(k)=f361-1/4*f362-f363;
        
        f411=(1+2*deta*c41(k))*e41(2*k-1:2*k)'*Tao*e41(2*k-1:2*k);
        f412=(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k))'*Tao*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k));
        f413=mu*exp(-v*(k-1)*T);
        f41(k)=f411-1/4*f412-f413;
        
        f421=(1+2*deta*c42(k))*e42(2*k-1:2*k)'*Tao*e42(2*k-1:2*k);
        f422=(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k))'*Tao*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k));
        f423=mu*exp(-v*(k-1)*T);
        f42(k)=f421-1/4*f422-f423;
        
        f431=(1+2*deta*c43(k))*e43(2*k-1:2*k)'*Tao*e43(2*k-1:2*k);
        f432=(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k))'*Tao*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k));
        f433=mu*exp(-v*(k-1)*T);
        f43(k)=f431-1/4*f432-f433;
        
        f441=(1+2*deta*c44(k))*e44(2*k-1:2*k)'*Tao*e44(2*k-1:2*k);
        f442=(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k))'*Tao*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k));
        f443=mu*exp(-v*(k-1)*T);
        f44(k)=f441-1/4*f442-f443;
        
        f451=(1+2*deta*c45(k))*e45(2*k-1:2*k)'*Tao*e45(2*k-1:2*k);
        f452=(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k))'*Tao*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k));
        f453=mu*exp(-v*(k-1)*T);
        f45(k)=f451-1/4*f452-f453;
        
        f461=(1+2*deta*c46(k))*e46(2*k-1:2*k)'*Tao*e46(2*k-1:2*k);
        f462=(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k))'*Tao*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k));
        f463=mu*exp(-v*(k-1)*T);
        f46(k)=f461-1/4*f462-f463;
        
        f511=(1+2*deta*c51(k))*e51(2*k-1:2*k)'*Tao*e51(2*k-1:2*k);
        f512=(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k))'*Tao*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k));
        f513=mu*exp(-v*(k-1)*T);
        f51(k)=f511-1/4*f512-f513;
        
        f521=(1+2*deta*c52(k))*e52(2*k-1:2*k)'*Tao*e52(2*k-1:2*k);
        f522=(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k))'*Tao*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k));
        f523=mu*exp(-v*(k-1)*T);
        f52(k)=f521-1/4*f522-f523;
        
        f531=(1+2*deta*c53(k))*e53(2*k-1:2*k)'*Tao*e53(2*k-1:2*k);
        f532=(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k))'*Tao*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k));
        f533=mu*exp(-v*(k-1)*T);
        f53(k)=f531-1/4*f532-f533;
        
        f541=(1+2*deta*c54(k))*e54(2*k-1:2*k)'*Tao*e54(2*k-1:2*k);
        f542=(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k))'*Tao*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k));
        f543=mu*exp(-v*(k-1)*T);
        f54(k)=f541-1/4*f542-f543;
        
        f551=(1+2*deta*c55(k))*e55(2*k-1:2*k)'*Tao*e55(2*k-1:2*k);
        f552=(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k))'*Tao*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k));
        f553=mu*exp(-v*(k-1)*T);
        f55(k)=f551-1/4*f552-f553;
        
        f561=(1+2*deta*c56(k))*e56(2*k-1:2*k)'*Tao*e56(2*k-1:2*k);
        f562=(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k))'*Tao*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k));
        f563=mu*exp(-v*(k-1)*T);
        f56(k)=f561-1/4*f562-f563;
        
        f611=(1+2*deta*c61(k))*e61(2*k-1:2*k)'*Tao*e61(2*k-1:2*k);
        f612=(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k))'*Tao*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k));
        f613=mu*exp(-v*(k-1)*T);
        f61(k)=f611-1/4*f612-f613;
        
        f621=(1+2*deta*c62(k))*e62(2*k-1:2*k)'*Tao*e62(2*k-1:2*k);
        f622=(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k))'*Tao*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k));
        f623=mu*exp(-v*(k-1)*T);
        f62(k)=f621-1/4*f622-f623;
        
        f631=(1+2*deta*c63(k))*e63(2*k-1:2*k)'*Tao*e63(2*k-1:2*k);
        f632=(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k))'*Tao*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k));
        f633=mu*exp(-v*(k-1)*T);
        f63(k)=f631-1/4*f632-f633;
        
        f641=(1+2*deta*c64(k))*e64(2*k-1:2*k)'*Tao*e64(2*k-1:2*k);
        f642=(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k))'*Tao*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k));
        f643=mu*exp(-v*(k-1)*T);
        f64(k)=f641-1/4*f642-f643;
        
        f651=(1+2*deta*c65(k))*e65(2*k-1:2*k)'*Tao*e65(2*k-1:2*k);
        f652=(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k))'*Tao*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k));
        f653=mu*exp(-v*(k-1)*T);
        f65(k)=f651-1/4*f652-f653;
        
        f661=(1+2*deta*c66(k))*e66(2*k-1:2*k)'*Tao*e66(2*k-1:2*k);
        f662=(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k))'*Tao*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k));
        f663=mu*exp(-v*(k-1)*T);
        f66(k)=f661-1/4*f662-f663;
        
       %% f_{i0}(t) in (24)
        f101=(1+2*deta*c1(k))*e10(2*k-1:2*k)'*Tao*e10(2*k-1:2*k);
        f102=(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k))'*Tao*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k));
        f103=mu*exp(-v*(k-1)*T);
        f10(k)=f101-f102-f103;
        
        f201=(1+2*deta*c2(k))*e20(2*k-1:2*k)'*Tao*e20(2*k-1:2*k);
        f202=(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k))'*Tao*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k));
        f203=mu*exp(-v*(k-1)*T);
        f20(k)=f201-f202-f203;
        
        f301=(1+2*deta*c3(k))*e30(2*k-1:2*k)'*Tao*e30(2*k-1:2*k);
        f302=(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k))'*Tao*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k));
        f303=mu*exp(-v*(k-1)*T);
        f30(k)=f301-f302-f303;
        
        f401=(1+2*deta*c4(k))*e40(2*k-1:2*k)'*Tao*e40(2*k-1:2*k);
        f402=(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k))'*Tao*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k));
        f403=mu*exp(-v*(k-1)*T);
        f40(k)=f401-f402-f403;
        
        f501=(1+2*deta*c5(k))*e50(2*k-1:2*k)'*Tao*e50(2*k-1:2*k);
        f502=(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k))'*Tao*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k));
        f503=mu*exp(-v*(k-1)*T);
        f50(k)=f501-f502-f503;
        
        f601=(1+2*deta*c6(k))*e60(2*k-1:2*k)'*Tao*e60(2*k-1:2*k);
        f602=(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k))'*Tao*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k));
        f603=mu*exp(-v*(k-1)*T);
        f60(k)=f601-f602-f603;
        
       %% trigger flag of ETM-a for edge (i,j), if trigger triISij=1, else triISij=-1
        triIS11(k)=1;
        triIS12(k)=1;
        triIS13(k)=1;
        triIS14(k)=1;
        triIS15(k)=1;
        triIS16(k)=1;
        
        triIS21(k)=1;
        triIS22(k)=1;
        triIS23(k)=1;
        triIS24(k)=1;
        triIS25(k)=1;
        triIS26(k)=1;
        
        triIS31(k)=1;
        triIS32(k)=1;
        triIS33(k)=1;
        triIS34(k)=1;
        triIS35(k)=1;
        triIS36(k)=1;
        
        triIS41(k)=1;
        triIS42(k)=1;
        triIS43(k)=1;
        triIS44(k)=1;
        triIS45(k)=1;
        triIS46(k)=1;
        
        triIS51(k)=1;
        triIS52(k)=1;
        triIS53(k)=1;
        triIS54(k)=1;
        triIS55(k)=1;
        triIS56(k)=1;
        
        triIS61(k)=1;
        triIS62(k)=1;
        triIS63(k)=1;
        triIS64(k)=1;
        triIS65(k)=1;
        triIS66(k)=1;
        
       %% trigger flag of ETM-b for edge (i,0), if trigger triISi0=1, else triISi0=-1
        triIS10(k)=1;
        triIS20(k)=1;
        triIS30(k)=1;
        triIS40(k)=1;
        triIS50(k)=1;
        triIS60(k)=1;
        
       %% trigger number of ETM-a for edge (i,j) 
        tricount11=tricount11+1;
        tricount12=tricount12+1;
        tricount13=tricount13+1;
        tricount14=tricount14+1;
        tricount15=tricount15+1;
        tricount16=tricount16+1;
        
        tricount21=tricount21+1;
        tricount22=tricount22+1;
        tricount23=tricount23+1;
        tricount24=tricount24+1;
        tricount25=tricount25+1;
        tricount26=tricount26+1;
        
        tricount31=tricount31+1;
        tricount32=tricount32+1;
        tricount33=tricount33+1;
        tricount34=tricount34+1;
        tricount35=tricount35+1;
        tricount36=tricount36+1;
        
        tricount41=tricount41+1;
        tricount42=tricount42+1;
        tricount43=tricount43+1;
        tricount44=tricount44+1;
        tricount45=tricount45+1;
        tricount46=tricount46+1;
        
        tricount51=tricount51+1;
        tricount52=tricount52+1;
        tricount53=tricount53+1;
        tricount54=tricount54+1;
        tricount55=tricount55+1;
        tricount56=tricount56+1;
        
        tricount61=tricount61+1;
        tricount62=tricount62+1;
        tricount63=tricount63+1;
        tricount64=tricount64+1;
        tricount65=tricount65+1;
        tricount66=tricount66+1;
        
       %% trigger number of ETM-b for edge (i,0) 
        tricount10=tricount10+1;
        tricount20=tricount20+1;
        tricount30=tricount30+1;
        tricount40=tricount40+1;
        tricount50=tricount50+1;
        tricount60=tricount60+1;
                
    else if k<kmax
       %% e^{S(t-t_{\sigma}^{ij})}
        GTstep_S11=expm(A0*(k-triT11)*T);
        GTstep_S12=expm(A0*(k-triT12)*T);
        GTstep_S13=expm(A0*(k-triT13)*T);
        GTstep_S14=expm(A0*(k-triT14)*T);
        GTstep_S15=expm(A0*(k-triT15)*T);
        GTstep_S16=expm(A0*(k-triT16)*T);
        
        GTstep_S21=expm(A0*(k-triT21)*T);
        GTstep_S22=expm(A0*(k-triT22)*T);
        GTstep_S23=expm(A0*(k-triT23)*T);
        GTstep_S24=expm(A0*(k-triT24)*T);
        GTstep_S25=expm(A0*(k-triT25)*T);
        GTstep_S26=expm(A0*(k-triT26)*T);
        
        GTstep_S31=expm(A0*(k-triT31)*T);
        GTstep_S32=expm(A0*(k-triT32)*T);
        GTstep_S33=expm(A0*(k-triT33)*T);
        GTstep_S34=expm(A0*(k-triT34)*T);
        GTstep_S35=expm(A0*(k-triT35)*T);
        GTstep_S36=expm(A0*(k-triT36)*T);
        
        GTstep_S41=expm(A0*(k-triT41)*T);
        GTstep_S42=expm(A0*(k-triT42)*T);
        GTstep_S43=expm(A0*(k-triT43)*T);
        GTstep_S44=expm(A0*(k-triT44)*T);
        GTstep_S45=expm(A0*(k-triT45)*T);
        GTstep_S46=expm(A0*(k-triT46)*T);
        
        GTstep_S51=expm(A0*(k-triT51)*T);
        GTstep_S52=expm(A0*(k-triT52)*T);
        GTstep_S53=expm(A0*(k-triT53)*T);
        GTstep_S54=expm(A0*(k-triT54)*T);
        GTstep_S55=expm(A0*(k-triT55)*T);
        GTstep_S56=expm(A0*(k-triT56)*T);
        
        GTstep_S61=expm(A0*(k-triT61)*T);
        GTstep_S62=expm(A0*(k-triT62)*T);
        GTstep_S63=expm(A0*(k-triT63)*T);
        GTstep_S64=expm(A0*(k-triT64)*T);
        GTstep_S65=expm(A0*(k-triT65)*T);
        GTstep_S66=expm(A0*(k-triT66)*T);
        
       %% e^{S(t-t_{\sigma}^{i0})}
        GTstep_S10=expm(A0*(k-triT10)*T);
        GTstep_S20=expm(A0*(k-triT20)*T);
        GTstep_S30=expm(A0*(k-triT30)*T);
        GTstep_S40=expm(A0*(k-triT40)*T);
        GTstep_S50=expm(A0*(k-triT50)*T);
        GTstep_S60=expm(A0*(k-triT60)*T);
        
       %% judge whether the ETM-a for edge (i,j) is triggered
        e11(2*k-1:2*k)=GTstep_S11*kesi11last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);
        f111=(1+2*deta*c11(k))*e11(2*k-1:2*k)'*Tao*e11(2*k-1:2*k);
        f112=(GTstep_S11*kesi11last(2*(k-1)-1:2*(k-1))-sign(huaA(1,1))*GTstep_S11*kesi11last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S11*kesi11last(2*(k-1)-1:2*(k-1))-sign(huaA(1,1))*GTstep_S11*kesi11last(2*(k-1)-1:2*(k-1)));
        f113=mu*exp(-v*(k-1)*T);
        f11(k)=f111-1/4*f112-f113; %% f_{11}(t) in (23)    
        if f11(k)>0
                 kesi11last(2*k-1:2*k)=kesi1(2*k-1:2*k);  %% if tigger, \vartheta_{i}(t_{\sigma}^{ij}) is updated to \vartheta_{i}(t)
                 triT11=k;  %% triggering instant t_{\sigma}^{ij} of ETM-a for edge (i,j) 
                 triIS11(k)=1; %% trigger flag of ETM-a for edge (i,j), if trigger triISij=1, else triISij=-1
                 tricount11=tricount11+1; %% trigger number of ETM-a for edge (i,j) 
                
            else
               kesi11last(2*k-1:2*k)=kesi11last(2*(k-1)-1:2*(k-1)); %% if not tigger, \vartheta_{i}(t_{\sigma}^{ij}) remains unchanged
               triIS11(k)=-1; %% trigger flag of ETM-a for edge (i,j), if trigger triISij=1, else triISij=-1
        end
        
        e12(2*k-1:2*k)=GTstep_S12*kesi12last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);
        f121=(1+2*deta*c12(k))*e12(2*k-1:2*k)'*Tao*e12(2*k-1:2*k);
        f122=(GTstep_S12*kesi12last(2*(k-1)-1:2*(k-1))-sign(huaA(1,2))*GTstep_S21*kesi21last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S12*kesi12last(2*(k-1)-1:2*(k-1))-sign(huaA(1,2))*GTstep_S21*kesi21last(2*(k-1)-1:2*(k-1)));
        f123=mu*exp(-v*(k-1)*T);
        f12(k)=f121-1/4*f122-f123;
        if f12(k)>0
                 kesi12last(2*k-1:2*k)=kesi1(2*k-1:2*k);
                 triT12=k;
                 triIS12(k)=1;
                 tricount12=tricount12+1;
                
            else
               kesi12last(2*k-1:2*k)=kesi12last(2*(k-1)-1:2*(k-1));
               triIS12(k)=-1; 
        end
        
        e13(2*k-1:2*k)=GTstep_S13*kesi13last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);
        f131=(1+2*deta*c13(k))*e13(2*k-1:2*k)'*Tao*e13(2*k-1:2*k);
        f132=(GTstep_S13*kesi13last(2*(k-1)-1:2*(k-1))-sign(huaA(1,3))*GTstep_S31*kesi31last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S13*kesi13last(2*(k-1)-1:2*(k-1))-sign(huaA(1,3))*GTstep_S31*kesi31last(2*(k-1)-1:2*(k-1)));
        f133=mu*exp(-v*(k-1)*T);
        f13(k)=f131-1/4*f132-f133;
        if f13(k)>0
                 kesi13last(2*k-1:2*k)=kesi1(2*k-1:2*k);
                 triT13=k;
                 triIS13(k)=1;
                 tricount13=tricount13+1;
                
            else
               kesi13last(2*k-1:2*k)=kesi13last(2*(k-1)-1:2*(k-1));
               triIS13(k)=-1; 
        end
        
        e14(2*k-1:2*k)=GTstep_S14*kesi14last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);
        f141=(1+2*deta*c14(k))*e14(2*k-1:2*k)'*Tao*e14(2*k-1:2*k);
        f142=(GTstep_S14*kesi14last(2*(k-1)-1:2*(k-1))-sign(huaA(1,4))*GTstep_S41*kesi41last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S14*kesi14last(2*(k-1)-1:2*(k-1))-sign(huaA(1,4))*GTstep_S41*kesi41last(2*(k-1)-1:2*(k-1)));
        f143=mu*exp(-v*(k-1)*T);
        f14(k)=f141-1/4*f142-f143;
        if f14(k)>0
                 kesi14last(2*k-1:2*k)=kesi1(2*k-1:2*k);
                 triT14=k;
                 triIS14(k)=1;
                 tricount14=tricount14+1;
                
            else
               kesi14last(2*k-1:2*k)=kesi14last(2*(k-1)-1:2*(k-1));
               triIS14(k)=-1; 
        end
        
        e15(2*k-1:2*k)=GTstep_S15*kesi15last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);  
        f151=(1+2*deta*c15(k))*e15(2*k-1:2*k)'*Tao*e15(2*k-1:2*k);
        f152=(GTstep_S15*kesi15last(2*(k-1)-1:2*(k-1))-sign(huaA(1,5))*GTstep_S51*kesi51last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S15*kesi15last(2*(k-1)-1:2*(k-1))-sign(huaA(1,5))*GTstep_S51*kesi51last(2*(k-1)-1:2*(k-1)));
        f153=mu*exp(-v*(k-1)*T);
        f15(k)=f151-1/4*f152-f153;
        if f15(k)>0
                 kesi15last(2*k-1:2*k)=kesi1(2*k-1:2*k);
                 triT15=k;
                 triIS15(k)=1;
                 tricount15=tricount15+1;
                
            else
               kesi15last(2*k-1:2*k)=kesi15last(2*(k-1)-1:2*(k-1));
               triIS15(k)=-1; 
        end
        
        e16(2*k-1:2*k)=GTstep_S16*kesi16last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);
        f161=(1+2*deta*c16(k))*e16(2*k-1:2*k)'*Tao*e16(2*k-1:2*k);
        f162=(GTstep_S16*kesi16last(2*(k-1)-1:2*(k-1))-sign(huaA(1,6))*GTstep_S61*kesi61last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S16*kesi16last(2*(k-1)-1:2*(k-1))-sign(huaA(1,6))*GTstep_S61*kesi61last(2*(k-1)-1:2*(k-1)));
        f163=mu*exp(-v*(k-1)*T);
        f16(k)=f161-1/4*f162-f163;
        if f16(k)>0
                 kesi16last(2*k-1:2*k)=kesi1(2*k-1:2*k);
                 triT16=k;
                 triIS16(k)=1;
                 tricount16=tricount16+1;
                
            else
               kesi16last(2*k-1:2*k)=kesi16last(2*(k-1)-1:2*(k-1));
               triIS16(k)=-1; 
        end
        
        e21(2*k-1:2*k)=GTstep_S21*kesi21last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f211=(1+2*deta*c21(k))*e21(2*k-1:2*k)'*Tao*e21(2*k-1:2*k);
        f212=(GTstep_S21*kesi21last(2*(k-1)-1:2*(k-1))-sign(huaA(2,1))*GTstep_S12*kesi12last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S21*kesi21last(2*(k-1)-1:2*(k-1))-sign(huaA(2,1))*GTstep_S12*kesi12last(2*(k-1)-1:2*(k-1)));
        f213=mu*exp(-v*(k-1)*T);
        f21(k)=f211-1/4*f212-f213;
        if f21(k)>0
                 kesi21last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT21=k;
                 triIS21(k)=1;
                 tricount14=tricount14+1;
                
            else
               kesi21last(2*k-1:2*k)=kesi21last(2*(k-1)-1:2*(k-1));
               triIS21(k)=-1; 
        end
        
        e22(2*k-1:2*k)=GTstep_S22*kesi22last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f221=(1+2*deta*c22(k))*e22(2*k-1:2*k)'*Tao*e22(2*k-1:2*k);
        f222=(GTstep_S22*kesi22last(2*(k-1)-1:2*(k-1))-sign(huaA(2,2))*GTstep_S22*kesi22last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S22*kesi22last(2*(k-1)-1:2*(k-1))-sign(huaA(2,2))*GTstep_S22*kesi22last(2*(k-1)-1:2*(k-1)));
        f223=mu*exp(-v*(k-1)*T);
        f22(k)=f221-1/4*f222-f223;
        if f22(k)>0
                 kesi22last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT22=k;
                 triIS22(k)=1;
                 tricount22=tricount22+1;
                
            else
               kesi22last(2*k-1:2*k)=kesi22last(2*(k-1)-1:2*(k-1));
               triIS22(k)=-1; 
        end
        
        e23(2*k-1:2*k)=GTstep_S23*kesi23last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f231=(1+2*deta*c23(k))*e23(2*k-1:2*k)'*Tao*e23(2*k-1:2*k);
        f232=(GTstep_S23*kesi23last(2*(k-1)-1:2*(k-1))-sign(huaA(2,3))*GTstep_S32*kesi32last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S23*kesi23last(2*(k-1)-1:2*(k-1))-sign(huaA(2,3))*GTstep_S32*kesi32last(2*(k-1)-1:2*(k-1)));
        f233=mu*exp(-v*(k-1)*T);
        f23(k)=f231-1/4*f232-f233;
        if f23(k)>0
                 kesi23last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT23=k;
                 triIS23(k)=1;
                 tricount23=tricount23+1;
                
            else
               kesi23last(2*k-1:2*k)=kesi23last(2*(k-1)-1:2*(k-1));
               triIS23(k)=-1; 
        end
        
        e24(2*k-1:2*k)=GTstep_S24*kesi24last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f241=(1+2*deta*c24(k))*e24(2*k-1:2*k)'*Tao*e24(2*k-1:2*k);
        f242=(GTstep_S24*kesi24last(2*(k-1)-1:2*(k-1))-sign(huaA(2,4))*GTstep_S42*kesi42last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S24*kesi24last(2*(k-1)-1:2*(k-1))-sign(huaA(2,4))*GTstep_S42*kesi42last(2*(k-1)-1:2*(k-1)));
        f243=mu*exp(-v*(k-1)*T);
        f24(k)=f241-1/4*f242-f243;
        if f24(k)>0
                 kesi24last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT24=k;
                 triIS24(k)=1;
                 tricount24=tricount24+1;
                
            else
               kesi24last(2*k-1:2*k)=kesi24last(2*(k-1)-1:2*(k-1));
               triIS24(k)=-1; 
        end
        
        e25(2*k-1:2*k)=GTstep_S25*kesi25last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f251=(1+2*deta*c25(k))*e25(2*k-1:2*k)'*Tao*e25(2*k-1:2*k);
        f252=(GTstep_S25*kesi25last(2*(k-1)-1:2*(k-1))-sign(huaA(2,5))*GTstep_S52*kesi52last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S25*kesi25last(2*(k-1)-1:2*(k-1))-sign(huaA(2,5))*GTstep_S52*kesi52last(2*(k-1)-1:2*(k-1)));
        f253=mu*exp(-v*(k-1)*T);
        f25(k)=f251-1/4*f252-f253;
        if f25(k)>0
                 kesi25last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT25=k;
                 triIS25(k)=1;
                 tricount25=tricount25+1;
                
            else
               kesi25last(2*k-1:2*k)=kesi25last(2*(k-1)-1:2*(k-1));
               triIS25(k)=-1; 
        end
        
        e26(2*k-1:2*k)=GTstep_S26*kesi26last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f261=(1+2*deta*c26(k))*e26(2*k-1:2*k)'*Tao*e26(2*k-1:2*k);
        f262=(GTstep_S26*kesi26last(2*(k-1)-1:2*(k-1))-sign(huaA(2,6))*GTstep_S62*kesi62last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S26*kesi26last(2*(k-1)-1:2*(k-1))-sign(huaA(2,6))*GTstep_S62*kesi62last(2*(k-1)-1:2*(k-1)));
        f263=mu*exp(-v*(k-1)*T);
        f26(k)=f261-1/4*f262-f263;
        if f26(k)>0
                 kesi26last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT26=k;
                 triIS26(k)=1;
                 tricount26=tricount26+1;
                
            else
               kesi26last(2*k-1:2*k)=kesi26last(2*(k-1)-1:2*(k-1));
               triIS26(k)=-1; 
        end
        
        e31(2*k-1:2*k)=GTstep_S31*kesi31last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f311=(1+2*deta*c31(k))*e31(2*k-1:2*k)'*Tao*e31(2*k-1:2*k);
        f312=(GTstep_S31*kesi31last(2*(k-1)-1:2*(k-1))-sign(huaA(3,1))*GTstep_S13*kesi13last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S31*kesi31last(2*(k-1)-1:2*(k-1))-sign(huaA(3,1))*GTstep_S13*kesi13last(2*(k-1)-1:2*(k-1)));
        f313=mu*exp(-v*(k-1)*T);
        f31(k)=f311-1/4*f312-f313;
        if f31(k)>0
                 kesi31last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT31=k;
                 triIS31(k)=1;
                 tricount31=tricount31+1;
                
            else
               kesi31last(2*k-1:2*k)=kesi31last(2*(k-1)-1:2*(k-1));
               triIS31(k)=-1; 
        end
        
        e32(2*k-1:2*k)=GTstep_S32*kesi32last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f321=(1+2*deta*c32(k))*e32(2*k-1:2*k)'*Tao*e32(2*k-1:2*k);
        f322=(GTstep_S32*kesi32last(2*(k-1)-1:2*(k-1))-sign(huaA(3,2))*GTstep_S23*kesi23last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S32*kesi32last(2*(k-1)-1:2*(k-1))-sign(huaA(3,2))*GTstep_S23*kesi23last(2*(k-1)-1:2*(k-1)));
        f323=mu*exp(-v*(k-1)*T);
        f32(k)=f321-1/4*f322-f323;
        if f32(k)>0
                 kesi32last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT32=k;
                 triIS32(k)=1;
                 tricount32=tricount32+1;
                
            else
               kesi32last(2*k-1:2*k)=kesi32last(2*(k-1)-1:2*(k-1));
               triIS32(k)=-1; 
        end
        
        e33(2*k-1:2*k)=GTstep_S33*kesi33last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f331=(1+2*deta*c33(k))*e33(2*k-1:2*k)'*Tao*e33(2*k-1:2*k);
        f332=(GTstep_S33*kesi33last(2*(k-1)-1:2*(k-1))-sign(huaA(3,3))*GTstep_S33*kesi33last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S33*kesi33last(2*(k-1)-1:2*(k-1))-sign(huaA(3,3))*GTstep_S33*kesi33last(2*(k-1)-1:2*(k-1)));
        f333=mu*exp(-v*(k-1)*T);
        f33(k)=f331-1/4*f332-f333;
        if f33(k)>0
                 kesi33last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT33=k;
                 triIS33(k)=1;
                 tricount33=tricount33+1;
                
            else
               kesi33last(2*k-1:2*k)=kesi33last(2*(k-1)-1:2*(k-1));
               triIS33(k)=-1; 
        end
        
        e34(2*k-1:2*k)=GTstep_S34*kesi34last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f341=(1+2*deta*c34(k))*e34(2*k-1:2*k)'*Tao*e34(2*k-1:2*k);
        f342=(GTstep_S34*kesi34last(2*(k-1)-1:2*(k-1))-sign(huaA(3,4))*GTstep_S43*kesi43last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S34*kesi34last(2*(k-1)-1:2*(k-1))-sign(huaA(3,4))*GTstep_S43*kesi43last(2*(k-1)-1:2*(k-1)));
        f343=mu*exp(-v*(k-1)*T);
        f34(k)=f341-1/4*f342-f343;
        if f34(k)>0
                 kesi34last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT34=k;
                 triIS34(k)=1;
                 tricount34=tricount34+1;
                
            else
               kesi34last(2*k-1:2*k)=kesi34last(2*(k-1)-1:2*(k-1));
               triIS34(k)=-1; 
        end
        
        e35(2*k-1:2*k)=GTstep_S35*kesi35last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f351=(1+2*deta*c35(k))*e35(2*k-1:2*k)'*Tao*e35(2*k-1:2*k);
        f352=(GTstep_S35*kesi35last(2*(k-1)-1:2*(k-1))-sign(huaA(3,5))*GTstep_S53*kesi53last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S35*kesi35last(2*(k-1)-1:2*(k-1))-sign(huaA(3,5))*GTstep_S53*kesi53last(2*(k-1)-1:2*(k-1)));
        f353=mu*exp(-v*(k-1)*T);
        f35(k)=f351-1/4*f352-f353;
        if f35(k)>0
                 kesi35last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT35=k;
                 triIS35(k)=1;
                 tricount35=tricount35+1;
                
            else
               kesi35last(2*k-1:2*k)=kesi35last(2*(k-1)-1:2*(k-1));
               triIS35(k)=-1; 
        end
        
        e36(2*k-1:2*k)=GTstep_S36*kesi36last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f361=(1+2*deta*c36(k))*e36(2*k-1:2*k)'*Tao*e36(2*k-1:2*k);
        f362=(GTstep_S36*kesi36last(2*(k-1)-1:2*(k-1))-sign(huaA(3,6))*GTstep_S63*kesi63last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S36*kesi36last(2*(k-1)-1:2*(k-1))-sign(huaA(3,6))*GTstep_S63*kesi63last(2*(k-1)-1:2*(k-1)));
        f363=mu*exp(-v*(k-1)*T);
        f36(k)=f361-1/4*f362-f363;
        if f36(k)>0
                 kesi36last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT36=k;
                 triIS36(k)=1;
                 tricount36=tricount36+1;
                
            else
               kesi36last(2*k-1:2*k)=kesi36last(2*(k-1)-1:2*(k-1));
               triIS36(k)=-1; 
        end
        
        e41(2*k-1:2*k)=GTstep_S41*kesi41last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f411=(1+2*deta*c41(k))*e41(2*k-1:2*k)'*Tao*e41(2*k-1:2*k);
        f412=(GTstep_S41*kesi41last(2*(k-1)-1:2*(k-1))-sign(huaA(4,1))*GTstep_S14*kesi14last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S41*kesi41last(2*(k-1)-1:2*(k-1))-sign(huaA(4,1))*GTstep_S14*kesi14last(2*(k-1)-1:2*(k-1)));
        f413=mu*exp(-v*(k-1)*T);
        f41(k)=f411-1/4*f412-f413;
        if f41(k)>0
                 kesi41last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT41=k;
                 triIS41(k)=1;
                 tricount41=tricount41+1;
                
            else
               kesi41last(2*k-1:2*k)=kesi41last(2*(k-1)-1:2*(k-1));
               triIS41(k)=-1; 
        end
        
        e42(2*k-1:2*k)=GTstep_S42*kesi42last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f421=(1+2*deta*c42(k))*e42(2*k-1:2*k)'*Tao*e42(2*k-1:2*k);
        f422=(GTstep_S42*kesi42last(2*(k-1)-1:2*(k-1))-sign(huaA(4,2))*GTstep_S24*kesi24last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S42*kesi42last(2*(k-1)-1:2*(k-1))-sign(huaA(4,2))*GTstep_S24*kesi24last(2*(k-1)-1:2*(k-1)));
        f423=mu*exp(-v*(k-1)*T);
        f42(k)=f421-1/4*f422-f423;
        if f42(k)>0
                 kesi42last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT42=k;
                 triIS42(k)=1;
                 tricount42=tricount42+1;
                
            else
               kesi42last(2*k-1:2*k)=kesi42last(2*(k-1)-1:2*(k-1));
               triIS42(k)=-1; 
        end
        
        e43(2*k-1:2*k)=GTstep_S43*kesi43last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f431=(1+2*deta*c43(k))*e43(2*k-1:2*k)'*Tao*e43(2*k-1:2*k);
        f432=(GTstep_S43*kesi43last(2*(k-1)-1:2*(k-1))-sign(huaA(4,3))*GTstep_S34*kesi34last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S43*kesi43last(2*(k-1)-1:2*(k-1))-sign(huaA(4,3))*GTstep_S34*kesi34last(2*(k-1)-1:2*(k-1)));
        f433=mu*exp(-v*(k-1)*T);
        f43(k)=f431-1/4*f432-f433;
        if f43(k)>0
                 kesi43last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT43=k;
                 triIS43(k)=1;
                 tricount43=tricount43+1;
                
            else
               kesi43last(2*k-1:2*k)=kesi43last(2*(k-1)-1:2*(k-1));
               triIS43(k)=-1; 
        end
        
        e44(2*k-1:2*k)=GTstep_S44*kesi44last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f441=(1+2*deta*c44(k))*e44(2*k-1:2*k)'*Tao*e44(2*k-1:2*k);
        f442=(GTstep_S44*kesi44last(2*(k-1)-1:2*(k-1))-sign(huaA(4,4))*GTstep_S44*kesi44last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S44*kesi44last(2*(k-1)-1:2*(k-1))-sign(huaA(4,4))*GTstep_S44*kesi44last(2*(k-1)-1:2*(k-1)));
        f443=mu*exp(-v*(k-1)*T);
        f44(k)=f441-1/4*f442-f443;
        if f44(k)>0
                 kesi44last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT44=k;
                 triIS44(k)=1;
                 tricount44=tricount44+1;
                
            else
               kesi44last(2*k-1:2*k)=kesi44last(2*(k-1)-1:2*(k-1));
               triIS44(k)=-1; 
        end
        
        e45(2*k-1:2*k)=GTstep_S45*kesi45last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f451=(1+2*deta*c45(k))*e45(2*k-1:2*k)'*Tao*e45(2*k-1:2*k);
        f452=(GTstep_S45*kesi45last(2*(k-1)-1:2*(k-1))-sign(huaA(4,5))*GTstep_S54*kesi54last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S45*kesi45last(2*(k-1)-1:2*(k-1))-sign(huaA(4,5))*GTstep_S54*kesi54last(2*(k-1)-1:2*(k-1)));
        f453=mu*exp(-v*(k-1)*T);
        f45(k)=f451-1/4*f452-f453;
        if f45(k)>0
                 kesi45last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT45=k;
                 triIS45(k)=1;
                 tricount45=tricount45+1;
                
            else
               kesi45last(2*k-1:2*k)=kesi45last(2*(k-1)-1:2*(k-1));
               triIS45(k)=-1; 
        end
        
        e46(2*k-1:2*k)=GTstep_S46*kesi46last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f461=(1+2*deta*c46(k))*e46(2*k-1:2*k)'*Tao*e46(2*k-1:2*k);
        f462=(GTstep_S46*kesi46last(2*(k-1)-1:2*(k-1))-sign(huaA(4,6))*GTstep_S64*kesi64last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S46*kesi46last(2*(k-1)-1:2*(k-1))-sign(huaA(4,6))*GTstep_S64*kesi64last(2*(k-1)-1:2*(k-1)));
        f463=mu*exp(-v*(k-1)*T);
        f46(k)=f461-1/4*f462-f463;
        if f46(k)>0
                 kesi46last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT46=k;
                 triIS46(k)=1;
                 tricount46=tricount46+1;
                
            else
               kesi46last(2*k-1:2*k)=kesi46last(2*(k-1)-1:2*(k-1));
               triIS46(k)=-1; 
        end
        
        e51(2*k-1:2*k)=GTstep_S51*kesi51last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);
        f511=(1+2*deta*c51(k))*e51(2*k-1:2*k)'*Tao*e51(2*k-1:2*k);
        f512=(GTstep_S51*kesi51last(2*(k-1)-1:2*(k-1))-sign(huaA(5,1))*GTstep_S15*kesi15last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S51*kesi51last(2*(k-1)-1:2*(k-1))-sign(huaA(5,1))*GTstep_S15*kesi15last(2*(k-1)-1:2*(k-1)));
        f513=mu*exp(-v*(k-1)*T);
        f51(k)=f511-1/4*f512-f513;
        if f51(k)>0
                 kesi51last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT51=k;
                 triIS51(k)=1;
                 tricount51=tricount51+1;
                
            else
               kesi51last(2*k-1:2*k)=kesi51last(2*(k-1)-1:2*(k-1));
               triIS51(k)=-1; 
        end
        
        e52(2*k-1:2*k)=GTstep_S52*kesi52last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);
        f521=(1+2*deta*c52(k))*e52(2*k-1:2*k)'*Tao*e52(2*k-1:2*k);
        f522=(GTstep_S52*kesi52last(2*(k-1)-1:2*(k-1))-sign(huaA(5,2))*GTstep_S25*kesi25last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S52*kesi52last(2*(k-1)-1:2*(k-1))-sign(huaA(5,2))*GTstep_S25*kesi25last(2*(k-1)-1:2*(k-1)));
        f523=mu*exp(-v*(k-1)*T);
        f52(k)=f521-1/4*f522-f523;
        if f52(k)>0
                 kesi52last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT52=k;
                 triIS52(k)=1;
                 tricount52=tricount52+1;
                
            else
               kesi52last(2*k-1:2*k)=kesi52last(2*(k-1)-1:2*(k-1));
               triIS52(k)=-1; 
        end
        
        e53(2*k-1:2*k)=GTstep_S53*kesi53last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);
        f531=(1+2*deta*c53(k))*e53(2*k-1:2*k)'*Tao*e53(2*k-1:2*k);
        f532=(GTstep_S53*kesi53last(2*(k-1)-1:2*(k-1))-sign(huaA(5,3))*GTstep_S35*kesi35last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S53*kesi53last(2*(k-1)-1:2*(k-1))-sign(huaA(5,3))*GTstep_S35*kesi35last(2*(k-1)-1:2*(k-1)));
        f533=mu*exp(-v*(k-1)*T);
        f53(k)=f531-1/4*f532-f533;
        if f53(k)>0
                 kesi53last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT53=k;
                 triIS53(k)=1;
                 tricount53=tricount53+1;
                
            else
               kesi53last(2*k-1:2*k)=kesi53last(2*(k-1)-1:2*(k-1));
               triIS53(k)=-1; 
        end
        
        e54(2*k-1:2*k)=GTstep_S54*kesi54last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);
        f541=(1+2*deta*c54(k))*e54(2*k-1:2*k)'*Tao*e54(2*k-1:2*k);
        f542=(GTstep_S54*kesi54last(2*(k-1)-1:2*(k-1))-sign(huaA(5,4))*GTstep_S45*kesi45last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S54*kesi54last(2*(k-1)-1:2*(k-1))-sign(huaA(5,4))*GTstep_S45*kesi45last(2*(k-1)-1:2*(k-1)));
        f543=mu*exp(-v*(k-1)*T);
        f54(k)=f541-1/4*f542-f543;
        if f54(k)>0
                 kesi54last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT54=k;
                 triIS54(k)=1;
                 tricount54=tricount54+1;
                
            else
               kesi54last(2*k-1:2*k)=kesi54last(2*(k-1)-1:2*(k-1));
               triIS54(k)=-1; 
        end
        
        e55(2*k-1:2*k)=GTstep_S55*kesi55last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);
        f551=(1+2*deta*c55(k))*e55(2*k-1:2*k)'*Tao*e55(2*k-1:2*k);
        f552=(GTstep_S55*kesi55last(2*(k-1)-1:2*(k-1))-sign(huaA(5,5))*GTstep_S55*kesi55last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S55*kesi55last(2*(k-1)-1:2*(k-1))-sign(huaA(5,5))*GTstep_S55*kesi55last(2*(k-1)-1:2*(k-1)));
        f553=mu*exp(-v*(k-1)*T);
        f55(k)=f551-1/4*f552-f553;
        if f55(k)>0
                 kesi55last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT55=k;
                 triIS55(k)=1;
                 tricount55=tricount55+1;
                
            else
               kesi55last(2*k-1:2*k)=kesi55last(2*(k-1)-1:2*(k-1));
               triIS55(k)=-1; 
        end
        
        e56(2*k-1:2*k)=GTstep_S56*kesi56last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);
        f561=(1+2*deta*c56(k))*e56(2*k-1:2*k)'*Tao*e56(2*k-1:2*k);
        f562=(GTstep_S56*kesi56last(2*(k-1)-1:2*(k-1))-sign(huaA(5,6))*GTstep_S65*kesi65last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S56*kesi56last(2*(k-1)-1:2*(k-1))-sign(huaA(5,6))*GTstep_S65*kesi65last(2*(k-1)-1:2*(k-1)));
        f563=mu*exp(-v*(k-1)*T);
        f56(k)=f561-1/4*f562-f563;
        if f56(k)>0
                 kesi56last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT56=k;
                 triIS56(k)=1;
                 tricount56=tricount56+1;
                
            else
               kesi56last(2*k-1:2*k)=kesi56last(2*(k-1)-1:2*(k-1));
               triIS56(k)=-1; 
        end
        
        e61(2*k-1:2*k)=GTstep_S61*kesi61last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k);
        f611=(1+2*deta*c61(k))*e61(2*k-1:2*k)'*Tao*e61(2*k-1:2*k);
        f612=(GTstep_S61*kesi61last(2*(k-1)-1:2*(k-1))-sign(huaA(6,1))*GTstep_S16*kesi16last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S61*kesi61last(2*(k-1)-1:2*(k-1))-sign(huaA(6,1))*GTstep_S16*kesi16last(2*(k-1)-1:2*(k-1)));
        f613=mu*exp(-v*(k-1)*T);
        f61(k)=f611-1/4*f612-f613;
        if f61(k)>0
                 kesi61last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT61=k;
                 triIS61(k)=1;
                 tricount61=tricount61+1;
                
            else
               kesi61last(2*k-1:2*k)=kesi61last(2*(k-1)-1:2*(k-1));
               triIS61(k)=-1; 
        end
        
        e62(2*k-1:2*k)=GTstep_S62*kesi62last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k);
        f621=(1+2*deta*c62(k))*e62(2*k-1:2*k)'*Tao*e62(2*k-1:2*k);
        f622=(GTstep_S62*kesi62last(2*(k-1)-1:2*(k-1))-sign(huaA(6,2))*GTstep_S26*kesi26last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S62*kesi62last(2*(k-1)-1:2*(k-1))-sign(huaA(6,2))*GTstep_S26*kesi26last(2*(k-1)-1:2*(k-1)));
        f623=mu*exp(-v*(k-1)*T);
        f62(k)=f621-1/4*f622-f623;
        if f62(k)>0
                 kesi62last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT62=k;
                 triIS62(k)=1;
                 tricount62=tricount62+1;
                
            else
               kesi62last(2*k-1:2*k)=kesi62last(2*(k-1)-1:2*(k-1));
               triIS62(k)=-1; 
        end
        
        e63(2*k-1:2*k)=GTstep_S63*kesi63last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k);
        f631=(1+2*deta*c63(k))*e63(2*k-1:2*k)'*Tao*e63(2*k-1:2*k);
        f632=(GTstep_S63*kesi63last(2*(k-1)-1:2*(k-1))-sign(huaA(6,3))*GTstep_S36*kesi36last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S63*kesi63last(2*(k-1)-1:2*(k-1))-sign(huaA(6,3))*GTstep_S36*kesi36last(2*(k-1)-1:2*(k-1)));
        f633=mu*exp(-v*(k-1)*T);
        f63(k)=f631-1/4*f632-f633;
        if f63(k)>0
                 kesi63last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT63=k;
                 triIS63(k)=1;
                 tricount63=tricount63+1;
                
            else
               kesi63last(2*k-1:2*k)=kesi63last(2*(k-1)-1:2*(k-1));
               triIS63(k)=-1; 
        end
        
        e64(2*k-1:2*k)=GTstep_S64*kesi64last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k);
        f641=(1+2*deta*c64(k))*e64(2*k-1:2*k)'*Tao*e64(2*k-1:2*k);
        f642=(GTstep_S64*kesi64last(2*(k-1)-1:2*(k-1))-sign(huaA(6,4))*GTstep_S46*kesi46last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S64*kesi64last(2*(k-1)-1:2*(k-1))-sign(huaA(6,4))*GTstep_S46*kesi46last(2*(k-1)-1:2*(k-1)));
        f643=mu*exp(-v*(k-1)*T);
        f64(k)=f641-1/4*f642-f643;
        if f64(k)>0
                 kesi64last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT64=k;
                 triIS64(k)=1;
                 tricount64=tricount64+1;
                
            else
               kesi64last(2*k-1:2*k)=kesi64last(2*(k-1)-1:2*(k-1));
               triIS64(k)=-1; 
        end
        
        e65(2*k-1:2*k)=GTstep_S65*kesi65last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k);
        f651=(1+2*deta*c65(k))*e65(2*k-1:2*k)'*Tao*e65(2*k-1:2*k);
        f652=(GTstep_S65*kesi65last(2*(k-1)-1:2*(k-1))-sign(huaA(6,5))*GTstep_S56*kesi56last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S65*kesi65last(2*(k-1)-1:2*(k-1))-sign(huaA(6,5))*GTstep_S56*kesi56last(2*(k-1)-1:2*(k-1)));
        f653=mu*exp(-v*(k-1)*T);
        f65(k)=f651-1/4*f652-f653;
        if f65(k)>0
                 kesi65last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT65=k;
                 triIS65(k)=1;
                 tricount65=tricount65+1;
                
            else
               kesi65last(2*k-1:2*k)=kesi65last(2*(k-1)-1:2*(k-1));
               triIS65(k)=-1; 
        end
        
        e66(2*k-1:2*k)=GTstep_S66*kesi66last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k); 
        f661=(1+2*deta*c66(k))*e66(2*k-1:2*k)'*Tao*e66(2*k-1:2*k);
        f662=(GTstep_S66*kesi66last(2*(k-1)-1:2*(k-1))-sign(huaA(6,6))*GTstep_S66*kesi66last(2*(k-1)-1:2*(k-1)))'*Tao*(GTstep_S66*kesi66last(2*(k-1)-1:2*(k-1))-sign(huaA(6,6))*GTstep_S66*kesi66last(2*(k-1)-1:2*(k-1)));
        f663=mu*exp(-v*(k-1)*T);
        f66(k)=f661-1/4*f662-f663;
        if f66(k)>0
                 kesi66last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT66=k;
                 triIS66(k)=1;
                 tricount66=tricount66+1;
                
            else
               kesi66last(2*k-1:2*k)=kesi66last(2*(k-1)-1:2*(k-1));
               triIS66(k)=-1; 
        end
        
       %% judge whether the ETM-b for edge (i,0) is triggered
        e10(2*k-1:2*k)=GTstep_S10*kesi10last(2*(k-1)-1:2*(k-1))-kesi1(2*k-1:2*k);
        f101=(1+2*deta*c1(k))*e10(2*k-1:2*k)'*Tao*e10(2*k-1:2*k);
        f102=(GTstep_S10*kesi10last(2*(k-1)-1:2*(k-1))-D(1,1)*x0(2*k-1:2*k))'*Tao*(GTstep_S10*kesi10last(2*(k-1)-1:2*(k-1))-D(1,1)*x0(2*k-1:2*k));
        f103=mu*exp(-v*(k-1)*T);
        f10(k)=f101-f102-f103; %% f_{i0}(t) in (24)
        if f10(k)>0
                 kesi10last(2*k-1:2*k)=kesi1(2*k-1:2*k); %% if tigger, \vartheta_{i}(t_{\sigma}^{i0}) is updated to \vartheta_{i}(t)
                 triT10=k; %% triggering instant t_{\sigma}^{i0} of ETM-b for edge (i,0) 
                 triIS10(k)=1;  %% trigger flag of ETM-b for edge (i,0), if trigger triISi0=1, else triISi0=-1
                 tricount10=tricount10+1; %% trigger number of ETM-b for edge (i,0) 
                
            else
               kesi10last(2*k-1:2*k)=kesi10last(2*(k-1)-1:2*(k-1)); %% if not tigger, \vartheta_{i}(t_{\sigma}^{i0}) remains unchanged
               triIS10(k)=-1; %% trigger number of ETM-b for edge (i,0) 
        end

        e20(2*k-1:2*k)=GTstep_S20*kesi20last(2*(k-1)-1:2*(k-1))-kesi2(2*k-1:2*k);
        f201=(1+2*deta*c2(k))*e20(2*k-1:2*k)'*Tao*e20(2*k-1:2*k);
        f202=(GTstep_S20*kesi20last(2*(k-1)-1:2*(k-1))-D(2,2)*x0(2*k-1:2*k))'*Tao*(GTstep_S20*kesi20last(2*(k-1)-1:2*(k-1))-D(2,2)*x0(2*k-1:2*k));
        f203=mu*exp(-v*(k-1)*T);
        f20(k)=f201-f202-f203;
        if f20(k)>0
                 kesi20last(2*k-1:2*k)=kesi2(2*k-1:2*k);
                 triT20=k;
                 triIS20(k)=1;
                 tricount20=tricount20+1;
                
            else
               kesi20last(2*k-1:2*k)=kesi20last(2*(k-1)-1:2*(k-1));
               triIS20(k)=-1; 
        end
        
        e30(2*k-1:2*k)=GTstep_S30*kesi30last(2*(k-1)-1:2*(k-1))-kesi3(2*k-1:2*k);
        f301=(1+2*deta*c3(k))*e30(2*k-1:2*k)'*Tao*e30(2*k-1:2*k);
        f302=(GTstep_S30*kesi30last(2*(k-1)-1:2*(k-1))-D(3,3)*x0(2*k-1:2*k))'*Tao*(GTstep_S30*kesi30last(2*(k-1)-1:2*(k-1))-D(3,3)*x0(2*k-1:2*k));
        f303=mu*exp(-v*(k-1)*T);
        f30(k)=f301-f302-f303;
        if f30(k)>0
                 kesi30last(2*k-1:2*k)=kesi3(2*k-1:2*k);
                 triT30=k;
                 triIS30(k)=1;
                 tricount30=tricount30+1;
                
            else
               kesi30last(2*k-1:2*k)=kesi30last(2*(k-1)-1:2*(k-1));
               triIS30(k)=-1; 
        end
        e40(2*k-1:2*k)=GTstep_S40*kesi40last(2*(k-1)-1:2*(k-1))-kesi4(2*k-1:2*k);
        f401=(1+2*deta*c4(k))*e40(2*k-1:2*k)'*Tao*e40(2*k-1:2*k);
        f402=(GTstep_S40*kesi40last(2*(k-1)-1:2*(k-1))-D(4,4)*x0(2*k-1:2*k))'*Tao*(GTstep_S40*kesi40last(2*(k-1)-1:2*(k-1))-D(4,4)*x0(2*k-1:2*k));
        f403=mu*exp(-v*(k-1)*T);
        f40(k)=f401-f402-f403;
        if f40(k)>0
                 kesi40last(2*k-1:2*k)=kesi4(2*k-1:2*k);
                 triT40=k;
                 triIS40(k)=1;
                 tricount40=tricount40+1;
                
            else
               kesi40last(2*k-1:2*k)=kesi40last(2*(k-1)-1:2*(k-1));
               triIS40(k)=-1; 
        end
        e50(2*k-1:2*k)=GTstep_S50*kesi50last(2*(k-1)-1:2*(k-1))-kesi5(2*k-1:2*k);  
        f501=(1+2*deta*c5(k))*e50(2*k-1:2*k)'*Tao*e50(2*k-1:2*k);
        f502=(GTstep_S50*kesi50last(2*(k-1)-1:2*(k-1))-D(5,5)*x0(2*k-1:2*k))'*Tao*(GTstep_S50*kesi50last(2*(k-1)-1:2*(k-1))-D(5,5)*x0(2*k-1:2*k));
        f503=mu*exp(-v*(k-1)*T);
        f50(k)=f501-f502-f503;
        if f50(k)>0
                 kesi50last(2*k-1:2*k)=kesi5(2*k-1:2*k);
                 triT50=k;
                 triIS50(k)=1;
                 tricount50=tricount50+1;
                
            else
               kesi50last(2*k-1:2*k)=kesi50last(2*(k-1)-1:2*(k-1));
               triIS50(k)=-1; 
        end
        e60(2*k-1:2*k)=GTstep_S60*kesi60last(2*(k-1)-1:2*(k-1))-kesi6(2*k-1:2*k); 
        f601=(1+2*deta*c6(k))*e60(2*k-1:2*k)'*Tao*e60(2*k-1:2*k);
        f602=(GTstep_S60*kesi60last(2*(k-1)-1:2*(k-1))-D(6,6)*x0(2*k-1:2*k))'*Tao*(GTstep_S60*kesi60last(2*(k-1)-1:2*(k-1))-D(6,6)*x0(2*k-1:2*k));
        f603=mu*exp(-v*(k-1)*T);
        f60(k)=f601-f602-f603;
        if f60(k)>0
                 kesi60last(2*k-1:2*k)=kesi6(2*k-1:2*k);
                 triT60=k;
                 triIS60(k)=1;
                 tricount60=tricount60+1;
                
            else
               kesi60last(2*k-1:2*k)=kesi60last(2*(k-1)-1:2*(k-1));
               triIS60(k)=-1; 
        end   
            
       %% y_{i}(t) in (1) 
        y1(k)=C1*x1(2*k-1:2*k); 
        y2(k)=C2*x2(2*k-1:2*k);
        y3(k)=C3*x3(3*k-2:3*k);
        y4(k)=C4*x4(3*k-2:3*k);
        y5(k)=C5*x5(4*k-3:4*k);
        y6(k)=C6*x6(4*k-3:4*k);
        y0(k)=C0*x0(2*k-1:2*k);
         
       %% updata e^{S(t-t_{\sigma}^{ij})}
        GTstep_S11=expm(A0*(k-triT11)*T);
        GTstep_S12=expm(A0*(k-triT12)*T);
        GTstep_S13=expm(A0*(k-triT13)*T);
        GTstep_S14=expm(A0*(k-triT14)*T);
        GTstep_S15=expm(A0*(k-triT15)*T);
        GTstep_S16=expm(A0*(k-triT16)*T);
        
        GTstep_S21=expm(A0*(k-triT21)*T);
        GTstep_S22=expm(A0*(k-triT22)*T);
        GTstep_S23=expm(A0*(k-triT23)*T);
        GTstep_S24=expm(A0*(k-triT24)*T);
        GTstep_S25=expm(A0*(k-triT25)*T);
        GTstep_S26=expm(A0*(k-triT26)*T);
        
        GTstep_S31=expm(A0*(k-triT31)*T);
        GTstep_S32=expm(A0*(k-triT32)*T);
        GTstep_S33=expm(A0*(k-triT33)*T);
        GTstep_S34=expm(A0*(k-triT34)*T);
        GTstep_S35=expm(A0*(k-triT35)*T);
        GTstep_S36=expm(A0*(k-triT36)*T);
        
        GTstep_S41=expm(A0*(k-triT41)*T);
        GTstep_S42=expm(A0*(k-triT42)*T);
        GTstep_S43=expm(A0*(k-triT43)*T);
        GTstep_S44=expm(A0*(k-triT44)*T);
        GTstep_S45=expm(A0*(k-triT45)*T);
        GTstep_S46=expm(A0*(k-triT46)*T);
        
        GTstep_S51=expm(A0*(k-triT51)*T);
        GTstep_S52=expm(A0*(k-triT52)*T);
        GTstep_S53=expm(A0*(k-triT53)*T);
        GTstep_S54=expm(A0*(k-triT54)*T);
        GTstep_S55=expm(A0*(k-triT55)*T);
        GTstep_S56=expm(A0*(k-triT56)*T);
        
        GTstep_S61=expm(A0*(k-triT61)*T);
        GTstep_S62=expm(A0*(k-triT62)*T);
        GTstep_S63=expm(A0*(k-triT63)*T);
        GTstep_S64=expm(A0*(k-triT64)*T);
        GTstep_S65=expm(A0*(k-triT65)*T);
        GTstep_S66=expm(A0*(k-triT66)*T);
        
       %% updata e^{S(t-t_{\sigma}^{i0})}
        GTstep_S10=expm(A0*(k-triT10)*T);
        GTstep_S20=expm(A0*(k-triT20)*T);
        GTstep_S30=expm(A0*(k-triT30)*T);
        GTstep_S40=expm(A0*(k-triT40)*T);
        GTstep_S50=expm(A0*(k-triT50)*T);
        GTstep_S60=expm(A0*(k-triT60)*T);
        
       %% u_{i}(t) in (27a)
        u1(k)=K11*eita1(2*k-1:2*k)+K21*kesi1(2*k-1:2*k);
        u2(k)=K12*eita2(2*k-1:2*k)+K22*kesi2(2*k-1:2*k);
        u3(k)=K13*eita3(3*k-2:3*k)+K23*kesi3(2*k-1:2*k);
        u4(k)=K14*eita4(3*k-2:3*k)+K24*kesi4(2*k-1:2*k);
        u5(k)=K15*eita5(4*k-3:4*k)+K25*kesi5(2*k-1:2*k);
        u6(k)=K16*eita6(4*k-3:4*k)+K26*kesi6(2*k-1:2*k);
        
       %% The iteration value of the variables at the next time
       %% x_{i}(t) in (1) and (2)   
        x1(2*(k+1)-1:2*(k+1))=x1(2*k-1:2*k)+T*(A1*x1(2*k-1:2*k)+B1*u1(k));
        x2(2*(k+1)-1:2*(k+1))=x2(2*k-1:2*k)+T*(A2*x2(2*k-1:2*k)+B2*u2(k));
        x3(3*(k+1)-2:3*(k+1))=x3(3*k-2:3*k)+T*(A3*x3(3*k-2:3*k)+B3*u3(k));
        x4(3*(k+1)-2:3*(k+1))=x4(3*k-2:3*k)+T*(A4*x4(3*k-2:3*k)+B4*u4(k));
        x5(4*(k+1)-3:4*(k+1))=x5(4*k-3:4*k)+T*(A5*x5(4*k-3:4*k)+B5*u5(k));
        x6(4*(k+1)-3:4*(k+1))=x6(4*k-3:4*k)+T*(A6*x6(4*k-3:4*k)+B6*u6(k));
        x0(2*(k+1)-1:2*(k+1))=x0(2*k-1:2*k)+T*(A0*x0(2*k-1:2*k));            
       
       %% compensator \vartheta_{i}(t) in (19a)  
        kesi1(2*(k+1)-1:2*(k+1))=kesi1(2*k-1:2*k)+T*(A0*kesi1(2*k-1:2*k))+T*K*(c11(k)*abs(huaA(1,1))*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k))+c12(k)*abs(huaA(1,2))*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k))+c13(k)*abs(huaA(1,3))*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k))+c14(k)*abs(huaA(1,4))*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k))+c15(k)*abs(huaA(1,5))*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k))+c16(k)*abs(huaA(1,6))*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k))+c1(k)*huaB(1,1)*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k)));
        kesi2(2*(k+1)-1:2*(k+1))=kesi2(2*k-1:2*k)+T*(A0*kesi2(2*k-1:2*k))+T*K*(c21(k)*abs(huaA(2,1))*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k))+c22(k)*abs(huaA(2,2))*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k))+c23(k)*abs(huaA(2,3))*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k))+c24(k)*abs(huaA(2,4))*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k))+c25(k)*abs(huaA(2,5))*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k))+c26(k)*abs(huaA(2,6))*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k))+c2(k)*huaB(2,2)*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k)));
        kesi3(2*(k+1)-1:2*(k+1))=kesi3(2*k-1:2*k)+T*(A0*kesi3(2*k-1:2*k))+T*K*(c31(k)*abs(huaA(3,1))*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k))+c32(k)*abs(huaA(3,2))*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k))+c33(k)*abs(huaA(3,3))*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k))+c34(k)*abs(huaA(3,4))*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k))+c35(k)*abs(huaA(3,5))*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k))+c36(k)*abs(huaA(3,6))*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k))+c3(k)*huaB(3,3)*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k)));
        kesi4(2*(k+1)-1:2*(k+1))=kesi4(2*k-1:2*k)+T*(A0*kesi4(2*k-1:2*k))+T*K*(c41(k)*abs(huaA(4,1))*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k))+c42(k)*abs(huaA(4,2))*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k))+c43(k)*abs(huaA(4,3))*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k))+c44(k)*abs(huaA(4,4))*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k))+c45(k)*abs(huaA(4,5))*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k))+c46(k)*abs(huaA(4,6))*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k))+c4(k)*huaB(4,4)*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k)));
        kesi5(2*(k+1)-1:2*(k+1))=kesi5(2*k-1:2*k)+T*(A0*kesi5(2*k-1:2*k))+T*K*(c51(k)*abs(huaA(5,1))*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k))+c52(k)*abs(huaA(5,2))*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k))+c53(k)*abs(huaA(5,3))*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k))+c54(k)*abs(huaA(5,4))*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k))+c55(k)*abs(huaA(5,5))*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k))+c56(k)*abs(huaA(5,6))*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k))+c5(k)*huaB(5,5)*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k)));
        kesi6(2*(k+1)-1:2*(k+1))=kesi6(2*k-1:2*k)+T*(A0*kesi6(2*k-1:2*k))+T*K*(c61(k)*abs(huaA(6,1))*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k))+c62(k)*abs(huaA(6,2))*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k))+c63(k)*abs(huaA(6,3))*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k))+c64(k)*abs(huaA(6,4))*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k))+c65(k)*abs(huaA(6,5))*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k))+c66(k)*abs(huaA(6,6))*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k))+c6(k)*huaB(6,6)*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k)));          
 
       %% state observer \psi_{i}(t) in (27b)
        eita1(2*(k+1)-1:2*(k+1))=eita1(2*k-1:2*k)+T*(A1*eita1(2*k-1:2*k))+T*(B1*u1(k))+T*F1*(C1*eita1(2*k-1:2*k)-y1(k));
        eita2(2*(k+1)-1:2*(k+1))=eita2(2*k-1:2*k)+T*(A2*eita2(2*k-1:2*k))+T*(B2*u2(k))+T*F2*(C2*eita2(2*k-1:2*k)-y2(k));
        eita3(3*(k+1)-2:3*(k+1))=eita3(3*k-2:3*k)+T*(A3*eita3(3*k-2:3*k))+T*(B3*u3(k))+T*F3*(C3*eita3(3*k-2:3*k)-y3(k));
        eita4(3*(k+1)-2:3*(k+1))=eita4(3*k-2:3*k)+T*(A4*eita4(3*k-2:3*k))+T*(B4*u4(k))+T*F4*(C4*eita4(3*k-2:3*k)-y4(k));
        eita5(4*(k+1)-3:4*(k+1))=eita5(4*k-3:4*k)+T*(A5*eita5(4*k-3:4*k))+T*(B5*u5(k))+T*F5*(C5*eita5(4*k-3:4*k)-y5(k));
        eita6(4*(k+1)-3:4*(k+1))=eita6(4*k-3:4*k)+T*(A6*eita6(4*k-3:4*k))+T*(B6*u6(k))+T*F6*(C6*eita6(4*k-3:4*k)-y6(k));
        
         %% c_{ij}(t) in (19b)
        c11(k+1)= c11(k)+T*l11*(abs(huaA(1,1))*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k))'*Tao*(GTstep_S11*kesi11last(2*k-1:2*k)-sign(huaA(1,1))*GTstep_S11*kesi11last(2*k-1:2*k)));
        c12(k+1)= c12(k)+T*l12*(abs(huaA(1,2))*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k))'*Tao*(GTstep_S12*kesi12last(2*k-1:2*k)-sign(huaA(1,2))*GTstep_S21*kesi21last(2*k-1:2*k)));
        c13(k+1)= c13(k)+T*l13*(abs(huaA(1,3))*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k))'*Tao*(GTstep_S13*kesi13last(2*k-1:2*k)-sign(huaA(1,3))*GTstep_S31*kesi31last(2*k-1:2*k)));
        c14(k+1)= c14(k)+T*l14*(abs(huaA(1,4))*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k))'*Tao*(GTstep_S14*kesi14last(2*k-1:2*k)-sign(huaA(1,4))*GTstep_S41*kesi41last(2*k-1:2*k)));
        c15(k+1)= c15(k)+T*l15*(abs(huaA(1,5))*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k))'*Tao*(GTstep_S15*kesi15last(2*k-1:2*k)-sign(huaA(1,5))*GTstep_S51*kesi51last(2*k-1:2*k)));
        c16(k+1)= c16(k)+T*l16*(abs(huaA(1,6))*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k))'*Tao*(GTstep_S16*kesi16last(2*k-1:2*k)-sign(huaA(1,6))*GTstep_S61*kesi61last(2*k-1:2*k)));
        
        c21(k+1)= c21(k)+T*l21*(abs(huaA(2,1))*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k))'*Tao*(GTstep_S21*kesi21last(2*k-1:2*k)-sign(huaA(2,1))*GTstep_S12*kesi12last(2*k-1:2*k)));
        c22(k+1)= c22(k)+T*l22*(abs(huaA(2,2))*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k))'*Tao*(GTstep_S22*kesi22last(2*k-1:2*k)-sign(huaA(2,2))*GTstep_S22*kesi22last(2*k-1:2*k)));
        c23(k+1)= c23(k)+T*l23*(abs(huaA(2,3))*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k))'*Tao*(GTstep_S23*kesi23last(2*k-1:2*k)-sign(huaA(2,3))*GTstep_S32*kesi32last(2*k-1:2*k)));
        c24(k+1)= c24(k)+T*l24*(abs(huaA(2,4))*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k))'*Tao*(GTstep_S24*kesi24last(2*k-1:2*k)-sign(huaA(2,4))*GTstep_S42*kesi42last(2*k-1:2*k)));
        c25(k+1)= c25(k)+T*l25*(abs(huaA(2,5))*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k))'*Tao*(GTstep_S25*kesi25last(2*k-1:2*k)-sign(huaA(2,5))*GTstep_S52*kesi52last(2*k-1:2*k)));
        c26(k+1)= c26(k)+T*l26*(abs(huaA(2,6))*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k))'*Tao*(GTstep_S26*kesi26last(2*k-1:2*k)-sign(huaA(2,6))*GTstep_S62*kesi62last(2*k-1:2*k)));
       
        c31(k+1)= c31(k)+T*l31*(abs(huaA(3,1))*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k))'*Tao*(GTstep_S31*kesi31last(2*k-1:2*k)-sign(huaA(3,1))*GTstep_S13*kesi13last(2*k-1:2*k)));
        c32(k+1)= c32(k)+T*l32*(abs(huaA(3,2))*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k))'*Tao*(GTstep_S32*kesi32last(2*k-1:2*k)-sign(huaA(3,2))*GTstep_S23*kesi23last(2*k-1:2*k)));
        c33(k+1)= c33(k)+T*l33*(abs(huaA(3,3))*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k))'*Tao*(GTstep_S33*kesi33last(2*k-1:2*k)-sign(huaA(3,3))*GTstep_S33*kesi33last(2*k-1:2*k)));
        c34(k+1)= c34(k)+T*l34*(abs(huaA(3,4))*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k))'*Tao*(GTstep_S34*kesi34last(2*k-1:2*k)-sign(huaA(3,4))*GTstep_S43*kesi43last(2*k-1:2*k)));
        c35(k+1)= c35(k)+T*l35*(abs(huaA(3,5))*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k))'*Tao*(GTstep_S35*kesi35last(2*k-1:2*k)-sign(huaA(3,5))*GTstep_S53*kesi53last(2*k-1:2*k)));
        c36(k+1)= c36(k)+T*l36*(abs(huaA(3,6))*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k))'*Tao*(GTstep_S36*kesi36last(2*k-1:2*k)-sign(huaA(3,6))*GTstep_S63*kesi63last(2*k-1:2*k)));
        
        c41(k+1)= c41(k)+T*l41*(abs(huaA(4,1))*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k))'*Tao*(GTstep_S41*kesi41last(2*k-1:2*k)-sign(huaA(4,1))*GTstep_S14*kesi14last(2*k-1:2*k)));
        c42(k+1)= c42(k)+T*l42*(abs(huaA(4,2))*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k))'*Tao*(GTstep_S42*kesi42last(2*k-1:2*k)-sign(huaA(4,2))*GTstep_S24*kesi24last(2*k-1:2*k)));
        c43(k+1)= c43(k)+T*l43*(abs(huaA(4,3))*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k))'*Tao*(GTstep_S43*kesi43last(2*k-1:2*k)-sign(huaA(4,3))*GTstep_S34*kesi34last(2*k-1:2*k)));
        c44(k+1)= c44(k)+T*l44*(abs(huaA(4,4))*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k))'*Tao*(GTstep_S44*kesi44last(2*k-1:2*k)-sign(huaA(4,4))*GTstep_S44*kesi44last(2*k-1:2*k)));
        c45(k+1)= c45(k)+T*l45*(abs(huaA(4,5))*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k))'*Tao*(GTstep_S45*kesi45last(2*k-1:2*k)-sign(huaA(4,5))*GTstep_S54*kesi54last(2*k-1:2*k)));
        c46(k+1)= c46(k)+T*l46*(abs(huaA(4,6))*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k))'*Tao*(GTstep_S46*kesi46last(2*k-1:2*k)-sign(huaA(4,6))*GTstep_S64*kesi64last(2*k-1:2*k)));
        
        c51(k+1)= c51(k)+T*l51*(abs(huaA(5,1))*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k))'*Tao*(GTstep_S51*kesi51last(2*k-1:2*k)-sign(huaA(5,1))*GTstep_S15*kesi15last(2*k-1:2*k)));
        c52(k+1)= c52(k)+T*l52*(abs(huaA(5,2))*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k))'*Tao*(GTstep_S52*kesi52last(2*k-1:2*k)-sign(huaA(5,2))*GTstep_S25*kesi25last(2*k-1:2*k)));
        c53(k+1)= c53(k)+T*l53*(abs(huaA(5,3))*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k))'*Tao*(GTstep_S53*kesi53last(2*k-1:2*k)-sign(huaA(5,3))*GTstep_S35*kesi35last(2*k-1:2*k)));
        c54(k+1)= c54(k)+T*l54*(abs(huaA(5,4))*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k))'*Tao*(GTstep_S54*kesi54last(2*k-1:2*k)-sign(huaA(5,4))*GTstep_S45*kesi45last(2*k-1:2*k)));
        c55(k+1)= c55(k)+T*l55*(abs(huaA(5,5))*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k))'*Tao*(GTstep_S55*kesi55last(2*k-1:2*k)-sign(huaA(5,5))*GTstep_S55*kesi55last(2*k-1:2*k)));
        c56(k+1)= c56(k)+T*l56*(abs(huaA(5,6))*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k))'*Tao*(GTstep_S56*kesi56last(2*k-1:2*k)-sign(huaA(5,6))*GTstep_S65*kesi65last(2*k-1:2*k)));
        
        c61(k+1)= c61(k)+T*l61*(abs(huaA(6,1))*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k))'*Tao*(GTstep_S61*kesi61last(2*k-1:2*k)-sign(huaA(6,1))*GTstep_S16*kesi16last(2*k-1:2*k)));
        c62(k+1)= c62(k)+T*l62*(abs(huaA(6,2))*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k))'*Tao*(GTstep_S62*kesi62last(2*k-1:2*k)-sign(huaA(6,2))*GTstep_S26*kesi26last(2*k-1:2*k)));
        c63(k+1)= c63(k)+T*l63*(abs(huaA(6,3))*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k))'*Tao*(GTstep_S63*kesi63last(2*k-1:2*k)-sign(huaA(6,3))*GTstep_S36*kesi36last(2*k-1:2*k)));
        c64(k+1)= c64(k)+T*l64*(abs(huaA(6,4))*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k))'*Tao*(GTstep_S64*kesi64last(2*k-1:2*k)-sign(huaA(6,4))*GTstep_S46*kesi46last(2*k-1:2*k)));
        c65(k+1)= c65(k)+T*l65*(abs(huaA(6,5))*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k))'*Tao*(GTstep_S65*kesi65last(2*k-1:2*k)-sign(huaA(6,5))*GTstep_S56*kesi56last(2*k-1:2*k)));
        c66(k+1)= c66(k)+T*l66*(abs(huaA(6,6))*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k))'*Tao*(GTstep_S66*kesi66last(2*k-1:2*k)-sign(huaA(6,6))*GTstep_S66*kesi66last(2*k-1:2*k)));
        
       %% c_{i0}(t) in (19c)
        c1(k+1)= c1(k)+T*l1*(huaB(1,1)*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k))'*Tao*(GTstep_S10*kesi10last(2*k-1:2*k)-D(1,1)*x0(2*k-1:2*k)));
        c2(k+1)= c2(k)+T*l2*(huaB(2,2)*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k))'*Tao*(GTstep_S20*kesi20last(2*k-1:2*k)-D(2,2)*x0(2*k-1:2*k)));
        c3(k+1)= c3(k)+T*l3*(huaB(3,3)*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k))'*Tao*(GTstep_S30*kesi30last(2*k-1:2*k)-D(3,3)*x0(2*k-1:2*k)));
        c4(k+1)= c4(k)+T*l4*(huaB(4,4)*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k))'*Tao*(GTstep_S40*kesi40last(2*k-1:2*k)-D(4,4)*x0(2*k-1:2*k)));
        c5(k+1)= c5(k)+T*l5*(huaB(5,5)*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k))'*Tao*(GTstep_S50*kesi50last(2*k-1:2*k)-D(5,5)*x0(2*k-1:2*k)));
        c6(k+1)= c6(k)+T*l6*(huaB(6,6)*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k))'*Tao*(GTstep_S60*kesi60last(2*k-1:2*k)-D(6,6)*x0(2*k-1:2*k)));
 
        else            
           %% y_{i}(t) in (1) 
            y1(k)=C1*x1(2*k-1:2*k);
            y2(k)=C2*x2(2*k-1:2*k);
            y3(k)=C3*x3(3*k-2:3*k);
            y4(k)=C4*x4(3*k-2:3*k);
            y5(k)=C5*x5(4*k-3:4*k);
            y6(k)=C6*x6(4*k-3:4*k);
            y0(k)=C0*x0(2*k-1:2*k);
        end
    end
end


%% The compensator estimation errors $\varsigma_{i}(t)$ of six follower agents.
figure(1); 
subplot(2,1,1);
plot(0:T:T*(kmax-1),kesi1(1:2:length(kesi1))-D(1,1)*x0(1:2:length(x0)),'k--','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi2(1:2:length(kesi1))-D(2,2)*x0(1:2:length(x0)),'r--','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi3(1:2:length(kesi1))-D(3,3)*x0(1:2:length(x0)),'b--','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi4(1:2:length(kesi1))-D(4,4)*x0(1:2:length(x0)),'-.g','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi5(1:2:length(kesi1))-D(5,5)*x0(1:2:length(x0)),'-.m','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi6(1:2:length(kesi1))-D(6,6)*x0(1:2:length(x0)),'-.c','linewidth',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$\varsigma_{i1}(t)$','Interpreter','latex','FontSize',25); 
h=legend('Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',20);
set(h,'Orientation','horizon')
subplot(2,1,2);
plot(0:T:T*(kmax-1),kesi1(2:2:length(kesi1))-D(1,1)*x0(2:2:length(x0)),'k--','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi2(2:2:length(kesi1))-D(2,2)*x0(2:2:length(x0)),'r--','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi3(2:2:length(kesi1))-D(3,3)*x0(2:2:length(x0)),'b--','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi4(2:2:length(kesi1))-D(4,4)*x0(2:2:length(x0)),'-.g','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi5(2:2:length(kesi1))-D(5,5)*x0(2:2:length(x0)),'-.m','linewidth',3);
hold on;
plot(0:T:T*(kmax-1),kesi6(2:2:length(kesi1))-D(6,6)*x0(2:2:length(x0)),'-.c','linewidth',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$\varsigma_{i2}(t)$','Interpreter','latex','FontSize',25); 

%% The adaptive coupling weights $c_{ij}(t)$
figure(2); 
plot(0:T:T*(kmax-1),c12(1:length(c12)),'k--','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c13(1:length(c12)),'r--','LineWidth',3,'MarkerFaceColor','r','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c14(1:length(c12)),'b--','LineWidth',3,'MarkerFaceColor','b','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c15(1:length(c12)),'g--','LineWidth',3,'MarkerFaceColor','g','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c16(1:length(c12)),'-.m','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c23(1:length(c12)),'-.c','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c45(1:length(c12)),'-.b','LineWidth',3,'MarkerFaceColor','c','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c56(1:length(c12)),'-.g','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$c_{ij}(t)$','Interpreter','latex','FontSize',25); 
h=legend('$c_{12}$','$c_{13}$','$c_{14}$','$c_{15}$','$c_{16}$','$c_{23}$','$c_{45}$','$c_{56}$');
set(h,'Fontsize',25);
set(h,'Orientation','horizon')
set(h,'Interpreter','latex')

%% The adaptive coupling weights $c_{i0}(t)$
figure(3); 
plot(0:T:T*(kmax-1),c1(1:length(c1)),'k--','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),c6(1:length(c6)),'-.r','LineWidth',3,'MarkerFaceColor','r','MarkerSize',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$c_{i0}(t)$','Interpreter','latex','FontSize',25); 
h=legend('$c_{10}$','$c_{60}$');
set(h,'Fontsize',25);
set(h,'Orientation','horizon')
set(h,'Interpreter','latex')

%% The triggering instants of all edges.
figure(4); 
scatter(0:T:T*(kmax-1),1*triIS12,25,'ko');
hold on;
scatter(0:T:T*(kmax-1),2*triIS21,25,'r*');
hold on;
scatter(0:T:T*(kmax-1),3*triIS13,25,'gd');
hold on;
scatter(0:T:T*(kmax-1),4*triIS31,25,'b>');
hold on;
scatter(0:T:T*(kmax-1),5*triIS14,25,'ms');
hold on;
scatter(0:T:T*(kmax-1),6*triIS41,25,'ch');
hold on;
scatter(0:T:T*(kmax-1),7*triIS15,25,'k*');
hold on;
scatter(0:T:T*(kmax-1),8*triIS51,25,'rd');
hold on;
scatter(0:T:T*(kmax-1),9*triIS16,25,'g>');
hold on;
scatter(0:T:T*(kmax-1),10*triIS61,25,'bs');
hold on;
scatter(0:T:T*(kmax-1),11*triIS23,25,'mh');
hold on;
scatter(0:T:T*(kmax-1),12*triIS32,25,'co');
hold on;
scatter(0:T:T*(kmax-1),13*triIS45,25,'kd');
hold on;
scatter(0:T:T*(kmax-1),14*triIS54,25,'r>');
hold on;
scatter(0:T:T*(kmax-1),15*triIS56,25,'cs');
hold on;
scatter(0:T:T*(kmax-1),16*triIS65,25,'bh');
hold on;
scatter(0:T:T*(kmax-1),17*triIS10,25,'mo');
hold on;
scatter(0:T:T*(kmax-1),18*triIS60,25,'g*');
hold on;
axis([0,10,0,19]);
box on;
set(gca,'yTick',[1:1:18]);
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel({'trigger instants of all edges'},'Interpreter','latex','FontSize',25);
h=legend('$e_{12}$','$e_{21}$','$e_{13}$','$e_{31}$','$e_{14}$','$e_{41}$','$e_{15}$','$e_{51}$','$e_{16}$','$e_{61}$','$e_{23}$','$e_{32}$','$e_{45}$','$e_{54}$','$e_{56}$','$e_{65}$','$e_{10}$','$e_{60}$');
set(h,'Fontsize',25);
set(h,'Interpreter','latex')

%% The control inputs $u_{i}(t)$ of six follower agents.
figure(5); 
plot(0:T:T*(kmax-1),u1(1:length(u1)),'k--','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),u2(1:length(u1)),'r--','LineWidth',3,'MarkerFaceColor','r','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),u3(1:length(u1)),'b--','LineWidth',3,'MarkerFaceColor','b','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),u4(1:length(u1)),'g--','LineWidth',3,'MarkerFaceColor','g','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),u5(1:length(u1)),'-.m','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),u6(1:length(u1)),'-.c','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$u_{i}(t)$','Interpreter','latex','FontSize',25); 
h=legend('$u_{1}$','$u_{2}$','$u_{3}$','$u_{4}$','$u_{5}$','$u_{6}$');
set(h,'Fontsize',25);
set(h,'Orientation','horizon')
set(h,'Interpreter','latex')


%% The outputs $y_{i}(t)$ and tracking errors $\zeta_{i}(t)$ of six follower agents
figure(6); 
subplot(2,1,1);
plot(0:T:T*(kmax-1),y0(1:length(y0)),'-.r','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y1(1:length(y1)),'k--','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y2(1:length(y2)),'r--','LineWidth',3,'MarkerFaceColor','r','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y3(1:length(y3)),'b--','LineWidth',3,'MarkerFaceColor','b','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y4(1:length(y4)),'-.g','LineWidth',3,'MarkerFaceColor','g','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y5(1:length(y5)),'-.m','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y6(1:length(y6)),'-.c','LineWidth',3,'MarkerFaceColor','c','MarkerSize',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$y_i(t)$','Interpreter','latex','FontSize',25); 
h=legend('Leader','Agent1','Agent2','Agent3','Agent4','Agent5','Agent6');
set(h,'Fontsize',20);
set(h,'Orientation','horizon')
subplot(2,1,2);
plot(0:T:T*(kmax-1),y1(1:length(y1))-D(1,1)*y0(1:length(y0)),'k--','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y2(1:length(y2))-D(2,2)*y0(1:length(y0)),'r--','LineWidth',3,'MarkerFaceColor','k','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y3(1:length(y3))-D(3,3)*y0(1:length(y0)),'b--','LineWidth',3,'MarkerFaceColor','r','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y4(1:length(y4))-D(4,4)*y0(1:length(y0)),'-.g','LineWidth',3,'MarkerFaceColor','b','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y5(1:length(y5))-D(5,5)*y0(1:length(y0)),'-.m','LineWidth',3,'MarkerFaceColor','g','MarkerSize',3);
hold on;
plot(0:T:T*(kmax-1),y6(1:length(y6))-D(6,6)*y0(1:length(y0)),'-.c','LineWidth',3,'MarkerFaceColor','m','MarkerSize',3);
hold on;
grid on;
set(gca,'FontSize',20,'FontName','Arial');
xlabel('time $(s)$','Interpreter','latex','FontSize',25);
ylabel('$\zeta_i(t)$','Interpreter','latex','FontSize',25); 





