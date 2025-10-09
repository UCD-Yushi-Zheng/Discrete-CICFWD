clear;clc;close all;

%---------------------------------------------------------------------------%
a1=4;
b1=3/2;
d1=1;
a2=-2;
b2=1;
d2=1;
a=1;
b=2;
d=0;

syms x
f=exp(-x^2/sqrt(0.2));
Hz_x1=20;
N=1001;

[C,u,v]=CICFWD(a1,b1,d1,a2,b2,d2,a,b,d,f,Hz_x1,N,x);

figure;
contourf(u,v,abs(C));
colorbar;
set(gca,'FontSize',20)
xlabel('u')
ylabel('v')
