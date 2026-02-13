clear;clc;close all;

%---------------------------Parameter use for CICFWD calculations------------------------------%
a=1;
b=2;
d=0;
a1=3;
b1=5;
d1=5;
a2=-1;
b2=1/2;
d2=-4;
W1=13;
W2=3.5;
Wt=4;

%---------------------------Input signal------------------------------%
syms x
f=exp(-x^2/sqrt(0.2));
Hz_x1=10;
N=1501;

%---------------------------CICFWD based on global bandwidth------------------------------%
[C,u,v]=CICFWD_global(a1,b1,d1,a2,b2,d2,a,b,d,f,Hz_x1,N,x);

figure;
contourf(u,v,abs(C));
colorbar;
set(gca,'FontSize',20)
xlabel('u')
ylabel('v')
xlim([-5,5]);
ylim([-15,15]);


%---------------------------CICFWD based on local bandwidth------------------------------%
f=exp(-x^2/sqrt(0.2));
Hz_x1=10;
N=1501;

[C,u,v]=CICFWD_local(a1,b1,d1,a2,b2,d2,a,b,d,f,Hz_x1,N,x,W1,W2,Wt);

figure;
contourf(u,abs(b)*v,abs(C));
colorbar;
set(gca,'FontSize',20)
xlabel('u')
ylabel('v')
xlim([-5,5]);
ylim([-15,15]);
