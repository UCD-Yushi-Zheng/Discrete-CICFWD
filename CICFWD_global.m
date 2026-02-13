
%---------------------------numerical CICFWD computation algorithm based on global bandwidth------------------------------%
function [result,u,v] = CICFWD_global(a1,b1,d1,a2,b2,d2,a,b,d,f,Hz_x1,N,x)

Hz_x2=Hz_x1*abs(b1)/abs(b2);

x0=-(N-1)/2:1:(N-1)/2;
T1=1./(Hz_x1);
T2=1./(Hz_x2);
x1=x0*T1;
x2=x0*T2;

u=x0*(abs(b1)/(N*T1));

Tu=u(2)-u(1);

y1=double(subs(f,x,x1));
y2=double(subs(f,x,x2));

F1 = direct_method1d(y1,T1, a1, b1, d1)*T1;
F2 = direct_method1d(y2,T2, a2, b2, d2)*T2;

v=x0*(abs(b)/(N*Tu*2));

[result,~]=BDFofLCT_function(F1,F2,Tu,a,b,d);

result=result*Tu*2;



end