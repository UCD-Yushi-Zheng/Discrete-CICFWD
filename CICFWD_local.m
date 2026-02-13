%---------------------------numerical CICFWD computation algorithm based on local bandwidth------------------------------%
function [result,ux,uf_mod] = CICFWD_local(a1,b1,d1,a2,b2,d2,a,b,d,f,Hz_x1,N,x,W1,W2,Wt)

Hz_x2=Hz_x1*abs(b1)/abs(b2);

x0=-(N-1)/2:1:(N-1)/2;
T1=1./(Hz_x1);
T2=1./(Hz_x2);
x1=x0*T1;
x2=x0*T2;

ux=x0*(abs(b1)/(N*T1));

Tu=ux(2)-ux(1);

y1=double(subs(f,x,x1));
y2=double(subs(f,x,x2));

%-------------------------------------------------------------------------------%
center=zeros(1,length(ux));
for n=1:1:length(ux)
    if (W2-W1)/4<=ux(n) && ux(n)<=(W1-W2)/4
        center(n)=(d1/b1+2*a/b)*ux(n);
    elseif  (W1-W2)/4<=ux(n) && (W2+W1)/4>=ux(n)
        center(n)=(d1/(2*b1)+d2/(2*b2))*ux(n)+(d1/(4*b1)-d2/(4*b2)+a/b)*(W1-W2)/2;
    elseif (W2-W1)/4>=ux(n) && -(W1+W2)/4<=ux(n)
        center(n)=-((d1/(2*b1)+d2/(2*b2))*(-ux(n))+(d1/(4*b1)-d2/(4*b2)+a/b)*(W1-W2)/2);
    end
end

W=zeros(1,length(ux));
for n=1:1:length(ux)
    if (W2-W1)/4 <= ux(n) && ux(n) <= (W1-W2)/4
        W(n)= 2*W2*(d1/(4*b1)-d2/(4*b2)) + Wt/(2*abs(b2)) + Wt/(2*abs(b1));
    else
        if (d1/(4*b1)-d2/(4*b2)+a/b)*(W1+W2-4*abs(ux(n)))+Wt/(2*abs(b2))+Wt/(2*abs(b1))>0
            W(n)=(d1/(4*b1)-d2/(4*b2)+a/b)*(W1+W2-4*abs(ux(n)))+Wt/(2*abs(b2))+Wt/(2*abs(b1));
        end
    end
end

%-------------------------------------------------------------------------------%

uf=(x0/(N*Tu*2));

Tuf=uf(2)-uf(1);
L=3*(max(center+W));
Nuf=L/Tuf;
uf_mod=(-(Nuf-1)/2:1:(Nuf-1)/2)*Tuf;

F1 = direct_method1d(y1,T1, a1, b1, d1)*T1;
F2 = direct_method1d(y2,T2, a2, b2, d2)*T2;

%---------------------------------------------------------------------------------%
len = N;
fext1 = [zeros(1,length(F1)),F1,zeros(1,length(F1))];
fext2 = [zeros(1,length(F2)),F2,zeros(1,length(F2))];
for n=1:1:len % Caculate the instantaneous auto-correlation function (IAF)
    for j = -(len-1)/2:1:(len-1)/2
        acorr(n,j+(len-1)/2+1) = fext1(n+len+j)*conj(fext2(n+len-j));
    end
end

gacorr=acorr';
w=zeros(length(uf_mod),len);
T=2*Tu;

for n=1:length(ux) % Calculate the lct of each line in the IAF
    centerA=center(n)-(ceil(center(n)/(1/T))-1)*(1/T);
    centerB=center(n)-(ceil(center(n)/(1/T)))*(1/T);
    if abs(centerA)<abs(centerB)
        center0=centerA;
    elseif abs(centerB)<abs(centerA)
        center0=centerB;
    elseif abs(centerB)==abs(centerA)
        center0=center;
    end
    idx=find(abs(uf-center0)==min(abs(uf-center0)));
    nSamp = length(gacorr(:,n)');   
    row = -(nSamp-1)/2:1:(nSamp-1)/2; 
    chirp1 = exp(1i*0.5*2*pi*((row*T).^2)*a/b);    
    FTchirp=fftshift(fft(gacorr(:,n)'.*chirp1));
    if center0 < 0
        F_mod=[FTchirp(end-((len-1)/2-(idx-1))+1:end), FTchirp(1:end-((len-1)/2-(idx-1)))];
    elseif center0>0
        F_mod=[FTchirp((len-1)/2-(len-idx)+1:end ), FTchirp( 1: (len-1)/2-(len-idx))];
    elseif center0==0
        F_mod=FTchirp;
    end
    ind=find(abs(uf_mod-center(n))==min(abs(uf_mod-center(n))));
    F_shift=zeros(1,length(uf_mod));
    F_shift(ind-(len-1)/2:ind+(len-1)/2)=F_mod;
    lct = F_shift;
    w(:,n)=lct;
end
result = w*T;


end