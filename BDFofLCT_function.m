function [w,gacorr] = BDFofLCT_function(F1,F2,T,a,b,d)

len = length(F1);
fext1 = [zeros(1,length(F1)),F1,zeros(1,length(F1))];
fext2 = [zeros(1,length(F2)),F2,zeros(1,length(F2))];
for i=1:1:len % Caculate the instantaneous auto-correlation function (IAF)
    for j = -(len-1)/2:1:(len-1)/2
        % disp(['A-fft Trial ' num2str(i) ' of ' num2str(len)])
        acorr(i,j+(len-1)/2+1) = fext1(i+len+j)*conj(fext2(i+len-j));
    end
end

gacorr=acorr';
w=gacorr;
for i=1:len % Calculate the lct of each line in the IAF 
    % disp(['A-fft Trial ' num2str(i) ' of ' num2str(len)])
    w(:,i) = direct_method1d(gacorr(:,i)',2*T, a, b, d);
end
end