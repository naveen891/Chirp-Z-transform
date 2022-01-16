%Note that we have taken a direct example from internet from reference
%mentioned in report
x=fir1(30,125/500,rectwin(31));
fs=1000; 
f1=100;
f2=150; 
M=1024;
w=exp(-j*2*pi*(f2-f1)/(M*fs));
a=exp(j*2*pi*f1/fs);

z0=fft(x);
z=direct_chirp(x,M,w,a);
z1=czt(x, M, w, a);
fz=((0:length(z)-1)'*(f2-f1)/length(z)) + f1;

figure();
subplot(2, 1, 1) 
plot(fz,abs(z));
axis([f1 f2 0 1.2])
title('Czt obtained from our function');

subplot(2, 1, 2) 
plot(fz,abs(z1));
axis([f1 f2 0 1.2])
title('Czt obtained from inbuilt-function');

%-------------------------------------------------

M=64;
m=0:M-1;

x=sin(2*pi*m/15);
FFT=fft(x);
CZT=direct_chirp(x, M, exp(-2j*pi/M), 1);
figure();
subplot(2, 1, 1)
plot(m, abs(FFT));
title("FFT");
subplot(2, 1, 2);
plot(m,abs(CZT));
title("CZT for circular contour");

%--------------------------------------------

var=1000:100:8000
tczt=zeros(numel(var));
tdirectczt=zeros(numel(var));

for k = 1:numel(var)
    
        f=0:var(k)-1;
        x=sin(2*pi*f/15);
        tic;
        CZT=czt(x, var(k), exp(-2j*pi/var(k)), 1);
        tczt(k) = toc ;
        
        tic;
        DirectCZT=direct_chirp(x, var(k), exp(-2j*pi/var(k)), 1);
        tdirectczt(k) = toc ;
    
end
figure ;
subplot(2,1,1) ;
plot(var, tczt) ;
title("Time taken for inbuilt function") ;
xlabel("N");
ylabel("t");

subplot(2,1,2) ;
plot(var, tdirectczt) ;
title("Time taken for  our function") ;
xlabel("N");
ylabel("t");

%The value of the functions obtained is almost the same.
%-----------------------------
%Other than the time complexity part, we were asked to  include any one of the application.
%Due to lack of resources about implementing the applications in 
%the internet or in the paper and lack of time we failed to implement it.
