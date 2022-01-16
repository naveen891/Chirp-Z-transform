function cz=direct_chirp(x, M, w, a)

    x=x.'; 
    [N, q] = size(x); 

L=2^ceil(log(N+M-1)/log(2));

maxl=max(M-1,N-1);
n=((-N+1):maxl).';
nmat=(n.^ 2)./ 2;

wmat=w.^(nmat); 

nn=(0: N-1)';
amat=a.^( -nn );
awmat=amat.*wmat(N+nn);
y=x.*awmat;
v=1./wmat(1:(M-1+N));

dft1=fft(y, L);
dft2=fft(v, L ); 
conv=dft1.*dft2;
g=ifft(conv);

var=N:N+M-1;
cz=wmat(var).*g(var);

end
