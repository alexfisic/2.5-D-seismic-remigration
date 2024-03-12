close all;
clear all;

load migrado_presal.dat;

[nt nx]=size(migrado_presal);

U=migrado_presal;

dx=50.0;
dt=0.004;

U2=zeros(nt,nx); % Fourier de tec
U3=zeros(nt,nx); % Seção tecVA

U4=zeros(nt,nx); % Fourier de U (dados sismicos)
U5=zeros(nt,nx); % Hilbert de U4 (freq)
U6=zeros(nt,nx); % Hilbert de U5 (tempo)
U7=zeros(nt,nx); % Secao magnitude da reflexao
U8=zeros(nt,nx); % Secao fase instantanea
U9=zeros(nt,nx); % Secao frequencia instantanea
tec=zeros(nt,nx);
M=6;
ii=complex(0,1);

for ix=1:nx;

for it=1:nt;
    inf = it-(M/2); sup = it+(M/2);
    if(inf > 0 && sup <= nt);
    for i=inf:sup;
    tec(it,ix) = tec(it,ix) + U(i,ix)* U(i,ix);
    end;
    tec(it,ix)=sqrt(tec(it,ix)/M);
    else; end;
end;

end;

for ix=1:nx;
    
    U2(:,ix)=fft(tec(:,ix));
    U4(:,ix)=fft(U(:,ix));

    for it = 1:round(nt/2 + 2); % Transformada de Hilbert
    U2(it,ix) = -ii*U2(it,ix);
    U5(it,ix)= -ii*U4(it,ix);
    end;

    for it = round(nt/2 + 3):nt; % Transformada de Hilbert
    U2(it,ix) = ii*U2(it,ix);
    U5(it,ix) = ii*U4(it,ix);
    end;

end;

for ix=1:nx;
    U3(:,ix)=real(ifft(U2(:,ix)));
    U6(:,ix)=ifft(U5(:,ix));
end;

for ix=1:nx;
a = real(U(:,ix));  c1=diff(a); c=zeros(nt,1); c(1,1)=1.0; c(2:nt)=c1(1:nt-1);
b = real(U6(:,ix)); d1=diff(b); d=zeros(nt,1); d(1,1)=1.0; d(2:nt)=d1(1:nt-1);
    U7(:,ix)=sqrt(a.^2 + b.^2);
    U8(:,ix) = unwrap(atan(b./a));
    U9(:,ix) = (a.*d - b.*c)./(a.^2 + b.^2); 
end;

for ix=1:nx;U9(:,ix)=U9(:,ix)/max(abs(U9(:,ix)));end;

figure,imagesc([1:nx]*dx,[1:nt]*dt,-tec), colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Seção tec');

figure,imagesc([1:nx]*dx,[1:nt]*dt,U), colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Seção migrada');

figure,imagesc([1:nx]*dx,[1:nt]*dt,-U3);colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Seção tecVA');

figure,imagesc([1:nx]*dx,[1:nt]*dt,U7); axis tight;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Seção magnitude de reflexão');

figure,imagesc([1:nx]*dx,[1:nt]*dt,U8); axis tight;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Seção fase instantânea');

figure,imagesc([1:nx]*dx,[1:nt]*dt,U9);axis tight;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Seção frequência instantânea');


