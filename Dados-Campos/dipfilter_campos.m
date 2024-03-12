close all; clear all;

load remigrado_0258_7dias_NaN.dat;

Q=remigrado_0258_7dias_NaN;
Q(isnan(Q))=0.0;

giga(:,1:1771)=Q(:,1:1771);
%giga(isnan(giga))=0.0;

nt=1750;
nx=1771;

dt=0.004;
dx=12.52;

data=zeros(nt,nx);

%for ix=1:nx;giga(:,ix)=giga(:,ix)/max(abs(giga(:,ix)));end;

data=giga;

%clear Q; clear giga;

pin=[5e-4 0.0 5e-4];

%=============== Transformada FK =====================%

for ix=1:nx;
    data2(:,ix)=fft(data(:,ix));
end;

for it=1:nt;
    data3(it,:)=fft(data2(it,:)); 
end;

%=====================================================%

dw=(2*pi)/(nt*dt);
fN=(2*pi)/(2*dt);       % Frequencia Nyquist

dk=(2*pi)/(nx*dx);
kN=(2*pi)/(2*dx);       % Frequencia Nyquist espacial

kx=zeros(1,nx);         % Vetor frequencia espacial
w=zeros(1,nt);          % Vetor frequencia temporal

kx(1)=0.0;
kx(round(nx/2)+2)=-(nx/2-1)*dk;
w(1)=0.0;

for ix=2:round(nx/2)+1;
   kx(ix)= kx(ix-1)+dk;
end;

for ix=round(nx/2)+3:nx;
    kx(ix)=kx(ix-1)+dk;
end;

for it=2:round(nt/2)+1;
    w(it)=w(it-1)+dw;
end;

ii=complex(0,1);

%====================================================%
for ix=1:nx;
a=real(data3(:,ix));
b=imag(data3(:,ix));
data4(:,ix)=sqrt(a.^2 + b.^2);
end;

data6(:,1:(nx/2))=data4(:,(nx/2)+1:nx);
data6(:,(nx/2)+1:nx)=data4(:,1:(nx/2));

for i=1:nx;data6(:,i)=data6(:,i)/max(abs(data6(:,i)));end;
% ============== Filtro de mergulho ================= % 

dipdata=zeros(size(data)); 
dip=zeros(size(data));

for ip=1:3;

for it=2:round(nt/2)+1;
kapa=w(it)*pin(ip);
ikapa=1+round(kapa/dk);
if(ikapa <= round(nx/2)+1);
dipdata(it,ikapa)=data3(it,ikapa); else; end;
end;

dip=dip+dipdata;

end;

data5=dip; ep=1e-03;

%=============== Transformada iFK =====================%

for ix=1:nx;
    dip(:,ix)=real(ifft(dip(:,ix)));
end;

for it=1:nt;
    dip(it,:)=real(ifft(dip(it,:))); 
end;

%====================== Figuras ======================%

for ix=1:nx;dip(:,ix)=dip(:,ix)/max(abs(dip(:,ix)));end;

for ix=1:nx;data(:,ix)=data(:,ix)/max(abs(data(:,ix)));end;

data=data-dip;

figure,imagesc(real(data)); colormap gray;

figure,imagesc(real(dip)); colormap gray;

figure;
subplot(1,2,1);
imagesc([-kN:dk:kN-dk],[1:nt/6]*(dw/(2*pi)),data6(1:nt/6,:)); colorbar;
subplot(1,2,2);
imagesc(data);


