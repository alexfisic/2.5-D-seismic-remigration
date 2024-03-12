close all; clear all;

nx=200; 
nz=100;
nt=800;
dx=25;
dt=0.005;
p1=.00017;
v=2000/1;
v2=v*v;
t0=0.375; t02=2*t0;
x0=0; xc=2500;

g=ricker(100,10,dt);
ng=length(g);

data1=zeros(ng+nt-1,nx);
data=zeros(nt,nx);
data_esp=data;
data0=data;
p=zeros(1,nx);          % Vetor vagarosidade

for ix=1:nx;
    x=x0+(ix-1)*dx;
    t=sqrt(t02*t02+(x-xc)*(x-xc)*(4/v2)); % Tem um valor 4 aqui esquisito!!!
    p(ix)=x/(t*v2);
    t1=p1*abs(x-xc);
    it=1+round(t/dt);
    it1=1+round(t1/dt); 
    it0=1+round(5*t0/dt); it01=1+round(sqrt(t02^2+2500*2500/v2)/dt);
    
    if(it <= nt);
    data(it,ix)=1; else;end;
    if(it1 <= nt); 
    data0(it1,ix)=1;else;end;
    data0(it0,ix)=1;  data0(it01,100)=1;
    data1(:,ix)=conv(data(:,ix)+data0(:,ix),g); 
    
end;

data(1:nt,:)=data1(1:nt,:);
data_esp(:,100)=data(:,100);

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
kx(nx/2+2)=-(nx/2-1)*dk;
w(1)=0.0;

for ix=2:nx/2+1;
   kx(ix)= kx(ix-1)+dk;
end;

for ix=nx/2+3:nx;
    kx(ix)=kx(ix-1)+dk;
end;

for it=2:nt/2+1;
    w(it)=w(it-1)+dw;
end;

ii=complex(0,1);

%====================================================%

for ix=1:nx;
a=real(data3(:,ix));
b=imag(data3(:,ix));
data4(:,ix)=sqrt(a.^2 + b.^2);
end;

data6(:,1:nx/2)=data4(:,(nx/2)+1:nx);
data6(:,(nx/2)+1:nx)=data4(:,1:nx/2);

for i=1:nx;data6(:,i)=data6(:,i)/max(abs(data6(:,i)));end;

% ============== Filtro de mergulho ================= % 

dipdata=zeros(size(data6)); pin=p(110);

for it=1:round(nt/12)+1;
kapa=w(it)*pin;
ikapa=1+round(kapa/dk);
if(ikapa <= (nx/2)+1);
dipdata(it,ikapa)=data3(it,ikapa); else; end;
end;

data5=dipdata; ep=1e-03;

%for i=1:ikapa;dipdata(1:74,i)=dipdata(1:74,i)/(ep+max(abs(dipdata(:,i))));end;

%=============== Transformada iFK =====================%

for ix=1:nx;
    dipdata(:,ix)=real(ifft(dipdata(:,ix)));
end;

for it=1:nt;
    dipdata(it,:)=real(ifft(dipdata(it,:))); 
end;

%====================== Figuras ======================%

figure,imagesc(real(data3)); colormap prism;

figure;
subplot(1,2,1);
imagesc([-kN:dk:kN-dk],[1:nt/2]*(dw/(2*pi)),data6(1:nt/2,:)); colorbar;
subplot(1,2,2);
imagesc(data);
