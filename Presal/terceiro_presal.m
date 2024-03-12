close all; clear all;

load kir_presal.dat;
load vel_migrat_presal.dat;
load velo_smooth_fourier_presal2.dat;
load vel_migrat_2009.dat;

q=reshape(kir_presal,1750,720);

R=vel_migrat_presal';
T=1000*vel_migrat_2009;
%T=velo_smooth_fourier_presal2';

[nt nx]=size(q);

dx=50.0;
dt=0.004;
nz=160;
dz=dx;
ntrace=nx;

x0s=0.0; 
x0r=50.0;
h=(x0r-x0s)/2;

secout=ormsby(nt,nx,dt,dx,5,8,15,25,q);

for ix=1:nx;
p1=secout(:,ix);
p1=halfdif2(p1,nt,dt,1);
secout(1:nt,ix)=p1(1:nt);
end;

dado=secout;

clear q;

%%== Zeros de R e T ==%%
%for ix=1:200;
%indice=find(R(:,ix)==0);
%ni=length(indice);
%R(indice(1,1):indice(ni,1),ix)=4500;
%T(indice(1,1):indice(ni,1),ix)=4500;
%end;

%% ==== MIGRAÇÃO ==== %%

disp('Migra os dados sísmicos!');

migi=zeros(nt,nx);
migi2=zeros(nt,nx);
migic=zeros(nt,nx); % Migração cinemática

amigi=zeros(nt,nx);  % Seção alias do Kirchhoff
amigi2=zeros(nt,nx); % Seção alias do subgmigrado
amigic=zeros(nt,nx); % Seção alias da migração cinemática

dipmig=zeros(1,nx);  % Vetor de pseudo-mergulho Kirchhoff (difração)

for it=1:nt; tau=dt+(it-1)*dt;
for ix=1:nx; x=dx+(ix-1)*dx; 

vm=R(it,ix);
v2=T(it,ix)/2;

for ixtrace=1:ntrace;
    
xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

time_1=sqrt(tau^2/4+((x-xs)^2/vm^2));
time_2=sqrt(tau^2/4+((x-xr)^2/vm^2));
time=time_1+time_2;

ldip1=(dx/vm^2)*(((x-xs)/time_1)+((x-xr)/time_2));

dipmig(1,ixtrace)=ldip1;

time21=sqrt(tau^2/4+((x-xs)^2/v2^2));
time22=sqrt(tau^2/4+((x-xr)^2/v2^2));
timee=time21+time22;

ldip2=(dx/v2^2)*(((x-xs)/time21)+((x-xr)/time22));

itime=round(1 +(time/dt));
itime2=round(1+(timee/dt));

wds=(tau/2/dx)*sqrt((1/time_1)+(1/time_2))*((time_1/time_2)+(time_2/time_1));
wds2=(tau/2/dx)*sqrt((1/time21)+(1/time22))*((time21/time22)+(time22/time21));

if(itime<=nt && ixtrace_co<=ntrace && ldip1>0.041 | ldip1<-0.041);
amigi(it,ix)=amigi(it,ix)+(1/sqrt(2*pi))*wds*dado(itime,ixtrace_co)*dx;
amigic(it,ix)=amigic(it,ix)+(1/(sqrt(2*pi)*dx))*dado(itime,ixtrace_co)*dx;
else;end;

if(itime2<=nt && ixtrace_co<=ntrace && ldip2>0.12 | ldip2<-0.12);
amigi2(it,ix)=amigi2(it,ix)+(1/sqrt(2*pi))*wds2*dado(itime2,ixtrace_co)*dx;
else;end;

if(itime<=nt && itime2<=nt && ixtrace_co<=ntrace && ldip1<0.041 | ldip1>-0.041);
migi(it,ix)=migi(it,ix)+(1/sqrt(2*pi))*wds*dado(itime,ixtrace_co)*dx;
migic(it,ix)=migic(it,ix)+(1/(sqrt(2*pi)*dx))*dado(itime,ixtrace_co)*dx;
else;end;

if(itime2<=nt && ixtrace_co<=ntrace && ldip2<0.12 | ldip2>-0.12);
migi2(it,ix)=migi2(it,ix)+(1/sqrt(2*pi))*wds2*dado(itime2,ixtrace_co)*dx;
else;end;

end;

end;
end;

migi2_esp=zeros(size(migi));
migi2_esp(:,360)=migi2(:,360);

clear xs xr;

%secout=ormsby(nt,nx,dt,dx,5,10,15,25,real(migi2));

%migi2=secout;

migi=migi-amigi;
migic=migic-amigic;
migi2=migi2-amigi2;

%% ==== REMIGRAÇÃO ==== %%

disp('Remigra os dados sísmicos!');

migi3=zeros(nt,nx);
amigi3=zeros(nt,nx); % Seção alias da remigração

for ix=1:nx;
p2=migi2(:,ix);
p2=halfdif2(p2,nt,dt,1);
migi2(1:nt,ix)=p2(1:nt);
end;

for ix=1:nx; x=dx+(ix-1)*dx; ix 
for it=1:nt; tau=dt+(it-1)*dt;

vm=R(it,ix);
v2=T(it,ix)/2;

if(vm~=0 && v2~=0); % LAÇO DA VELOCIDADE NULA

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2; 
ixtrace_co=1+round(xtrace/dx);

factor=vm^2-v2^2;

if(factor>=0);
teta1=sqrt(tau^2/4+((x-xtrace)^2/factor));
teta2=sqrt(tau^2/4+((x-xtrace)^2/factor));
teta=teta1+teta2;
else;
teta1=real(sqrt(tau^2/4+((x-xtrace)^2/factor)));
teta2=real(sqrt(tau^2/4+((x-xtrace)^2/factor)));
teta=real(teta1+teta2);
end;

iteta=round(1+(teta/dt));

deno_velo=1-(((v2/vm)^2)*(tau/teta));
nume_velo=1;
novo=nume_velo/deno_velo;

w1=(tau/4/v2^(3/2))*sqrt(2/tau)*2;
w2=(1/(teta1+teta2))*sqrt((1/teta1)+(1/teta2));
wrm=w1*w2*(novo); %*(teta^(2.21));

ldipr=dipmig(1,ixtrace);

if(iteta<=nt && ixtrace_co<=nx && ldipr>0.13 | ldipr<-0.13);
amigi3(it,ix)=amigi3(it,ix)+(1/sqrt(2*pi))*wrm*migi2(iteta,ixtrace_co)*dx;
else;end;

if(iteta<=nt && ixtrace_co<=nx && ldipr<0.13 | ldipr>-0.13);
migi3(it,ix)=migi3(it,ix)+wrm*migi2(iteta,ixtrace_co)*dx;
else;end;

end;

else;migi3(it,ix)=0.0;end; % FIM DO LAÇO DA VELOCIDADE NULA

end;
end;

migi3=migi3-amigi3;

%% == TRANSFORMAÇÃO TEMPO-PROFUNDIDADE == %

MIG=zeros(nz,nx);
MIG2=zeros(nz,nx);

for ix=1:nx;
for it=1:nt;

vm=R(it,ix);
v2=T(it,ix)/2;

  t=dt+(it-1)*dt;
  z=0.5*vm*t; z2=0.5*v2*t;
  iz=1+round(z/dz); iz2=1+round(z2/dz);

if(iz>0 && iz<=nz && iz2> 0 && iz2<=nz);
MIG(iz,ix)=migi3(it,ix);MIG2(iz2,ix)=migi2(it,ix);
end;

end;
end;

f=real(MIG);
save ttd_presal.dat f -ascii;
clear f;

k=real(migi);
save migrado_presal.dat k -ascii;
clear k;

figure;
imagesc([1:ntrace]*dx,[1:nt]*dt,real(dado));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados sísmicos Pré-Sal');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal migrado - Kirchhoff');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migic));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal migrado - Kirchhoff (cinemático)');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi2));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal submigrado - Kirchhoff');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi3));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal remigrado');

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[1:nz]*dz,real(MIG));colormap gray;grid;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title(' Transformação tempo-profundidade - remigrado');
subplot(1,2,2);
imagesc([1:nx]*dx,[1:nz]*dz,real(MIG2)); colormap gray;grid;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title(' Transformação tempo-profundidade - velocidade errada');