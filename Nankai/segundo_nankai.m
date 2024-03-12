close all; clear all;

load nankai_gain.dat;
load vel_nankai_RMS.dat;
load vel_nankai_smooth.dat;

[nt nx]=size(nankai_gain);

dx=16.67;
dt=0.004;
nz=601;
dz=dx;

dado=nankai_gain;

%% ==== MIGRAÇÃO ==== %%

disp('Migra os dados sísmicos!');

migi=zeros(nt,nx);
migi2=zeros(nt,nx);

for ix=1:nx; x=dx+(ix-1)*dx; 
for it=1:nt; tau=dt+(it-1)*dt;

vm=vel_nankai_RMS(it,ix);
v2=vel_nankai_smooth(it,ix)/2;

for ixtrace=1:nx;xt=dx+(ixtrace-1)*dx;

time=sqrt(tau^2+(4*(x-xt)^2/vm^2));
time2=sqrt(tau^2+(4*(x-xt)^2/v2^2));

itime=round(1 +(time/dt));
itime2=round(1+(time2/dt));

wds=2*sqrt(time);

if(itime<=nt && itime2<=nt);
migi(it,ix)=migi(it,ix)+wds*dado(itime,ixtrace)*dx;
migi2(it,ix)=migi2(it,ix)+wds*dado(itime2,ixtrace)*dx;
else;end;

end;

end;
end;

%% ==== REMIGRAÇÃO ==== %%

disp('Remigra os dados sísmicos!');

migi3=zeros(nt,nx);

for ix=1:nx;
p2=migi2(:,ix);
p2=halfdif2(p2,nt,dt,1);
migi2(1:nt,ix)=p2(1:nt);
end;

for ix=1:nx; x=dx+(ix-1)*dx; ix 
for it=1:nt; tau=dt+(it-1)*dt;

vm=vel_nankai_RMS(it,ix);
v2=vel_nankai_smooth(it,ix)/2;

for ixtrace=1:nx;xz=dx+(ixtrace-1)*dx;

factor=vm^2-v2^2;

if(factor>=0);
teta=sqrt(tau^2+(4*(x-xz)^2/factor));
else;teta=0;end;
iteta=round(1+(teta/dt));

% wrm=1/(v2*teta); % Função-peso original

deno_velo=1-(((v2/vm)^2)*(tau/teta));
nume_velo=1;
novo=nume_velo/deno_velo;
wrm=(2/v2)*(tau/teta)*(1/sqrt(tau*teta))*(novo)*(teta^(1.21));

if(iteta<=nt);
migi3(it,ix)=migi3(it,ix)+wrm*migi2(iteta,ixtrace)*dx;
else;end;

end;

end;
end;

%% == TRANSFORMAÇÃO TEMPO-PROFUNDIDADE == %

MIG=zeros(nz,nx);
MIG2=zeros(nz,nx);

for ix=1:nx;
for it=1:nt;

vm=vel_nankai_RMS(it,ix);
v2=vel_nankai_smooth(it,ix)/2;

  t=dt+(it-1)*dt;
  z=0.5*vm*t; z2=0.5*v2*t;
  iz=1+round(z/dz); iz2=1+round(z2/dz);

if(iz<=nz&&iz2<=nz);
MIG(iz,ix)=migi3(it,ix);MIG2(iz2,ix)=migi2(it,ix);
end;

end;
end;

f=real(MIG);
save ttd_nankai.dat f -ascii;
clear f;

outi=real(migi3);
save nankai_remig.dat outi -ascii;
clear outi;

figure;
imagesc([1:nx]*dx,[1300:2000]*dt,-real(dado(1300:2000,:)));colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados sísmicos Nankai');

figure;
imagesc([1:nx]*dx,[1300:2000]*dt,real(migi(1300:2000,:)));colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Nankai migrado - Kirchhoff');
figure;
imagesc([1:nx]*dx,[1300:2000]*dt,real(migi2(1300:2000,:)));colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Nankai migrado - Kirchhoff');
figure;
imagesc([1:nx]*dx,[1300:2000]*dt,real(migi3(1300:2000,:)));colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Nankai remigrado');

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[200:601]*dz,real(MIG(200:601,:)));colormap gray;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title(' Transformação tempo-profundidade - remigrado');
subplot(1,2,2);
imagesc([1:nx]*dx,[200:601]*dz,real(MIG2(200:601,:)));colormap gray;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title(' Transformação tempo-profundidade - velocidade errada');

