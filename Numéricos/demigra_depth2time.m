close all; clear all;

load sin_2D_transp.dat;
load dado_2D_tempo.dat;

[nz nx]=size(sin_2D_transp);

ntrace=nx;
ntime=1000;

dx=25.;
dz=25.;
dt=0.004;

x0s=0.0;
x0r=50.;

h=(x0r-x0s)/2;

v1=2500;

saida=zeros(ntime,ntrace);

mig=sin_2D_transp;
dado=dado_2D_tempo;

for ixtrace=1:ntrace; ixtrace

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;

ixtrace_co=1+round(xtrace/dx);

for it=1:ntime;

tau=dt+(it-1)*dt;

limite=(v1/2)*tau;
nl=round(limite/dz);

linf=ixtrace_co-round(nl/1);
lup=ixtrace_co+round(nl/1);

if(linf<1);linf=1;else;end;
if(lup>ntrace);lup=ntrace;else;end;

if(limite^2>h^2); % Laço do tempo

for ix2=linf:lup; % Laço mais interno

x2=ix2*dx;    

if((x2-xtrace)^2<limite^2);
zeta=sqrt(limite^2-h^2)*sqrt(1-((x2-xtrace)^2/limite^2));
else;
zeta=0.0;    
end;

izeta=round(1+(zeta/dz));

t1=(1/v1)*sqrt(zeta^2+(x2-xs)^2);
t2=(1/v1)*sqrt(zeta^2+(x2-xr)^2);

wis=(1/2/v1)*(1/(t1+t2))*sqrt((1/t1)+(1/t2));

if(izeta<=nz && ixtrace_co<=ntrace);
saida(it,ixtrace_co)=saida(it,ixtrace_co)+wis*mig(izeta,ix2)*dx;    
else;end;

end; % Fim do laço mais interno

else;saida(it,ixtrace_co)=0.0;end; % Fim laço do tempo

end;

end;

figure;
imagesc([1:ntrace]*dx,[1:ntime]*dt,-real(saida));colormap gray; grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dado sísmico demigrado');

figure;
imagesc([1:ntrace]*dx,[1:ntime]*dt,-real(dado));colormap gray; grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dado sísmico original');

figure;
imagesc([1:nx]*dx,[1:nz]*dz,-real(mig)); colormap gray;grid;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title('Dado sísmico em profundidade');





