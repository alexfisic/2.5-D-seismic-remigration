close all; clear all;

disp('Modelagem, migração e remigração em Offset-Finito');

nx=400;
nz=200;
ntrace=300;

dx=25;
dz=25;

nt=500;   %1000;
dt=0.008; %0.004;

z0=2500;x0=0.0;
v10=2500; vm=v10;
v2=2000;
x0s=0.0; 
x0r=50.0;
h=(x0r-x0s)/2;

g=ricker(100,10,0.004);
ng=length(g);

dado=zeros(nt,ntrace);
dado2=zeros(nt+ng-1,ntrace);
dado_esp=zeros(size(dado));

%% ===== MODELAGEM === %%

disp('Gera dados sísmicos!');

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

for ix0=1:nx
x00=0.0+(ix0-1)*dx;

if(x00>=1200 && x00<=3800);
prof=z0+(950+50)*(exp(-(x00-2500).^2/500000));
elseif(x00>=4000 && x00<=4400);
prof=z0-0.1*x00;
else;    
prof=z0;
end;

if(x00>=1200 && x00<=3800);
gx=(950+50)*(exp(-(x00-2500).^2/500000))*(-(x00-2500)/250000);
elseif(x00>=4000 && x00<=4400);
gx=-0.1;
else;    
gx=0.0;
end;

v1=v10+0.03*x00; % Variação lateral de velocidade

%v1=v10;

t1=(1/v1)*sqrt(prof^2+(xs-x00)^2);
t2=(1/v1)*sqrt(prof^2+(xr-x00)^2);
t=t1+t2;

it=round(t/dt);

lls=sqrt(prof^2+(x00-xs)^2);
llr=sqrt(prof^2+(x00-xr)^2);

ll=lls+llr;

amp=(1/(lls*llr))*(1/sqrt((v1/lls)+(v1/llr)));
D=(v10/v1^2)*lls+(0.015/v1^2)*(lls^2+llr^2)+(v10/v1^2)*llr;
amp2=(1/v1)*(1/sqrt(lls*llr))*(1/sqrt(D));

oblq=0.5*((prof/lls)+(prof/llr));

if(it<=nt && ixtrace_co<=ntrace);
dado(it,ixtrace_co)=dado(it,ixtrace_co)+amp2*oblq*sqrt(1+gx^2)*dx; else; end;

end;

if(ixtrace_co<=ntrace)
p=conv(g,dado(:,ixtrace_co));
p=halfdif2(p,nt+ng-1,dt,1);
dado2(1:nt,ixtrace_co)=p(1:nt);else;end;

end;

dado(1:nt,:)=dado2(1:nt,:);
dado_esp(:,100)=dado(:,100); % Resposta ao impulso

clear xs xr;

%% ==== MIGRAÇÃO ==== %%

disp('Migra os dados sísmicos!');

migi=zeros(nt,nx);
migi2=zeros(nt,nx);
migic=zeros(nt,nx); % Migração cinemática

for ix=1:nx; x=dx+(ix-1)*dx; 
for it=1:nt; tau=dt+(it-1)*dt;

for ixtrace=1:ntrace;
    
xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

time_1=sqrt(tau^2/4+((x-xs)^2/vm^2));
time_2=sqrt(tau^2/4+((x-xr)^2/vm^2));
time=time_1+time_2;

time21=sqrt(tau^2/4+((x-xs)^2/v2^2));
time22=sqrt(tau^2/4+((x-xr)^2/v2^2));
timee=time21+time22;

itime=round(1 +(time/dt));
itime2=round(1+(timee/dt));

wds=(tau/2/dx)*sqrt((1/time_1)+(1/time_2))*((time_1/time_2)+(time_2/time_1));
wds2=(tau/2/dx)*sqrt((1/time21)+(1/time22))*((time21/time22)+(time22/time21));

if(itime<=nt && itime2<=nt && ixtrace_co<=ntrace);
migi(it,ix)=migi(it,ix)+(1/sqrt(2*pi))*wds*dado(itime,ixtrace_co)*dx;
migic(it,ix)=migic(it,ix)+(1/(sqrt(2*pi)*dx))*dado(itime,ixtrace_co)*dx;
migi2(it,ix)=migi2(it,ix)+(1/sqrt(2*pi))*wds2*dado(itime2,ixtrace_co)*dx;
else;end;

end;

end;
end;

migi2_esp=zeros(size(migi));
migi2_esp(:,40)=migi2(:,40);

clear xs xr;

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

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2; 
ixtrace_co=1+round(xtrace/dx);

factor=vm^2-v2^2;
extra=4*(h^2)*((1/vm^2)-(1/v2^2));
radicando1=tau^2/4+((x-xtrace)^2/factor)+extra;
radicando2=tau^2/4+((x-xtrace)^2/factor)+extra;

if(radicando1>=0 && radicando2>=0);
teta1=sqrt(radicando1);
teta2=sqrt(radicando2);
teta=teta1+teta2;
else;teta=0;teta1=0.5;teta2=0.5;end;

iteta=round(1+(teta/dt));

estreti=(tau/4/radicando1)+(tau/4/radicando2);

deno_velo=1-(((v2/vm)^2)*(tau/teta));
nume_velo=1;
novo=nume_velo/deno_velo;

w1=(tau/4/v2^(3/2))*sqrt(2/tau)*2;
w2=(1/(teta1+teta2))*sqrt((1/teta1)+(1/teta2));
wrm=w1*w2*(novo); %*sqrt(estreti);

if(iteta<=nt && ixtrace_co<=nx);
migi3(it,ix)=migi3(it,ix)+wrm*migi2(iteta,ixtrace_co)*dx;
else;end;

end;

end;
end;

%% == TRANSFORMAÇÃO TEMPO-PROFUNDIDADE == %

MIG=zeros(nz,nx);
MIG2=zeros(nz,nx);

for ix=1:nx;
for it=1:nt;

  t=dt+(it-1)*dt;
  z=0.5*vm*t; z2=0.5*v2*t;
  iz=1+round(z/dz); iz2=1+round(z2/dz);

if(iz<=nz && iz2<=nz);
MIG(iz,ix)=migi3(it,ix);MIG2(iz2,ix)=migi2(it,ix);
end;

end;
end;

f=real(MIG);
save ttd.dat f -ascii;
clear f;

figure;
imagesc([1:ntrace]*dx,[1:nt]*dt,-real(dado));grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados sísmicos sintético');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi));grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Sintético migrado - Kirchhoff');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi2));grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Sintético submigrado - Kirchhoff');

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi3));grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Sintético remigrado');

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[1:nz]*dz,real(MIG));colormap gray;grid;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title(' Transformação tempo-profundidade - remigrado');
subplot(1,2,2);
imagesc([1:nx]*dx,[1:nz]*dz,real(MIG2)); colormap gray;grid;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title(' Transformação tempo-profundidade - velocidade errada');