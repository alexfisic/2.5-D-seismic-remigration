close all; clear all;

nx=200;
nz=200;


dx=25;
dz=25;

nt=500;   %1000;
dt=0.008; %0.004;

z0=2500;x0=0.0;
v10=2500; vm=v10;
v2=2000;

g=ricker(100,10,0.004);
ng=length(g);

dado=zeros(nt,nx);
dado2=zeros(nt+ng-1,nx);
dado_esp=zeros(size(dado));

%% ===== MODELAGEM === %%

disp('Gera dados sísmicos!');

for ix=1:nx;
x=x0+(ix-1)*dx;

for ix0=1:nx
x00=0.0+(ix0-1)*dx;

if(x00>=1200 && x00<=3800);
prof=z0-(950+50)*(exp(-(x00-2500).^2/300000));
else;
prof=z0;
end;

%prof=z0-0.1*x00;

v1=v10+0.03*x00; % Variação lateral de velocidade

% v1=v10;

t=(2/v1)*sqrt(prof^2+(x-x00)^2);
it=round(t/dt);

ll=sqrt(prof^2+(x00-x)^2);
amp=(1/ll^2)*(1/sqrt(2/ll));
oblq=z0/(ll*v1);

if(it<=nt);
dado(it,ix)=dado(it,ix)+amp*oblq*1.0*dx; else; end;

end;

p=conv(g,dado(:,ix));
p=halfdif2(p,nt+ng-1,dt,1);
dado2(1:nt,ix)=p(1:nt);

end;

dado(1:nt,:)=dado2(1:nt,:);
dado_esp(:,100)=dado(:,100); % Resposta ao impulso

%% ==== MIGRAÇÃO ==== %%

disp('Migra os dados sísmicos!');

migi=zeros(nt,nx);
migi2=zeros(nt,nx);

for ix=1:nx; x=dx+(ix-1)*dx; 
for it=1:nt; tau=dt+(it-1)*dt;

for ixtrace=1:nx;xt=dx+(ixtrace-1)*dx;

time=sqrt(tau^2+(4*(x-xt)^2/vm^2));
time2=sqrt(tau^2+(4*(x-xt)^2/v2^2));

itime=round(1 +(time/dt));
itime2=round(1+(time2/dt));

wds=2*sqrt(time);

if(itime<=nt && itime2<=nt);
migi(it,ix)=migi(it,ix)+wds*dado_esp(itime,ixtrace)*dx;
migi2(it,ix)=migi2(it,ix)+wds*dado_esp(itime2,ixtrace)*dx;
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
wrm=(2/v2)*(tau/teta)*(1/sqrt(tau*teta))*(novo);

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

  t=dt+(it-1)*dt;
  z=0.5*vm*t; z2=0.5*v2*t;
  iz=1+round(z/dz); iz2=1+round(z2/dz);

if(iz<=nz&&iz2<=nz);
MIG(iz,ix)=migi3(it,ix);MIG2(iz2,ix)=migi2(it,ix);
end;

end;
end;

f=real(MIG);
save ttd.dat f -ascii;
clear f;

figure;
imagesc([1:nx]*dx,[1:nt]*dt,-real(dado));

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi));
figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi2));
figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi3));

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[1:nz]*dz,real(MIG));colormap gray;
subplot(1,2,2);
imagesc([1:nx]*dx,[1:nz]*dz,real(MIG2)); colormap gray;


