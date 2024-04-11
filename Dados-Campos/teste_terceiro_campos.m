close all; clear all;

load remigrado_0258_nl4.dat;
load vel_0258_RMS_new.dat;

q=remigrado_0258_nl4;

%q(isnan(q))=0.0;

R=vel_0258_RMS_new;
T=vel_0258_RMS_new;

[nt nx]=size(q);

%campos1(1:nt,1)=ones(nt,1);

%q=campos1;

dx=12.52;
dt=0.004;
nz=160;
dz=5.0;
ntrace=nx;

x0s=0.0; 
x0r=50.0;
h=(x0r-x0s)/2;

flagdip1=0.022;  % Vínculo de mergulho para migração
flagdip2=0.015;  % Vínculo de mergulho para submigração
flagdip3=0.022;  % Vínculo de mergulho para remigração

secout=q; %ormsby(nt,nx,dt,dx,3,6,20,40,q);

for ix=1:nx;
p1=secout(:,ix);
p1=halfdif2(p1,nt,dt,1);
secout(1:nt,ix)=p1(1:nt);
end;

dado=secout;

clear q;

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

vm=R(it,ix)/1;
v2=T(it,ix)/2;

if(vm~=0 && v2~=0); % LAÇO DA VELOCIDADE NULA

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

if(itime<=nt && ixtrace_co<=ntrace && ldip1>flagdip1 | ldip1<-flagdip1);
amigi(it,ix)=amigi(it,ix)+(1/sqrt(2*pi))*wds*dado(itime,ixtrace_co)*dx;
amigic(it,ix)=amigic(it,ix)+(1/(sqrt(2*pi)*dx))*dado(itime,ixtrace_co)*dx;
else;end;

if(itime2<=nt && ixtrace_co<=ntrace && ldip2>flagdip2 | ldip2<-flagdip2);
amigi2(it,ix)=amigi2(it,ix)+(1/sqrt(2*pi))*wds2*dado(itime2,ixtrace_co)*dx;
else;end;

if(itime<=nt && ixtrace_co<=ntrace && ldip1<flagdip1 | ldip1>-flagdip1);
migi(it,ix)=migi(it,ix)+(1/sqrt(2*pi))*wds*dado(itime,ixtrace_co)*dx;
migic(it,ix)=migic(it,ix)+(1/(sqrt(2*pi)*dx))*dado(itime,ixtrace_co)*dx;
else;end;

if(itime2<=nt && ixtrace_co<=ntrace && ldip2<flagdip2 | ldip2>-flagdip2);
migi2(it,ix)=migi2(it,ix)+(1/sqrt(2*pi))*wds2*dado(itime2,ixtrace_co)*dx;
else;end;

end;

else;migi3(it,ix)=0.0;end; % FIM DO LAÇO DA VELOCIDADE NULA

end;
end;

migi2_esp=zeros(size(migi));
migi2_esp(:,360)=migi2(:,360);

clear xs xr;

cronodata=real(migi2);

%secout=ormsby(nt,nx,dt,dx,3,6,12,24,real(migi2));

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

vm=R(it,ix)/2;
v2=T(it,ix)/4;

if(vm~=0 && v2~=0); % LAÇO DA VELOCIDADE NULA

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2; 
ixtrace_co=1+round(xtrace/dx);

factor=vm^2-v2^2;
extra=4*(h^2)*((1/v2^2)-(1/vm^2));
radicando1=tau^2/4+((x-xtrace)^2/factor)+extra;
radicando2=tau^2/4+((x-xtrace)^2/factor)+extra;

if(radicando1>0 && radicando2>0);
teta1=sqrt(radicando1);
teta2=sqrt(radicando2);
teta=teta1+teta2;
else;
teta1=real(sqrt(radicando1));
teta2=real(sqrt(radicando2));
teta=real(teta1+teta2)
end;

iteta=round(1+(teta/dt));

deno_velo=1-(((v2/vm)^2)*(tau/teta));
nume_velo=1;
novo=nume_velo/deno_velo;

w1=(tau/4/v2^(3/2))*sqrt(2/tau)*2;
w2=(1/(teta1+teta2))*sqrt((1/teta1)+(1/teta2));
wrm=w1*w2*(novo); %*(teta^(2.21));

ldipr=dipmig(1,ixtrace);

if(iteta<=nt && ixtrace_co<=nx && ldipr>flagdip3 | ldipr<-flagdip3);
amigi3(it,ix)=amigi3(it,ix)+(1/sqrt(2*pi))*wrm*migi2(iteta,ixtrace_co)*dx;
else;end;

if(iteta<=nt && ixtrace_co<=nx && ldipr<flagdip3 | ldipr>-flagdip3);
migi3(it,ix)=migi3(it,ix)+wrm*migi2(iteta,ixtrace_co)*dx;
else;end;

end;

else;migi3(it,ix)=0.0;end; % FIM DO LAÇO DA VELOCIDADE NULA

end;
end;

migi3=migi3-amigi3;

k=real(migi);
save migrado_campos.dat k -ascii;
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