close all; clear all;

perg=input('Afastamento: 1-Finito; 2-Nulo-->');

nx=200;
nz=200;

ntrace=nx;

dx=25;
dz=25;

nt=1000;
dt=0.004;

z0=2500;x0=0.0;
v1=2500;
v2=2000;
v10=2500;

x0s=0.0; 
x0r=50.0;
h=(x0r-x0s)/2;

g=ricker(100,10,0.004);
ng=length(g);

dado=zeros(nt,nx);
dado2=zeros(nt+ng-1,nx);
mvel=zeros(nz,nx);

%% ===== MODELAGEM === %%

if(perg==2); % AFASTAMENTO-NULO

for ix=1:nx;
x=x0+(ix-1)*dx;

for ix0=1:nx;

x00=0.0+(ix0-1)*dx;

if(x00>=1500 && x00<=3500);
prof=z0+(950+50)*(exp(-(x00-2500).^2/300000));
else;
prof=z0;
end;

iz=round(prof/dz);
mvel(1:iz,ix0)=v1;
mvel(iz+1:nz,ix0)=3300.;

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

elseif(perg==1); % AFASTAMENTO-FINITO

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

for ix0=1:nx
x00=0.0+(ix0-1)*dx;  

if(x00>=1500 && x00<=3500);
prof=z0+(950+50)*(exp(-(x00-2500).^2/300000));
else;
prof=z0;
end;

gx=1.0;

iz=round(prof/dz);
mvel(1:iz,ix0)=v1;
mvel(iz+1:nz,ix0)=3300.;

%v1=v10+0.03*x00; % Variação lateral de velocidade

v1=v10;

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

end; % FIM DAS ESCOLHA DO TIPO DE AFASTAMENTO

dado(1:nt,:)=dado2(1:nt,:);

%% ==== MIGRAÇÃO ==== %%

migi=zeros(nz,nx);
migi2=zeros(nz,nx);

if(perg==2); % AFASTAMENTO-NULO

for ix=1:nx; x=dx+(ix-1)*dx; 
for iz=1:nz; z=dz+(iz-1)*dz;

for ixtrace=1:nx;xt=dx+(ixtrace-1)*dx;

time=(2/v1)*sqrt(z^2+(x-xt)^2);
time2=(2/v2)*sqrt(z^2+(x-xt)^2);
itime=1+round(time/dt);
itime2=1+round(time2/dt);

l=sqrt(z^2+(x-xt)^2);
wds=2*z*sqrt(2/(v1*l));

if(itime<=nt && itime2<=nt);
migi(iz,ix)=migi(iz,ix)+wds*dado(itime,ixtrace)*dx;
migi2(iz,ix)=migi2(iz,ix)+wds*dado(itime2,ixtrace)*dx;
else;end;

end;

end;
end;

%% ==== REMIGRAÇÃO ==== %%

migi3=zeros(nz,nx);

for ix=1:nx;
p2=migi2(:,ix);
p2=halfdif2(p2,nz,dz,1);
migi2(1:nz,ix)=p2(1:nz);
end;

for ix=1:nx; x=dx+(ix-1)*dx; ix
for iz=1:nz; z=dz+(iz-1)*dz;

for ixtrace=1:nx;xz=dx+(ixtrace-1)*dx;

factor=1-(v2/v1)^2;
zeta=(v2/v1)*sqrt(z^2+((x-xz)^2/factor));
izeta=1+round(zeta/dz);

l1=sqrt(z^2+(x-xz)^2);
l2=sqrt(zeta^2+(x-xz)^2);

w1=2*z*sqrt(2/(v1*l1));
w2=(v2/(4*l2))*sqrt(2/(v2*l2));

w=w1*w2;

if(izeta<=nz);
migi3(iz,ix)=migi3(iz,ix)+w*migi2(izeta,ixtrace)*dx;
else;end;

end;

end;
end;

elseif(perg==1); % AFASTAMENTO-FINITO

migi3=zeros(nz,nx);

for ix=1:nx; x=dx+(ix-1)*dx; ix
for iz=1:nz; z=dz+(iz-1)*dz;

for ixtrace=1:ntrace;
    
xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

time_1=(1/v1)*sqrt(z^2+((x-xs)^2));
time_2=(1/v1)*sqrt(z^2+((x-xr)^2));
time=time_1+time_2;

time21=(1/v2)*sqrt(z^2+((x-xs)^2));
time22=(1/v2)*sqrt(z^2+((x-xr)^2));
timee=time21+time22;

itime=round(1 +(time/dt));
itime2=round(1+(timee/dt));

wds=(z/dx)*sqrt((1/time_1)+(1/time_2))*((time_1/time_2)+(time_2/time_1));
wds2=(z/dx)*sqrt((1/time21)+(1/time22))*((time21/time22)+(time22/time21));

if(itime<=nt && itime2<=nt && ixtrace_co<=ntrace);
migi(iz,ix)=migi(iz,ix)+(1/sqrt(2*pi))*wds*dado(itime,ixtrace_co)*dx;
migi2(iz,ix)=migi2(iz,ix)+(1/sqrt(2*pi))*wds2*dado(itime2,ixtrace_co)*dx;
else;end;

end; % Laço do traço

end; % Laço-z
end; % Laço-x

migi3(1:nz,1:nx)=0.0;
    
end;  % FIM DA ESCOLHA DO MÉTODO DE AFASTAMENTO  

kk=real(migi);kk1=kk';
save sin_2D_depth.dat kk1 -ascii;
clear kk, kk1;

kl=mvel';
save vel_2D_depth_transp.dat kl -ascii;
clear kl;

km=real(dado);
save dado_2D_tempo.dat km -ascii;
clear km;

figure;
imagesc([1:nx]*dx,[1:nt]*dt,-real(dado));colormap gray;

figure;
imagesc([1:nx]*dx,[1:nz]*dz,-real(migi));colormap gray;
figure;
imagesc([1:nx]*dx,[1:nz]*dz,-real(migi2));colormap gray;
figure;
imagesc([1:nx]*dx,[1:nz]*dz,real(migi3));colormap gray;

figure;
imagesc([1:nx]*dx,[1:nz]*dz,mvel),colorbar;

