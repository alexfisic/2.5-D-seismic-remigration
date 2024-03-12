disp('Este programa deve ser rodado após segundo ou terceiro');

[nt nx]=size(migi);

dx=25.0;
dt=0.008;
ntrace=300; % Se for escolhido terceiro, ntrace=300

saida=zeros(nt,ntrace);

x0s=0.0; 
x0r=50.0;
h=(x0r-x0s)/2;

%% ==== DEMIGRAÇÃO ==== %%

disp('DeMigra os dados sísmicos!');

dado_mig=zeros(nt,nx);
dado_mig=migi;

for ixtrace=1:ntrace; ixtrace

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

for it=1:nt; tau=dt+(it-1)*dt;
%for ix=1:nx; x=dx+(ix-1)*dx; 
   
limite=(vm/2)*tau;
nl=round(limite/dx);

linf=ixtrace_co-round(nl/2); % Para o modelo Campos, 10
lup=ixtrace_co+round(nl/2);  % Para o modelo Campos, 10

if(linf<1);linf=1;else;end;
if(lup>ntrace);lup=ntrace;else;end;

if(tau^2>4*h^2/vm^2); % LAÇO DO TEMPO

for ix2=linf:lup; % LAÇO MAIS INTERNO

x2=ix2*dx;

if((x2-xtrace)^2<limite^2);
time=sqrt(tau^2-(4*h^2/vm^2))*sqrt(1-((x2-xtrace)^2/limite^2));
else;
time=real(sqrt(tau^2-(4*h^2/vm^2))*sqrt(1-((x2-xtrace)^2/limite^2)));
end;

itime=round(1 +(time/dt));

aa=sqrt(tau^2-(4*h^2/vm^2));
bb=sqrt(1-((x2-xtrace)^2/limite^2));
estreti=tau*(bb/aa)+(aa/bb)*(8/(vm^2*tau^3))*(x2-xtrace)^2;

t1=sqrt(time^2/4+((x2-xs)^2/vm^2));
t2=sqrt(time^2/4+((x2-xr)^2/vm^2));

wis=(1/(2*vm))*(1/(t1+t2))*sqrt((1/t1)+(1/t2))*sqrt(estreti);

if(itime<=nt && ixtrace_co<=ntrace);
saida(it,ixtrace_co)=saida(it,ixtrace_co)+wis*dado_mig(itime,ix2)*dx;
else;end;

end; % FIM DO LAÇO MAIS INTERNO

else;saida(it,ixtrace_co)=1e-03;end; % FIM DO LAÇO DO TEMPO

%end;

end;
end;

figure;
imagesc([1:ntrace]*dx,[1:nt]*dt,-real(saida(1:nt,:)));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados sísmicos demigrado');

figure;
imagesc([1:ntrace]*dx,[1:nt]*dt,real(dado_mig));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados migrado de entrada');
