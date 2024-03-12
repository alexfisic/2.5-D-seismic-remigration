disp('Este programa demigra da profundidade para obter seção sísmica');

close all; clear all;

load linha_0258_depth.dat;
load velocidade_0258_depth.dat;

[nz nx]=size(linha_0258_depth');

dx=12.52;
dz=5.;
dt=0.004;
ntrace=nx;
nt=1750;

saida=zeros(nt,ntrace);
R=velocidade_0258_depth';

x0s=0.; 
x0r=50.;
h=(x0r-x0s)/2;
t0=dt;
x0=x0s;

%% ==== DEMIGRAÇÃO ==== %%

disp('DeMigra os dados sísmicos!');

time_proc_1=cputime; % Para computar o tempo de maquina utilizado no processo;

dado_mig=zeros(nz,nx);
dado_mig=linha_0258_depth';

for it=950:nt; tau=dt+(it-1)*dt; it
for ixtrace=1:ntrace; 

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

for iz=1:nz; % Laço da profundidade

vm=R(iz,ixtrace);

limite=(vm/2)*tau;
nl=round(limite/dx);

linf=ixtrace_co-round(nl/4); % Para o modelo Campos, 2 ou 10
lup=ixtrace_co+round(nl/4);  % Para o modelo Campos, 2 ou 10

if(linf<1);linf=1;else;end;
if(lup>ntrace);lup=ntrace;else;end;

if((vm^2*tau^2)/4>h^2); % LAÇO DO "TEMPO"

for ix2=linf:lup; % LAÇO MAIS INTERNO

x2=ix2*dx;

if((x2-xtrace)^2<limite^2);
zeta=sqrt((vm^2*tau^2/4)-h^2)*sqrt(1-((x2-xtrace)^2/limite^2));
else;
zeta=0.0;
end;

izeta=round(1 +(zeta/dz));

aa=sqrt((vm^2*tau^2/4)-h^2);
bb=sqrt(1-((x2-xtrace)^2/limite^2));
estreti=tau*(bb/aa)+(aa/bb)*(8/(vm^2*tau^3))*(x2-xtrace)^2;

t1=(1/vm)*sqrt(zeta^2+((x2-xs)^2/1));
t2=(1/vm)*sqrt(zeta^2+((x2-xr)^2/1));

wis=(1/(2*vm))*(1/(t1+t2))*sqrt((1/t1)+(1/t2))*sqrt(estreti);

if(izeta<=nz && ixtrace_co<=ntrace);
saida(it,ixtrace_co)=saida(it,ixtrace_co)+wis*dado_mig(izeta,ix2)*dx;
else;end;

end; % FIM DO LAÇO MAIS INTERNO

else;saida(it,ixtrace_co)=1e-03;end; % FIM DO LAÇO DO "TEMPO"

end; % Fim do laço da profundidade

end;
end;

time_proc=cputime-time_proc_1;
duracao=time_proc/3600;   
disp('A demigração durou (horas):');
disp(duracao);
duracao2=time_proc/60;  
disp('A demigração durou (minutos):');
disp(duracao2);

gig=real(saida);
save remigrado.dat gig -ascii;
clear gig;

figure;
imagesc([1:ntrace]*dx,[1:nt]*dt,real(saida(1:nt,:)));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados sísmicos Campos demigrado');

figure;
imagesc([1:ntrace]*dx,[1:nz]*dz,real(dado_mig));colormap gray;grid;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title('Dados sísmicos Campos de entrada - profundidade');

