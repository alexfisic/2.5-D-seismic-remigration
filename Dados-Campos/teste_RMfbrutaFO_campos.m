disp('Este prorama deve ser rodado após teste_terceiro_campos');

%% ==== REMIGRAÇÃO FORÇA BRUTA ==== %%

disp('Remigra os dados sísmicos com FORÇA BRUTA!');

time_proc_1=cputime; % Para computar o tempo de maquina utilizado no processo;

migi4=zeros(nt,nx);
amigi4=zeros(nt,nx);

for it=1:nt; tau=dt+(it-1)*dt; it
for ix=1:nx; x=dx+(ix-1)*dx;  

vm=R(it,ix)/2;
v2=T(it,ix)/4;

if(vm~=0 && v2~=0); % LAÇO DA VELOCIDADE NULA

for ixtrace=1:ntrace;
  
xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;
xtrace=(xs+xr)/2;

ixtrace_co=1+round(xtrace/dx);

tempo1=sqrt(tau^2/4+(((x-xs)^2)/vm^2));
tempo2=sqrt(tau^2/4+(((x-xr)^2)/vm^2));
tempo=tempo1+tempo2;

ldipm=(dx/vm^2)*(((x-xs)/tempo1)+((x-xr)/tempo2));

itempo=round(1+(tempo/dt));

% === Determinante de Beylkin (tempo) === %

b11=(1/vm^2)*(((x-xs)/tempo1)+((x-xr)/tempo2));
b12=(tau/2/vm)*((1/tempo1)+(1/tempo2));

b21=(-1/vm^2)*((1/tempo1)+(1/tempo2))+(1/vm^4)*(((x-xs)^2/tempo1^3)+...
    ((x-xr)^2/tempo2^3));
b22=(tau/2/vm^3)*(((x-xs)/tempo1^3)+((x-xr)/tempo2^3));

beylkin=b11*b22-b12*b21;

% ======================================= %

limite=(v2/2)*tempo;
nl=round(limite/dx);

linf=ixtrace_co-round(nl/10); % Melhor resultado: nl/2, nl/10
lup=ixtrace_co+round(nl/10);  % Melhor resultado: nl/2, nl/10

if(linf<1);linf=1;else;end;
if(lup>ntrace);lup=ntrace;else;end;

if(tempo^2>4*h^2/v2^2); % LAÇO DO TEMPO

for ix2=linf:lup; % LAÇO MAIS INTERNO

x2=ix2*dx;

if((x2-xtrace)^2<limite^2);
teta2=(1/1)*sqrt(tempo^2-(4*h^2/v2^2))*sqrt(1-((x2-xtrace)^2/limite^2));
else;
teta2=real(sqrt(tempo^2-(4*h^2/v2^2))*sqrt(1-((x2-xtrace)^2/limite^2)));
end;

iteta2=round(1+(teta2/dt));

aa=sqrt(tempo^2-(4*h^2/v2^2));
bb=sqrt(1-((x2-xtrace)^2/limite^2));
estreti=tempo*(bb/aa)+(aa/bb)*(8/(v2^2*tempo^3))*(x2-xtrace)^2;

ldip4=((-aa/bb)*4*(x2-xtrace)/(v2^2*tempo^2))*dx;

deno_velo=1-(((v2/vm)^2)*(tau/teta2));
nume_velo=1;
novo=nume_velo/deno_velo;

t1=sqrt(teta2^2/4+((x2-xs)^2/v2^2)); % TVT do ramo 1 da isócrona
t2=sqrt(teta2^2/4+((x2-xr)^2/v2^2)); % TVT do ramo 2 da isócrona

w1=(tau/2)*sqrt((1/tempo1)+(1/tempo2))*((tempo1/tempo2)+(tempo2/tempo1));
w2=(1/v2)*(1/(t1+t2))*sqrt((1/t1)+(1/t2));

wrm=(sqrt(2)/v2^(1/2)/dx/sqrt(2*pi))*w1*w2*(novo)*sqrt(estreti); %*(teta2^(2.21));
w_beylkin=((v2^(3/2))/4/1)*sqrt(estreti)*novo*beylkin;           %*(teta2^(2.21));

if(iteta2<=nt && ldipm>0.022 | ldipm<-0.022);
amigi4(it,ix)=amigi4(it,ix)+wrm*migi2(iteta2,ix2)*dx;
else;end;

if(iteta2<=nt && ldipm<0.022 | ldipm>-0.022);
migi4(it,ix)=migi4(it,ix)+wrm*migi2(iteta2,ix2)*dx;
else;end;

end; % FIM DO LAÇO MAIS INTERNO

else;migi4(it,ix)=0.0;end; % FIM DO LAÇO DO TEMPO

end;

else;migi4(it,ix)=0.0;end; % FIM DO LAÇO DA VELOCIDADE NULA

end;
end;

migi4=migi4-amigi4;

outi=real(migi4);
save campos_remigfb.dat outi -ascii;
clear outi;

time_proc=cputime-time_proc_1;
duracao=time_proc/3600;   
disp('A remigração força-bruta durou (horas):');
disp(duracao);
duracao2=time_proc/60;  
disp('A remigração força-bruta durou (minutos):');
disp(duracao2);

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi4(1:nt,:)));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Campos remigrado - força bruta');