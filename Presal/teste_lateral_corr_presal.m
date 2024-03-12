disp('Este prorama deve ser rodado após teste_terceiro_presal');

%% ==== REMIGRAÇÃO ==== %%

disp('Remigra os dados sísmicos!');
disp('Esse programa testa a correção de variação lateral de velocidade');

migi6=zeros(nt,nx);
migi_aux=zeros(size(migi6));
migi_esp=zeros(size(migi_aux));

for ix=1:nx;
p2=migi(:,ix);
p2=halfdif2(p2,nt,dt,1);
migi_aux(1:nt,ix)=p2(1:nt);
end;

migi_esp(:,360)=migi_aux(:,360);

for ix=1:nx; x=dx+(ix-1)*dx; ix 
for it=1:nt; tau=dt+(it-1)*dt;

vm=T(it,ix);    % Inversão dos campos de velocidades (sem variação)
v2=R(it,ix);    % Inversão dos campos de velocidades (com variação)

if(vm~=0 && v2~=0); % LAÇO DA VELOCIDADE NULA

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;

xtrace=(xs+xr)/2; 
ixtrace_co=1+round(xtrace/dx);

factor1=(v2^2)*((tau^2/4)+(h^2/vm^2));
factor2=(tau^2)+(4*h^2)*((1/vm^2)-(1/v2^2));

if(factor1>(x-xtrace)^2 | factor2>0); 
teta1=real(sqrt(factor2)*sqrt(1-((x-xtrace)^2/factor1)));
teta2=real(sqrt(factor2)*sqrt(1-((x-xtrace)^2/factor1)));
teta=real(teta1);
end;

iteta=round(1+(teta/dt));

deno_velo=1-(((vm/v2)^2)*(tau/teta));
nume_velo=1;
novo=nume_velo/deno_velo;

w1=(tau/4/v2^(3/2))*sqrt(2/tau)*2;
w2=(1/(teta1+teta2))*sqrt((1/teta1)+(1/teta2));
wrm=w1*w2*(novo); %*(teta^(2.21));

if(iteta<=nt && ixtrace_co<=nx);
migi6(it,ix)=migi6(it,ix)+migi_aux(iteta,ixtrace_co)*dx;
else;end;

end;

else;migi6(it,ix)=0.0;end; % FIM DO LAÇO DA VELOCIDADE NULA

end;
end;

for ix=1:nx;
    migi6(:,ix)=migi6(:,ix)/max(abs(migi6(:,ix)));
end;

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[1:nt]*dt,real(migi6(1:nt,:)));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal corrigido - variação lateral');
subplot(1,2,2);
imagesc([1:nx]*dx,[1:nt]*dt,real(migi(1:nt,:)));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal migrado - original');

