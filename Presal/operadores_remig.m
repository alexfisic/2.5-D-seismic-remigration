disp('Este programa deve ser rodado após teste_terceiro_campos');

v1=1500.00; % Velocidade lâmina d'água
v2=v1/2;    % Velocidade subamostrada

x0s=0.0; 
x0r=50.0;

tmig=zeros(1,nx);
tremig=zeros(1,nx);
dipmig=zeros(1,nx);
dipremig=zeros(1,nx);

ix1=360; xmig1=ix1*dx;
it1=145; tmig1=it1*dt;

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;

xtrace=(xs+xr)/2;
ixtrace_co=1+round(xtrace/dx);

t1=sqrt(tmig1^2/4+((xmig1-xs)^2/v1^2));
t2=sqrt(tmig1^2/4+((xmig1-xr)^2/v1^2));

tmig(1,ixtrace)=t1+t2;

dipmig(1,ixtrace)=(dx/v1^2)*(((xmig1-xs)/t1)+((xmig1-xr)/t2));

ldipd=dipmig(1,ixtrace);

if(ixtrace == 380); % Condição de mergulho

limite=(v2/2)*(t1+t2);
nl=round(limite/dx);

linf=ixtrace_co-round(nl); % Melhor resultado: nl/2, nl/10
lup=ixtrace_co+round(nl);  % Melhor resultado: nl/2, nl/10

if(linf<1);linf=1;else;end;
if(lup>ntrace);lup=ntrace;else;end;

if((t1+t2)^2>4*h^2/v2^2); % LAÇO DO TEMPO

for ix2=linf:lup; % LAÇO MAIS INTERNO

x2=ix2*dx;

if((x2-xtrace)^2<limite^2);
teta=sqrt((t1+t2)^2-(4*h^2/v2^2))*sqrt(1-((x2-xtrace)^2/limite^2));
else;
teta=real(sqrt((t1+t2)^2-(4*h^2/v2^2))*sqrt(1-((x2-xtrace)^2/limite^2)));
end;

aa=sqrt((t1+t2)^2-(4*h^2/v2^2));
bb=sqrt(1-((x2-xtrace)^2/limite^2));

dipremig(1,ix2)=((-aa/bb)*4*(x2-xtrace)/(v2^2*(t1+t2)^2))*dx;

tremig(1,ix2)=teta;

end; % FIM DO LAÇO MAIS INTERNO

else;end; % FIM DO LAÇO DO TEMPO

else;end; % Fim da condição de mergulho

end;

figure;
plot([1:nx]*dx,tremig,[1:nx]*dx,tmig);
xlabel('Distância (m)');ylabel('Tempo (s)');
set(gca,'ydir','reverse');
yyaxis right;
plot([1:nx]*dx,real(dipremig),[1:nx]*dx,-real(dipmig)); axis tight;grid;
ylabel('dTempo (s)');

