disp('Este prorama deve ser rodado após segundo');

%% ==== REMIGRAÇÃO FORÇA BRUTA ==== %%

disp('Remigra os dados sísmicos com FORÇA BRUTA!');

migi4=zeros(nt,nx);

for ix=1:nx; x=dx+(ix-1)*dx; ix  
for it=1:nt; tau=dt+(it-1)*dt;

for ixtrace=1:nx;
  
xz=dx+(ixtrace-1)*dx;

tempo=sqrt(tau^2+(4*((x-xz)^2)/vm^2));
itempo=round(1+(tempo/dt));

limite=(v2/2)*tempo;
nl=round(limite/dx);

linf=ixtrace-round(nl/2);
lup=ixtrace+round(nl/2);

if(linf<1);linf=1;else;end;
if(lup>nx);lup=nx;else;end;

for ix2=linf:lup; % LAÇO MAIS INTERNO

x2=ix2*dx;

if((x2-xz)^2<limite^2);
teta2=(1/1)*tempo*sqrt(1-((x2-xz)^2/limite^2));
else;teta2=0.0;end;

iteta2=round(1+(teta2/dt));

deno_velo=1-(((v2/vm)^2)*(tau/teta2));
nume_velo=1;
novo=nume_velo/deno_velo;
wrm=(2/v2)*(tau/teta2)*(1/sqrt(tau*teta2))*(novo);

if(iteta2<=nt);
migi4(it,ix)=migi4(it,ix)+wrm*migi2(iteta2,ix2)*dx;
else;end;

end; % FIM DO LAÇO MAIS INTERNO

end;

end;
end;

figure;
imagesc([1:nx]*dx,[1:nt]*dt,real(migi4));colormap gray;