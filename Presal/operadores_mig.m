disp('Este programa deve ser rodado após teste_terceiro_campos');

v1=1500.00/1; % Velocidade lâmina d'água

x0s=0.0; 
x0r=50.0;

tmig=zeros(1,nx);
dipmig=zeros(1,nx);

ix1=360; xmig1=ix1*dx;
it1=145; tmig1=it1*dt;

for ixtrace=1:ntrace;

xs=x0s+(ixtrace-1)*dx;
xr=x0r+(ixtrace-1)*dx;

t1=sqrt(tmig1^2/4+((xmig1-xs)^2/v1^2));
t2=sqrt(tmig1^2/4+((xmig1-xr)^2/v1^2));

tmig(1,ixtrace)=t1+t2;

dipmig(1,ixtrace)=(dx/v1^2)*(((xmig1-xs)/t1)+((xmig1-xr)/t2));

end;

figure;
plot([1:nx]*dx,tmig);
xlabel('Distância (m)');ylabel('Tempo (s)');
set(gca,'ydir','reverse');
yyaxis right;
plot([1:nx]*dx,-dipmig); axis tight;grid;
ylabel('dTempo (s)');


