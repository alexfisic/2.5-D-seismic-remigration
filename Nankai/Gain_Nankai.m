disp('AGC nos dados da Fossa de Nankai');

load nankaii.dat;

d=nankaii';

[nt nx]=size(d);

dt=0.004;dx=16.67;

for ix=1:nx;
d(:,ix)=AGC(d(:,ix),8,nt);
end;

figure;
imagesc([1:nx]*dx,[1:nt]*dt,d(1300:1750,:)),colormap gray;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Dados sísmicos Nankai');

save nankai_gain.dat d -ascii;