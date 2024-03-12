close all; clear all;

disp('Razão de escala temporal do dado sísmico do presal');
 
load campos3.dat;

sect = real(campos3);
[nt nx]=size(sect);

campos3(1:nt,1)=ones(nt,1);
sect=campos3;

dt = 0.004;
dx = 50;
df = 1/(nt*dt);

volume_H = zeros(nt,round(nt/5),nx);

for ix=1:nx;
   sect(:,ix)=sect(:,ix)/max(abs(sect(:,ix)));
end;

data = zeros(nt,nx);

spec_mean = zeros(nt,round(nt/5));

tg=[1:nt]*dt; 
dfw=1/(nt*dt);
nciclos=5; nf=nt;
ii=complex(0,1);

disp('Cria o volume de dados');

for ix=1:20; ix %nx; ix

x1=sect(:,ix);

% == Transformada de Ondaleta == %

for itau=1:nt; tau = dt + (itau-1)*dt; 

for ifi=1:round(nf/5); efe=dfw+(ifi-1)*dfw;

sigma=nciclos/(2*pi*efe); 

wmorlet = exp(-ii*2*pi*efe*((tg-tau)/1)).*exp(-(0.5*((1/sigma)^2)*(tg-tau).^2));

passante2 = wmorlet.*x1';
passante3=sum(passante2)*dt/sqrt(sigma);

volume_H(itau,ifi,ix) = abs(passante3);

end; % Fim laço frquência

end; % Fim laço-tau

end; % Fim laço-x

% == Calcula o espectro médio == %

disp('Calcula o espectro médio');

for it=1:nt;

for ifi=1:round(nf/5);

for ix=2:nx;
 spec_mean(it,ifi) = spec_mean(it,ifi) + volume_H(it,ifi,ix)*dx;
end;
spec_mean(it,ifi) = spec_mean(it,ifi)/((nx-1)*dx);
end;

end;

save espectro.dat spec_mean -ascii;

disp('Calcula razão de escala');

rescala=zeros(nt,1);

for it=1:nt;
    for ifi=1:round(nt/5);
        rescala(it,1) = rescala(it,1) + spec_mean(it,ifi)*dfw;
    end;
end;

rescala=rescala/((round(nt/5)-1)*dfw);

save razao.dat rescala -ascii;

figure;
subplot(1,2,1);
plot(rescala',[1:nt]*dt); axis tight;grid;
set(gca,'YDir','reverse');
subplot(1,2,2);
imagesc([1:round(nf/5)]*dfw,[1:nt]*dt,spec_mean(1:nt,1:round(nf/5)));grid;

