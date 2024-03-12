% == Transformada de ondaleta == %

x1=real(campos1(:,500));

tg=[1:nt]*dt; 
dfw=1/(nt*dt);
nciclos=5; nf=nt;
ii=complex(0,1);
WT=zeros(nt,nf);

for itau=1:nt; tau = dt + (itau-1)*dt;

for ifi=1:nf; efe=dfw+(ifi-1)*dfw;

sigma=nciclos/(2*pi*efe); 

wmorlet = exp(-ii*2*pi*efe*((tg-tau)/1)).*exp(-(0.5*((1/sigma)^2)*(tg-tau).^2));

passante2 = wmorlet.*x1';
passante3=sum(passante2)*dt/sqrt(sigma);

WT(itau,ifi) = abs(passante3);

end;

end;

figure;
subplot(1,2,1);
plot(x1',[1:nt]*dt); axis tight;grid;
set(gca,'YDir','reverse');
subplot(1,2,2);
imagesc([1:round(nf/5)]*dfw,[1:nt]*dt,WT(1:nt,1:round(nf/5)));grid;
