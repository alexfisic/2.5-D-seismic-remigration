% == Transformada de Gabor == %

tg=[1:nt]*dt;
alpha=10.1;
G=zeros(nt,nf);

for itau=1:nt; tau = dt + (itau-1)*dt;

peso = exp(-(alpha*pi*(tg-tau).^2));

passante = fft(peso.*x1');

a = real(passante);
b = imag(passante);

G(itau,:) = sqrt(a.^2 + b.^2)/nt;

end;

figure;
subplot(1,2,1);
plot(x1',[1:nt]*dt); axis tight;grid;
set(gca,'YDir','reverse');
subplot(1,2,2);
imagesc([1:round(nf/5)]*dfw,[1:nt]*dt,G(1:nt,1:round(nf/5)));grid;
