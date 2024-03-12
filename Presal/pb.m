function secout=pb(nt,nx,dt,dx,fc,secin);

peso = zeros(1,nt);

wc = 2*pi*fc;
t = [1:floor(nt/8)]*(dt/2);

peso = (wc/pi)*(sin(wc*t)./(wc*t)); 

peso(1,:)=peso(1,:)/max(abs(peso(1,:)));

np=length(peso);

for ix=1:nx;
   secin(:,ix)=secin(:,ix)/max(abs(secin(:,ix)));
end;

secout = zeros(nt,nx);

for ix=1:nx;

eta = conv(peso,secin(:,ix));
secout (1:nt,ix) = eta(1:nt)';

end;

for ix=1:nx;
    secout(:,ix)=secout(:,ix)/max(abs(secout(:,ix))); 
end;

