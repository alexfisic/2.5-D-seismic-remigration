function secout=ormsby(nt,nx,dt,dx,f1,f2,f3,f4,secin);

sect=secin;

peso = zeros(1,nt);

t = [1:floor(nt/8)]*(dt/2);

w1 = 2*pi*f1; w2 = 2*pi*f2; w3 = 2*pi*f3; w4 = 2*pi*f4;

p1 = (cos(w2*t)./(t.^2)); p2 = (cos(w1*t)./(t.^2));

p3 = (cos(w4*t)./(t.^2)); p4 = (cos(w3*t)./(t.^2));

peso = (1/(pi*(w2-w1)))*(p1 - p2) - (1/(pi*(w4-w3)))*(p3 - p4);

peso(1,:)=peso(1,:)/max(abs(peso(1,:)));

np=length(peso);

%for ix=1:nx;
%   sect(:,ix)=sect(:,ix)/max(abs(sect(:,ix)));
%end;

secout = zeros(nt,nx);

for ix=1:nx;

eta = conv(peso,sect(:,ix));
secout (1:nt,ix) = eta(1:nt)';

end;

%for ix=1:nx;
%    secout(:,ix)=secout(:,ix)/max(abs(secout(:,ix))); 
%end;
