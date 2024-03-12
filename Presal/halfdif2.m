function [cx]=halfdif2(cx,nx,dt,sign);

%sign=1;
dw=2*pi/((nx)*dt);

y=fft(cx);

for i=1:round(nx/2);
y(i)=2.0*y(i)*sqrt(sign*(i-1)*dw/(2*pi));
y(i+round(nx/2))=0;
end;

cx=ifft(y);
