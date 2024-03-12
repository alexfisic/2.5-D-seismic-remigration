close all; clear all;

%load urucu_velo_3.dat;
load SEG_vel_f.dat;

%velo_ml = urucu_velo_3';
velo_ml = SEG_vel_f;

velo_out = zeros(size(velo_ml));

[nz nx] = size(velo_ml);

%z = [10:10:3500];
dz = 10;
z = [1:364]*dz;
index = length(z);
T = z(index); %T=2*T;

n = input('Numero de termo da Serie de Fourier (n):');
c = zeros(1,n+1);
d = zeros(1,n+1);

for i=1:nx;

v = velo_ml(:,i)*1000;
v = v';
v7 = zeros(1,nz);

omega=(1*pi/T);
c(1) = a_fourier(0,T,v,dz,dz);
d(1) = a_fourier(1,T,v,dz,dz);

for j=3:n+1;
    c(j-1) = a_fourier(j-1,T,v,dz,dz);
    d(j-1) = b_fourier(j-1,T,v,dz,dz);
end;

for ik=1:nz;
    for iy=2:n+1;
      v7(1,ik) = v7(1,ik) + c(iy-1)*cos((iy-1)*omega*z(ik)) + d(iy-1)*sin((iy-1)*omega*z(ik));   
    end;
end;


velo_out(1:nz,i) = v7(1,1:nz);

end;

velo_out = velo_out/1000;

figure,imagesc(velo_ml);

figure,imagesc(velo_out);