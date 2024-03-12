disp('Este programa deve ser rodado após RMS_nankai');

load vel_campos1_RMS.dat;

velo_ml = vel_campos1_RMS;

velo_out = zeros(size(velo_ml));

[ntime nx] = size(velo_ml);

dz = 0.004;
z = [1:ntime]*dz;
index = length(z);
T = z(index); %T=2*T;

n = input('Numero de termo da Serie de Fourier (n):');
c = zeros(1,n+1);
d = zeros(1,n+1);

for i=1:nx; i

v = velo_ml(:,i);
%v = v'; % Sempre verificar quando é necessário transpor!
v7 = zeros(1,ntime);

omega=(2*pi/T);

c(1) = a_fourier(0,T,v,dz,dz);
d(1) = b_fourier(1,T,v,dz,dz);

for j=3:n+1;
    c(j-1) = a_fourier(j-1,T,v,dz,dz);
    d(j-1) = b_fourier(j-1,T,v,dz,dz);
end;

for ik=1:ntime;
    for iy=3:n+1;
      v7(1,ik) = v7(1,ik) + c(iy-1)*cos((iy-1)*omega*z(ik)) + d(iy-1)*sin((iy-1)*omega*z(ik));   
    end;
      v7(1,ik) = v7(1,ik) + d(1)*sin(omega*z(ik));
end;

velo_out(1:ntime,i) = v7(1,1:ntime)+c(1)/2;

c=c*0; d=d*0;

end;

figure,imagesc([1:nx]*dx,[1:nt]*dt,velo_out); colorbar;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Modelo de velicidades suavizado (Fourier) - tempo (m/s)');


save vel_campos1_smooth.dat velo_out -ascii;