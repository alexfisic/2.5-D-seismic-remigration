close all; clear all;

disp('Calcula os tempos TTI do modelo de velocidade Nankai');

load nankai_dvel.dat;

[nz nx]=size(nankai_dvel);

tti=zeros(nz,nx);
modevel=zeros(nz,nx);

modevel=nankai_dvel;

dx=16.67;
dz=dx;

for ix=1:nx; ix % Laço-x

   x=dx+(ix-1)*dx;
   t0=dz/modevel(1,ix);
   tti(1,ix)=t0;

   z1=4100.; 
   z2=4250.;
   z3=6067.;
   z4=8335.;

   nz1=round(z1/dz);
   nz2=round(z2/dz);
   nz3=round(z3/dz);
   nz4=round(z4/dz);

for iz=2:nz1;
   for i=1:iz;
tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix));
   end;
tti(iz,ix)=(tti(iz,ix)+t0)*2.0;
end;

for iz=nz1+1:nz2;
   for i=nz1:iz;
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0;
   end;
   tti(iz,ix)=(tti(iz,ix)+tti(nz1,ix))*1.0;
end;

for iz=nz2+1:nz3;
   for i=nz2:iz;
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0;
   end;
   tti(iz,ix)=(tti(iz,ix)+tti(nz2,ix));
end;

for iz=nz3+1:nz4;
   for i=nz3:iz;
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0;
   end;
   tti(iz,ix)=(tti(iz,ix)+tti(nz3,ix));
end;

for iz=nz4+1:nz;
   for i=nz4:iz;
   tti(iz,ix)=tti(iz,ix)+(dz/modevel(i,ix))*2.0;
   end;
   tti(iz,ix)=(tti(iz,ix)+tti(nz4,ix));
end;

end; % Fim do laço-x

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[1:nz]*dz,nankai_dvel);colorbar;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title('Modelo de velicidades - profundidade (m/s)');
subplot(1,2,2);
imagesc([1:nx]*dx,[1:nz]*dz,tti);colorbar;
xlabel('Distância (m)');ylabel('Profundidade (m)');
title('Tempo de trânsito intervalar - TTI (s)');
