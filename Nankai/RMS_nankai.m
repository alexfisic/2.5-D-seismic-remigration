disp('Este programa deve ser rodado após TTI_NANKAI');

v_nmo=zeros(nz,nx);
vN=zeros(nz,nx);

nt=2750;dt=0.004;
ntime=nt;

velRMS=zeros(nt,nx);

for ix=1:nx; % Laço-x

   x=dx+(ix-1)*dx;

   z1=4100.;
   z2=4250.;
   z3=6067.;
   z4=8335.;

   velRMS(1,ix)=modevel(1,ix);

   nz1=round(z1/dz);
   nz2=round(z2/dz);
   nz3=round(z3/dz);
   nz4=round(z4/dz);

for iz=2:nz1; % Inicio camada-1

   for i=1:iz;
   if(i==1);
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)^2);
   end;
   end;

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix));

itime1=round(tti(iz,ix)/dt);

if(itime1<=nt);
velRMS(itime1,ix)=vN(iz,ix);
else;
end;

for it=2:itime1;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime1,ix)-velRMS(1,ix))*((it*dt-dt)/(itime1*dt-dt));
end;

end; % Fim camada-1

for iz=nz1+1:nz2; % Inicio camada-2

   for i=1:iz;
   if(i==1);
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)^2);
   end;
   end;

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix));

itime2=round(tti(iz,ix)/dt);

if(itime2<=nt);
velRMS(itime2,ix)=vN(iz,ix);
else;
end;

for it=itime1+1:itime2;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime2,ix)-velRMS(1,ix))*((it*dt-dt)/(itime2*dt-dt));
end;

end; % Fim camada-2

for iz=nz2+1:nz3; % Inicio camada-3

   for i=1:iz;
   if(i==1);
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)^2);
   end;
   end;

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix));

itime3=round(tti(iz,ix)/dt);

if(itime3<=nt);
velRMS(itime3,ix)=vN(iz,ix);
else;
end;

for it=itime2+1:itime3;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime3,ix)-velRMS(1,ix))*((it*dt-dt)/(itime3*dt-dt));
end;

end; % Fim camada-3

for iz=nz3+1:nz4; % Inicio camada-4

   for i=1:iz;
   if(i==1);
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)^2);
   end;
   end;

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix));

itime4=round(tti(iz,ix)/dt);

if(itime4<=nt);
velRMS(itime4,ix)=vN(iz,ix);
else;
end;

for it=itime3+1:itime4;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime4,ix)-velRMS(1,ix))*((it*dt-dt)/(itime4*dt-dt));
end;

end; % Fim camada-4

for iz=nz4+1:nz; % inicio camada-5

for i=1:iz;
   if(i==1);
v_nmo(iz,ix)=v_nmo(iz,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(iz,ix)=v_nmo(iz,ix)+(tti(i,ix)-tti(i-1,ix))*(modevel(i,ix)^2);
   end;
end;

vN(iz,ix)=sqrt(v_nmo(iz,ix)/tti(iz,ix));

itime5=round(tti(iz,ix)/dt);

if(itime5<=nt);
velRMS(itime5,ix)=vN(iz,ix);
else;
end;

for it=itime4+1:ntime;
   velRMS(it,ix)=velRMS(itime4,ix)+(velRMS(itime5,ix)-velRMS(itime4,ix))*((it*dt-itime4*dt)/(itime5*dt-itime4*dt));
end;

end; %Fim camada-5

end; % Fim do laço-x

save vel_nankai_RMS.dat velRMS -ascii;

figure;
imagesc([1:nx]*dx,[1:nt]*dt,velRMS);colorbar;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Modelo de velocidades - tempo (m/s)');
