disp('Esse programa transforma as velocidades intervalares para RMS - tempo');

load velo1.dat;

[nt nx]=size(velo1);

dt=0.004;
dx=50.0;
ntime=nt;

modevel=velo1;

v_nmo=zeros(nt,nx);
vN=zeros(nt,nx);

velRMS=zeros(nt,nx);

perg=input('Transformção: 1-INTERVALAR; 2-CONTÍNUA-->');

if(perg==1); % Laço do método

for ix=1:nx; % Laço-x

   x=dx+(ix-1)*dx;

   t1=2.0;
   t2=3.0;
   t3=4.0;
   t4=5.0;

   nt1=round(t1/dt);
   nt2=round(t2/dt);
   nt3=round(t3/dt);
   nt4=round(t4/dt);

for it=2:nt1; % Inicio pseudocamada-1

   for i=1:it;
   if(i==1);
v_nmo(it,ix)=v_nmo(it,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(it,ix)=v_nmo(it,ix)+(i-(i-1))*dt*(modevel(i,ix)^2);
   end;
   end;

vN(it,ix)=sqrt(v_nmo(it,ix)/(it*dt));

itime1=it;

if(itime1<=nt);
velRMS(itime1,ix)=vN(it,ix);
else;
end;

velRMS(1,ix)=1500.0;

for it=2:itime1;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime1,ix)-velRMS(1,ix))*((it*dt-dt)/(itime1*dt-dt));
end;

end; % Fim pseudocamada-1

for it=nt1+1:nt2; % Inicio pseudocamada-2

   for i=1:it;
   if(i==1);
v_nmo(it,ix)=v_nmo(it,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(it,ix)=v_nmo(it,ix)+(i-(i-1))*dt*(modevel(i,ix)^2);
   end;
   end;

vN(it,ix)=sqrt(v_nmo(it,ix)/(it*dt));

itime2=it;

if(itime2<=nt);
velRMS(itime2,ix)=vN(it,ix);
else;
end;

velRMS(1,ix)=1500.0;

for it=itime1+1:itime2;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime2,ix)-velRMS(1,ix))*((it*dt-dt)/(itime2*dt-dt));
end;

end; % Fim pseudocamada-2

for it=nt2+1:nt3; % Inicio pseudocamada-3

   for i=1:it;
   if(i==1);
v_nmo(it,ix)=v_nmo(it,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(it,ix)=v_nmo(it,ix)+(i-(i-1))*dt*(modevel(i,ix)^2);
   end;
   end;

vN(it,ix)=sqrt(v_nmo(it,ix)/(it*dt));

itime3=it;

if(itime3<=nt);
velRMS(itime3,ix)=vN(it,ix);
else;
end;

velRMS(1,ix)=1500.0;

for it=itime2+1:itime3;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime3,ix)-velRMS(1,ix))*((it*dt-dt)/(itime3*dt-dt));
end;

end; % Fim pseudocamada-3

for it=nt3+1:nt4; % Inicio pseudocamada-4

   for i=1:it;
   if(i==1);
v_nmo(it,ix)=v_nmo(it,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(it,ix)=v_nmo(it,ix)+(i-(i-1))*dt*(modevel(i,ix)^2);
   end;
   end;

vN(it,ix)=sqrt(v_nmo(it,ix)/(it*dt));

itime4=it;

if(itime4<=nt);
velRMS(itime4,ix)=vN(it,ix);
else;
end;

velRMS(1,ix)=1500.0;

for it=itime3+1:itime4;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime4,ix)-velRMS(1,ix))*((it*dt-dt)/(itime4*dt-dt));
end;

end; % Fim pseudocamada-4

for it=nt4+1:nt; % inicio pseudocamada-5

for i=1:it;
   if(i==1);
v_nmo(it,ix)=v_nmo(it,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(it,ix)=v_nmo(it,ix)+(i-(i-1))*dt*(modevel(i,ix)^2);
   end;
end;

vN(it,ix)=sqrt(v_nmo(it,ix)/(i*dt));

itime5=it;

if(itime5<=nt);
velRMS(itime5,ix)=vN(it,ix);
else;
end;

velRMS(1,ix)=1500.0;

for it=itime4+1:ntime;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime5,ix)-velRMS(1,ix))*((it*dt-1*dt)/(itime5*dt-1*dt));
end;

end; %Fim pseudocamada-5

end; % Fim do laço-x

elseif(perg==2); % Laço do método

for ix=1:nx; % Laço-x

for it=2:nt; % Inicio pseudocamada

   for i=1:it;
   if(i==1);
v_nmo(it,ix)=v_nmo(it,ix)+(modevel(i,ix)^2)*dt;
   else;
v_nmo(it,ix)=v_nmo(it,ix)+(i-(i-1))*dt*(modevel(i,ix)^2);
   end;
   end;

vN(it,ix)=sqrt(v_nmo(it,ix)/(it*dt));

itime1=it;

if(itime1<=nt);
velRMS(itime1,ix)=vN(it,ix);
else;
end;

velRMS(1,ix)=modevel(1,ix);

for it=2:itime1;
   velRMS(it,ix)=velRMS(1,ix)+(velRMS(itime1,ix)-velRMS(1,ix))*((it*dt-dt)/(itime1*dt-dt));
end;

end; % Fim pseudocamada

end; % Fim do laço-x

end; % Fim do IF do método

save vel_campos1_RMS.dat velRMS -ascii;

figure;
subplot(1,2,1);
imagesc([1:nx]*dx,[1:nt]*dt,velo1);colorbar;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Modelo de velocidades INT - tempo (m/s)');
subplot(1,2,2);
imagesc([1:nx]*dx,[1:nt]*dt,velRMS);colorbar;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Modelo de velocidades RMS - tempo (m/s)');



