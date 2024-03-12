clear SEG_novo;
load SEG_vel_f.dat;

[nz nx]=size(SEG_vel_f);

dx=input('Amostragem anterior (dx)-->');
dxn=input('Nova amostragem (dxn)-->');
nxn=input('Nova amostragem em x (nxn)-->');
nzn=input('Nova amostragem em z (nzn)-->');

SEG_novo=zeros(nzn,nxn);

dz=dx;
dzn=dxn;

disp('Reamostragem!!');

for i=1:nz;
 z=dz+(i-1)*dz;
 iz=round(z/dzn);
for j=1:nx;
 x=dx+(j-1)*dx;
 ix=round(x/dxn);
SEG_novo(iz,ix)=SEG_vel_f(i,j);
end;
end;

disp('Bordas');

SEG_novo(1,:)=SEG_novo(2,:);
SEG_novo(nzn,:)=SEG_novo(nzn-2,:);
SEG_novo(:,1)=SEG_novo(:,2);
SEG_novo(:,nxn)=SEG_novo(:,nxn-2);

disp('Interpolação das colunas!!');

nznn=(nzn-1)/2;
nxnn=(nxn/2);

for i=1:nznn;
iz=2*i+1;
for j=1:nxn;

if(iz==nzn);iz=iz-2;
%===========================================================================%
if(SEG_novo(iz,j)==0.0 & SEG_novo(iz+1,j)~= 0.0);

SEG_novo(iz,j)=0.5*(SEG_novo(iz+1,j)+SEG_novo(iz-1,j));

elseif(SEG_novo(iz,j)==0.0 & SEG_novo(iz+1,j)==0.0 & SEG_novo(iz-1,j)~=0.0);

SEG_novo(iz,j)=0.5*(SEG_novo(iz+2,j)+SEG_novo(iz-1,j));

else;
end;
%===========================================================================%
else;
end;

end;
end;

disp('Figura');

figure,imagesc(SEG_novo);colormap(ocean);
