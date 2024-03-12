

load vel_rms.dat;

vin=reshape(vel_rms,460,150);

[ntime nx] = size(vin);

dt=0.004; dx=25;
x0=0;

VI = zeros(ntime,nx);
VS = zeros(size(VI));

VI=vin;

ap=1;

for iz=1:ntime;
    for ix=1:nx;
    x=x0+(ix-1)*dx;
    VS(iz,ix) = VI(iz,1)+(VI(iz,nx)-VI(iz,1))*((x-x0)/(nx*dx-x0));
    end;
end;


for ix=1:nx;
    for iz=1:ntime;
    t=dt+(iz-1)*dt;
    VS(iz,ix) = VI(1,ix)+(VI(ntime,ix)-VI(1,ix))*((t-0.0)/(ntime*dt-0.0));
    end;
end;


figure,imagesc(VS),colorbar;

figure,imagesc(vin),colorbar;


