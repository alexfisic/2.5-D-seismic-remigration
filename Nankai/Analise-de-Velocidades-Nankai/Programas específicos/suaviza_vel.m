

load vel_rms.dat;

vin=reshape(vel_rms,1590,150);

[ntime nx] = size(vin);

VI = zeros(ntime,nx);
VS = zeros(size(VI));

VI=vin;

ap=5;

for iz=1:ntime
    for ix=1:nx
        istart = 1 + (ix-1) - ap;
        iend = 1 + (ix-1) + ap;
        if istart<1;istart=1;end;
        if iend>nx;iend=nx;end;
        VS(iz,ix) = mean(VI(iz,istart:iend));
    end
end


for iz=1:ntime
    istart = 1 + (iz-1) - ap;
    iend = 1 + (iz-1) + ap;
    if istart<1;istart=1;end;
    if iend>ntime;iend=ntime;end;
    for ix=1:nx
        VS(iz,ix) = mean(VS(istart:iend,ix));
    end
end

VMIG=2*VS;

figure,imagesc(2*VS),colorbar;

figure,imagesc(vin),colorbar;

save vel_migrat.dat VMIG -ascii;


