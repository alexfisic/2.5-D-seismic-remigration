disp('Tenta equalizar as amplitudes após uma migração/remigração');

sectin=real(todos);

[nt nx]=size(todos);

for ix=1:nx;
    for it=1:nt; t=dt+(it-1)*dt;
        if(it<=1000);
        sectin(it,ix)=sectin(it,ix)*(t^2.21)*exp(-1.21*t);
        else;
        sectin(it,ix)=sectin(it,ix)*(t^2.21)*exp(-1.21*t);
        end;
    end;
end;

figure;
imagesc([1:nx]*dx,[42:nt]*dt,sectin(42:nt,:));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Seção remigrada - amplitudes equalizadas');