disp('Tenta equalizar as amplitudes após uma migração/remigração');

sectin=real(ic);

for ix=1:nx;
    sectin(:,ix)=sectin(:,ix)/max(abs(sectin(:,ix)));
end;

for ix=1:nx;
    for it=1:nt; t=dt+(it-1)*dt;
        if(it<=990);
        sectin(it,ix)=sectin(it,ix)*(t^0.21)*exp(-4.2*t);
        else;
        sectin(it,ix)=sectin(it,ix)*(t^0)*exp(-0.0*t);
        end;
    end;
end;

figure;
imagesc([1:nx]*dx,[1:nt]*dt,-sectin(1:nt,:));colormap gray;grid;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Seção remigrada - amplitudes equalizadas');