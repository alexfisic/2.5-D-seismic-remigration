disp('Este programa deve ser rodados após terceiro_presal');

perc=input('Percentagem de patamar de amplitude (0.1,0.2,0.3 -- decimal)-->');

ddsm_kir=zeros(size(dado));
antinan=1e-3;

for i=1:nx;
amp_max_kir=max(abs(migic(:,i)));
for j=1:nt;

%===== perc da amplitude máxima =====%

if (migic(j,i) > perc*amp_max_kir); 
division_kir=migi(j,i)/migic(j,i);
ddsm_kir(j,i)=-division_kir;
else;
ddsm_kir(j,i)=antinan;
end;

end;
end;  

figure,imagesc([1:nx]*dx,[1:nt]*dt,-real(ddsm_kir));grid;colorbar;
xlabel('Distância (m)');ylabel('Tempo (s)');
title('Pré-sal - Atributo');




