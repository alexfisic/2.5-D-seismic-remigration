disp('Este programa deve ser rodado após terceiro');

itempo=250;
scaling_data=((1.41*dx)/sqrt(itempo*dt));
scaling_remig=((1.41*1.41*dx^(5/2)))/((itempo*dt)^(3/2));
scaling_migi=((1.41*dx)/(itempo*dt));
scaling_migic=((1.41*dx)/sqrt(itempo*dt));

figure;
plot([1:ntrace]*dx,real(migi(itempo,1:ntrace))*scaling_migi,...
[1:ntrace]*dx,real(migic(itempo,1:ntrace))*scaling_migic,'k',...
[1:ntrace]*dx,scaling_remig*real(migi3(itempo,1:ntrace)),'g-V'),...
axis tight,grid;
xlabel('Distance (meters)');ylabel('Amplitude');
title('1 - Picked amplitudes: data (magenta), migrated (green), remigrated (green diamonds)'); 
yyaxis right;
plot([1:ntrace]*dx,real(dado(itempo,1:ntrace))*scaling_data);
ylabel('Amplitude');


