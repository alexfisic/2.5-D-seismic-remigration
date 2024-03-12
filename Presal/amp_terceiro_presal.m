disp('Este programa deve ser rodado após terceiro_presal');

itempo=1094;
scaling_data=1.41*sqrt(dx);
scaling_remig=2*1.41*(dx/(itempo*dt));
scaling_migi=2*sqrt(dx)/sqrt(itempo*dt);
scaling_migic=1.41*sqrt(dx);

figure;
plot([1:ntrace]*dx,real(migi(itempo,1:ntrace))*scaling_migi,'b',...
[1:ntrace]*dx,real(migic(itempo,1:ntrace))*scaling_migic,'k',...
[1:ntrace]*dx,scaling_remig*real(migi3(itempo,1:ntrace)),'g-V'),...
axis tight,grid;
xlabel('Distance (meters)');ylabel('Amplitude');
title('1 - Picked amplitudes: data (magenta), migrated (green), remigrated (green diamonds)'); 
yyaxis right;
plot([1:ntrace]*dx,real(dado(itempo,1:ntrace))*scaling_data);
ylabel('Amplitude');
