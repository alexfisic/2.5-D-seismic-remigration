disp('An�lise espectral do dado s�smico do presal');

load kir_presal.dat; 
sect = reshape(kir_presal,1750,720);
[nt nx]=size(sect);

dt = 0.004;
dx = 50;
df = 1/(nt*dt);

for ix=1:nx;
   sect(:,ix)=sect(:,ix)/max(abs(sect(:,ix)));
end;

data = zeros(nt,nx);
data2 = zeros(nt,nx);

spec_mean = zeros(nt,1);
vec_freq  = zeros(1,nx);

for ix=1:nx;

y = fft(sect(:,ix));

a = real(y); 
b = imag(y);

amp = sqrt(a.*a + b.*b);

d = amp.*amp;

data (1:nt,ix) = amp(1:nt); % ESPECTRO DE AMPLITUDE

data2(1:nt,ix) = d(1:nt);   % ESPECTRO DE POTÊNCIA

end;


for it=1:round(nt/8);

for ix=1:nx;
 spec_mean(it,1) = spec_mean(it,1) + data(it,ix)*dx;
end;
 
end;

spec_mean = spec_mean/((nx-1)*dx);

% Normalização do Vetor de Espectro de Freqüência

for ix=1:nx;

for it=1:nt/2;

vec_freq(1,ix) = vec_freq(1,ix) + data2(it,ix)*df;

end;

end;

% Vetor Freqüência Central

freq_central = zeros(1,nx);

for ix=1:nx;

for it=1:nt/2;

f = df + (it-1)*df;

freq_central(1,ix) = freq_central(1,ix) + f*data2(it,ix)*df;

end;

freq_central(1,ix) = freq_central(1,ix)/vec_freq(1,ix);

end;

% Vetor Banda de Freqüência

banda = zeros(1,nx);

for ix=1:nx;

for it=1:nt/2;

f = df + (it-1)*df;

banda(1,ix) = banda(1,ix) + ((f-freq_central(1,ix))^2)*data2(it,ix)*df;

end;

banda(1,ix) =sqrt(banda(1,ix)/vec_freq(1,ix));

end;

% Vetor Freqüência Dominante

dom_freq = zeros(1,nx);

for ix=1:nx;

for it=1:nt/2;

f = df + (it-1)*df;

dom_freq(1,ix) = dom_freq(1,ix) + (f^2)*data2(it,ix)*df;

end;

dom_freq(1,ix) = sqrt(dom_freq(1,ix)/vec_freq(1,ix));

end;

f_central = 0.0;
f_dominante = 0.0;
banda_dado = 0.0;

for ix=1:nx;
 f_central = f_central + freq_central(1,ix)*dx;
 f_dominante = f_dominante + dom_freq(1,ix)*dx;
 banda_dado = banda_dado + banda(1,ix)*dx;
end;

f_central = f_central/((nx-1)*dx);
f_dominante = f_dominante/((nx-1)*dx);
banda_dado = banda_dado/((nx-1)*dx);

disp('frequencia central='), disp(f_central)
disp('frequencia dominante='), disp(f_dominante)
disp('intervalo de banda do dado [f_min,f_max] ='), disp(f_central-banda_dado), disp(f_central+banda_dado)

%%%=========== FIGURAS ============%%%

% figure, plot([1:nt/4]*df,spec_mean(1:207,1)); grid;

figure, plot([1:nt/4]*df,spec_mean(1:nt/4,1),[f_central-banda_dado f_central-banda_dado],[0 12],'g',[f_central+banda_dado f_central+banda_dado],[0 12],'g',[f_central f_central],[0 12],'r',[f_dominante f_dominante],[0 12],'k'); grid;

figure, plot([1:nt/4]*df,spec_mean(1:nt/4,1),[1:nt/4]*df,amp(1:nt/4,:)); grid;

figure, plot([1:nx]*dx,freq_central,'b',[1:nx]*dx,banda,'g:d',[1:nx]*dx,dom_freq,'r');grid;


%% Como plotar 2 barras verticais junto do comando plot do Matlab %%

%%figure, plot([1:nt/4]*df,spec_mean(1:207,1),[f_central-banda_dado f_central-banda_dado],[0 12],[f_central+banda_dado f_central+banda_dado],[0 12]); 
