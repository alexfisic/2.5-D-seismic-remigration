disp('Esse programa deve ser rodado somente após o programa VIBROSEIS');

close all;
clear all; 

nt=256;
antinan=1e-06;
perc=.1;

x1=zeros(1,256);
taper=tuckey(nt,0.15);

perg=input('Escolha do sinal: 1-complicado; 2-classico; 3-cosseno; 4-chirp; 5-parabolico; 6-ricker:');

perg2=input('Tipo de metodo da transformada: 1-escalar; 2-vetor; 3-prod.convolucional:');

perg3=input('Realiza transformada de Gabor? 1-SIM; 2-NAO:');

perg4=input('Realiza transformada de ondaleta? 1-SIM; 2-NAO:');

if(perg==1); % Escolha do sinal

dt=1.01; 

t1=1:256;
x1(1:256)=cos(2*pi*(10+t1/7).*(t1/256))+cos(2*pi*((256/2.58)-t1/6).*t1/256);
t2=114:122;
x1(114:122)=x1(114:122)+cos(2*pi*t2*0.42);
t3=134:142;
x1(134:142)=x1(134:142)+cos(2*pi*t3*0.42);

elseif(perg==2);

dt=1.01; 

t1=1:128;
x1(1:128)=cos(2*pi*t1*(12/256));
t2=128:256;
x1(128:256)=cos(2*pi*t2*(50/256));
t3=20:30;
x1(20:30)=x1(20:30)+0.5*cos(2*pi*t3*(104/256));

elseif(perg==3);

dt=1.01; 

t1=1:128;
x1(1:128)=cos(2*pi*t1*(1/5));

t2=129:256;
x1(129:256)=cos(2*pi*t2*(2/5));


elseif(perg==4);

dt=0.0039; 

t1=[1:256]*dt;
arg1=(pi*(68/1)*t1)-(pi*(20/1)*(t1.^2));
arg2=2*pi*sin((5/1)*pi*t1)+(120/1)*pi*t1;
arg3=pi*(168/1)*t1+(28/1)*pi*(t1.^2);

x1(1:256)=cos(arg1)+cos(arg2)+cos(arg3);

elseif(perg==5);

dt=0.003962;

t1=[1:256]*dt;
arg1=144*pi*(t1-0.3).^3+30*pi*t1;
arg2=128*pi*t1-50*pi*t1.^2;

x1(1:256)=cos(arg1)+cos(arg2);

else;

dt=0.003962;
x1=ricker(nt,30,dt);

end; % Fim da escolha do sinal

x1=x1/max(x1);
x1=x1.*taper*1.0;

nf=nt;

U=zeros(256,nf);

pvetor=[1:256];

df=1/(256*dt);

h=fft(x1); h2=zeros(1,256);

%% === Informaçao de fase (TF) ===%%

ii=complex(0,1);
ang=angle(h); fase_reflec=unwrap(ang); 
arg=exp(ii*fase_reflec);

for i=1:256;h(1,i)=h(1,i)/abs(max(h)); end;

%%=== Transformada-S ===%%

if(perg2==1); % Tipo de método de transformada

for m=1:nf;
for p=1:256;
k = p + m; 
decay1=exp(-2*pi*(((p-256)^2)/m^2));
decay2=exp(-2*pi*(((1-p)^2)/m^2));
d=decay1;
if (k <= 256);  
h2(k)=h(p)*d;else;end; 
end;
U(:,m)=ifft(h2(1,:));
end;

elseif(perg2==2);

for m=1:nf;
for p=1:256;
k = p + m; 
if (k <= 256);  
h2(k)=h(p);else; end; 
end;

decay1=exp(-2*pi*(((pvetor-256).^2)/m^2));
d=decay1;

s=real(h2); t=imag(h2);
ii=complex(0,1);
amph=sqrt(s.^2+t.^2);
ang=angle(h2); fase_reflec=unwrap(ang); 
arg=exp(ii*fase_reflec);

U(:,m)=ifft(amph.*arg.*d);
end;

else;

for m=1:nf;
for p=1:256;
k = p + m; 
if (k <= 256);  
h2(k)=h(p);else; end; 
end;

vetor=(pvetor)*dt;
decay1=(abs(m*df)/sqrt(2*pi))*exp(-(0.5*(vetor.^2)*(m*df)^2));

dd=fft(decay1);

U(:,m)=ifft(h2.*dd);
end;

end; % Fim do método da transformada

%%=== Desespelhamento do sinal ===%

for p=1:nt;
U2(p,1:nf/2)=U(p,nf:-1:(nf/2)+1);
U2(p,(nf/2)+1:nf)=U(p,nf/2:-1:1);
end;

%%=== Espectro de Amplitude ===%

for m=1:nf;

a=real((U2(:,m)));  
b=imag((U2(:,m)));  

S(:,m)=sqrt(a.^2+b.^2)/nt; 
end;

%===== perc (%) da amplitude maxima =====%

for i=1:nt;
amp_max=max(S(i,:));
for j=1:nf;

if (S(i,j) >= perc*amp_max); 
ddS(i,j)=S(i,j)/amp_max;
else;
ddS(i,j)=0.0; 
end;

end;
end; 

%%=== Espectros de amplitude TVM ====%%

ampS=zeros(1,nf);
ampD=zeros(1,nf);

for i=1:nf;
  for j=1:nt;
      ampS(1,i)=ampS(1,i)+S(j,i)*dt;
      ampD(1,i)=ampD(1,i)+ddS(j,i)*dt;
  end; 
end; 

ampS=ampS/max(ampS);

%%=== Transformada de Gabor ====%%

if(perg3==1);

tg=[1:nt]*dt;
alpha=0.001;
G=zeros(nt,nf);

for itau=1:nt; tau = dt + (itau-1)*dt;

peso = exp(-(alpha*pi*(tg-tau).^2));

passante = fft(peso.*x1);

a = real(passante);
b = imag(passante);

G(itau,:) = sqrt(a.^2 + b.^2)/nt;

end;

else; end;

%%=== Transformada de ondaleta ====%%

if(perg4==1);

tg=[1:nt]*dt; dfw=1/(1.0*nt*dt*1.0);
nciclos=15;
ii=complex(0,1);
WT=zeros(nt,nf);

for itau=1:nt; tau = dt + (itau-1)*dt;

for ifi=1:nf; efe=dfw+(ifi-1)*dfw;

sigma=nciclos/(2*pi*efe); 

wmorlet = exp(-ii*2*pi*efe*((tg-tau)/1)).*exp(-(0.5*((1/sigma)^2)*(tg-tau).^2));

passante2 = wmorlet.*x1;
passante3=sum(passante2)*dt/sqrt(sigma);

WT(itau,ifi) = abs(passante3);

end;

end;

else; end;


%%=============================%%

figure;
subplot(2,1,1);
plot([1:256]*dt,x1); axis tight;grid;
subplot(2,1,2);
imagesc([1:256]*dt,[1:nf]*df,S');grid; 

figure;
imagesc([1:256]*dt,[1:nf]*df,ddS');grid; colorbar;

figure;
subplot(3,1,1);
plot([1:256]*df,ampS/nt); axis tight;grid;
subplot(3,1,2);
plot([1:256]*df,abs(h)/nt);axis tight; grid;
subplot(3,1,3);
plot([1:256]*df,ampD);axis tight; grid; 

if(perg3==1);
figure;
imagesc([1:256]*dt,[1:nf]*df,G');grid; colorbar;
else;end;

if(perg4==1);
figure;
imagesc([1:256]*dt,[1:nf]*dfw,WT');grid; colorbar;
else;end;

