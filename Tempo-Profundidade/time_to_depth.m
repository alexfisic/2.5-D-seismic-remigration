%                    MÓDULO 1: CONVERSÃO TEMPO-PROFUNDIDADE                                            
% Script para fazer a conversão tempo profundidade com raio normal a partir
% de uma seção empilhada após 'picking'.
% MARCELO MESQUITA - 2015.

close all;
clear all;
clc;

v1cam=1500;
v2cam=1800;
v3cam=2100;


load 'picktempo.dat';                       % Leitura do arquivo com 'picks' no dado real.
dcmp=25;                                    % Distancia entre tracos cmp.
ptrac=1;intrac=20;utrac=401;                % intervalo entre os tracos,primeiro traco, ultimo traco.
nt=(utrac-ptrac)+1;                         % Numero de tracos.
cmp_x=dcmp*[(ptrac-1):intrac:(utrac-1)];    % conversão cmp para x (offset).
t1=picktempo(1,[ptrac:intrac:utrac]);       % Leirura do vetor tempo duplo para a interface 2.

% Linha entre primeiro e último ponto.
sp_cmp_x=linspace(cmp_x(1),cmp_x(length(cmp_x)),dcmp*nt);
% =========================================================================
%               Tracamento de raios para a interface 1:
% =========================================================================
dg = 1; dt = 0.004; tmax=t1/2;          % Espacamento entre os nos; Intervalo de tempo; Tempo de trânsito simples.
ang_trac1=zeros(length(t1)+1);   %[[ang1] ang1_end]; % Angulo de partida de cada raio;
zs=0;                                   %Coordenadas iniciais de x e z da fonte.
xs = cmp_x; nz=1500; nx=xs(end);        % Vetor posição das fontes; Profundidade do modelo (z).
x=xs(1):dg:nx; z=0:dg:nz; v0 = v1cam;   % Tamanho do modelo (matriz); Velocidade inicial.
vmodel=ones(nz,nx)*v1cam;               % Preenchendo a matriz com v1.
rayvelmod(vmodel,dg); clear vmodel      % Inicializa o modelo de velocidade.

for i=1:length(tmax)                 % Laco para insercao de varias fontes.
   % if (xs > nt);break;end  (Retirar para aceitar modelos menores).
    tstep=0:dt:tmax(i);          
    theta=pi*ang_trac1(i)/180;                     % Mudan�a de Radianos para graus.
    r0=[xs(i),zs,sin(theta)/v0,cos(theta)/v0]';    % Valores iniciais do vetor r.
    [t,r]=ode45('drayvec',tstep,r0);               % Resolve as eq. pelo metodo RK de 4 ordem     
    coord_x1=r(:,1);                               % Coordenda 'x' do raio. 
    coord_z1=r(:,2);                               % Coordenada 'z' do raio.   
    xx1(i)=coord_x1(end);                          % Coordenada final de 'x'.
    zz1(i)=coord_z1(end);                          % Coordenada final de 'z'.
end
int1=csapi(xs,zz1,sp_cmp_x); % Interpolando os pontos da interface 1.

% =========================================================================
%               Tracamento de raios para a interface 2:
% =========================================================================
clear rayvelmod(vmodel,dg);

t2=picktempo(2,[ptrac:intrac:utrac]);       % Leirura do vetor tempo duplo para a interface 2.
tmax2=t2/2;                                 % Tempo de trânsito simples.
% Linha entre primeiro e último ponto.
sp_cmp_x=linspace(cmp_x(1),cmp_x(length(cmp_x)),dcmp*nt);
ang_trac2=zeros(length(t2)+1);   %[[ang2] ang2_end];% Angulo de partida de cada raio para a interface 2;
% Gerando modelo de velocidade 2:
vmodel2=zeros(nz,nx);
ind=zeros(nz,nx);
for i1=1:nx
    for i2=1:nz
        if(i2 <= int1(i1))
            vmodel2(i2,i1)=v1cam;    ind(i2,i1)=1;
        else
            vmodel2(i2,i1)=v2cam;    ind(i2,i1)=2;
        end
    end
end

% ========== Suavizando o modelo ============
suavx=300; suavz=300;                   % Fatores de suavização.
vsm2=smooth2a(vmodel2,suavx,suavz);     % Suavizando interfaces.   
rayvelmod(vsm2,dg); clear vsm2 vmodel2        % Inicializa o modelo de velocidade.
% ============================================

for i=1:length(tmax2)                   % Laco para inserçao de varias fontes.
%    if (xs > 458);break;end            (Retirar para aceitar modelos menores).
    tstep=0:dt:tmax2(i);          
        theta=pi*ang_trac2(i)/180;                     % Mudan�a de Radianos para graus.
        r0=[xs(i),zs,sin(theta)/v0,cos(theta)/v0]';    % Valores iniciais do vetor r.
        [t,r]=ode45('drayvec',tstep,r0);               % Resolve as equacoes pelo metodo RK de 4 ordem
                                                       % e nos da a trajetoria do raio.                                                       
        coord_x2=r(:,1);                               % Coordenda 'x' do raio. 
        coord_z2=r(:,2);                               % Coordenada 'z' do raio.
        xx2(i)=coord_x2(end);
        zz2(i)=coord_z2(end);                          % Interface 2.
end
int2=csapi(xs,zz2,sp_cmp_x);    % Vetor profundidade da interface 2 interpolado.

% =========================================================================
%               Tracamento de raios para a interface 3:
% =========================================================================

figure(1)
plot(int1,'LineWidth',2);set(gca,'ydir','reverse'); hold on
plot(int2,'LineWidth',2);
hold on;

%==========================================================================
%           Parâmetros estruturais reais do modelo
%==========================================================================

x1_real=[0 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
z1_real=[500 400 300 350 500 650 700 650 500 450 500];

x2_real=[0 1000 2000 3000 4000 5000 6000 7000 8000 10000];
z2_real=[1000 1100 1200 1250 1200 1100 1000 930 900 1000];

sp_real_x1=linspace(x1_real(1),x1_real(end),dcmp*nt);
sp_real_x2=linspace(x2_real(1),x2_real(end),dcmp*nt);

int1_real=csapi(x1_real,z1_real,sp_real_x1);
int2_real=csapi(x2_real,z2_real,sp_real_x2);

plot(int1_real,'r--','LineWidth',2);set(gca,'ydir','reverse');
plot(int2_real,'r--','LineWidth',2);set(gca,'ydir','reverse');
title('Interface real - tracejado'); xlabel('Distancia (m)'); ylabel('Profundidade (m)')

