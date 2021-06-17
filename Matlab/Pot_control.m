clc;
clear;
close all;
%% Propriedades no Detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potencia LED
PLed = 1000e-3; % mW
% campo de visao do fotodetector
FOV = 60; % graus
% area do fotodetector
Area = 5.8E-6; % m^2
% semi-angulo do Tx para meia potencia
theta = 70; % graus
% n Ordem dem Emissao Lambertiana 
n = -log10(2)/log10(cosd(theta));
% Ganho optico no concentrador
index = 1.5;
% G_Con = index^2/sin(FOV)^2;
TS = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes
c = physconst('LightSpeed')*1e-9; % tempo em nano segundos
% coeficiente de reflexão da parede, neste caso é considerado o mesmo para todas as paredes na sala.
rho = 0.8;

%% Posicoes Tx e Rx
% tamanhos iniciais da sala
xl = 5; yl = 5; zl = 3; 
% altura do Rx
%h = 2.15;

%% Taxa de transmissão
Rb = 200e6; % 200 Mbps

%% Tx1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Posições Tx e Rx
x_Tx = 2.269; y_Tx = 0.638; z_Tx = zl;
Pos_Tx = [x_Tx y_Tx z_Tx];

% x_Rx = 3.75; y_Rx = 4.75; z_Rx = 0;
x_Rx = 0; y_Rx = 0; z_Rx = 0;
Pos_Rx = [x_Rx y_Rx z_Rx];
% normaliza o locais para os Rx
M = 10;
Nx = xl*M; Ny = yl*M; Nz = round(zl*M);
x = linspace(0,xl,Nx);
y = linspace(0,yl,Ny);    
z = linspace(0,zl,Nz);
N = max([Nx Ny Nz]);
% base para os Rx
N_Tx = [0 0 -1];
N_Rx = -N_Tx;


%% HLOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vTxRx = Pos_Tx-Pos_Rx;
d = sqrt(sum((vTxRx).^2));
% ângulo entre Tx e Rx
cosphi = abs(dot(vTxRx,N_Rx))/d;
phi_los = acosd(cosphi);
if(phi_los<=FOV)
    costheta = abs(dot(vTxRx,N_Tx))/d;
    HLOS = (n+1)*Area*costheta^n/(2*pi*d^2);
else
    HLOS = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HNLOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = 0;
h2 = 0;
h3 = 0;
h4 = 0;
for kk=1:Ny
%%%%%%%%%%%%%%% plano X=0 %%%%%%%%%%%%%%%%%%
% fator de área
dA = zl*yl/(Ny*Nz);
% normal da parede
n1 = [1 0 0];
    for ll=1:Nz
        % ponto na parede Z,Y - Wall Point
        WP = [0 y(kk) z(ll)];
        vTxWP = Pos_Tx - WP; 
        % distancia do TX para a parede Z,Y(1)
        D1 = sqrt(sum((vTxWP).^2));
        % distancia TX e Rx
        % angulos Tx e incidencia na parede
        cos_phi = abs(dot(N_Tx,vTxWP))/D1;
        cos_alpha = abs(dot(vTxWP,n1))/D1;
        % distancia do WP para o Rx
        vWPRx = WP-Pos_Rx;
        D2 = sqrt(sum((vWPRx).^2));
        % angulos Rx e reflexão
        cos_psi = abs(dot(vWPRx,N_Rx))/D2;
        cos_beta = abs(dot(vWPRx,n1))/D2;
        if abs(acosd(cos_psi))<=FOV
            h1 = h1 +(n+1)*Area*rho*dA*...
            cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
        end
    end
end
%%%%%%%%%%%%%%% plano Y=0 %%%%%%%%%%%%%%%%%%
% fator de área
dA = zl*xl/(Nx*Nz);
% normal da parede
n2 = [0 1 0];
for kk=1:Nx
    for ll=1:Nz
        % ponto na parede Z,Y - Wall Point
        WP = [x(kk) 0 z(ll)];
        vTxWP = Pos_Tx - WP; 
        % distancia do TX para a parede Z,Y(1)
        D1 = sqrt(sum((vTxWP).^2));
        % distancia TX e Rx
        % angulos Tx e incidencia na parede
        cos_phi = abs(dot(N_Tx,vTxWP))/D1;
        cos_alpha = abs(dot(vTxWP,n2))/D1;
        % distancia do WP para o Rx
        vWPRx = WP-Pos_Rx;
        D2 = sqrt(sum((vWPRx).^2));
        % angulos Rx e reflexão
        cos_psi = abs(dot(vWPRx,N_Rx))/D2;
        cos_beta = abs(dot(vWPRx,n2))/D2;
        if abs(acosd(cos_psi))<=FOV
            h2  = h2 +(n+1)*Area*rho*dA*...
            cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
        end
    end
end
%%%%%%%%%%%%%%% plano X=xl %%%%%%%%%%%%%%%%%%
% fator de área
dA = zl*yl/(Ny*Nz);
% normal da parede
n3 = [-1 0 0];
for kk=1:Ny
    for ll=1:Nz
        % ponto na parede Z,Y - Wall Point
        WP = [xl y(kk) z(ll)];
        vTxWP = Pos_Tx - WP; 
        % distancia do TX para a parede Z,Y(1)
        D1 = sqrt(sum((vTxWP).^2));
        % distancia TX e Rx
        % angulos Tx e incidencia na parede
        cos_phi = abs(dot(N_Tx,vTxWP))/D1;
        cos_alpha = abs(dot(vTxWP,n3))/D1;
        % distancia do WP para o Rx
        vWPRx = WP-Pos_Rx;
        D2 = sqrt(sum((vWPRx).^2));
        % angulos Rx e reflexão
        cos_psi = abs(dot(vWPRx,N_Rx))/D2;
        cos_beta = abs(dot(vWPRx,n3))/D2;
        if abs(acosd(cos_psi))<=FOV
            h3 = h3 +(n+1)*Area*rho*dA*...
            cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
        end
    end
end        
%%%%%%%%%%%%%%% plano Y=yl %%%%%%%%%%%%%%%%%%
% fator de área
dA = zl*xl/(Nx*Nz);
% normal da parede
n4 = [0 -1 0];
for kk=1:Nx
    for ll=1:Nz
        % ponto na parede Z,Y - Wall Point
        WP = [x(kk) yl z(ll)];
        vTxWP = Pos_Tx - WP; 
        % distancia do TX para a parede Z,Y(1)
        D1 = sqrt(sum((vTxWP).^2));
        % distancia TX e Rx
        % angulos Tx e incidencia na parede
        cos_phi = abs(dot(N_Tx,vTxWP))/D1;
        cos_alpha = abs(dot(vTxWP,n4))/D1;
        % distancia do WP para o Rx
        vWPRx = WP-Pos_Rx;
        D2 = sqrt(sum((vWPRx).^2));
        % angulos Rx e reflexão
        cos_psi = abs(dot(vWPRx,N_Rx))/D2;
        cos_beta = abs(dot(vWPRx,n4))/D2;
        if abs(acosd(cos_psi))<=FOV
            h4 = h4 +(n+1)*Area*rho*dA*...
            cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
        end
    end
end        
HNLOS = h1+h2+h3+h4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Ruído
% ruido
% Dados retirados Fundamental Analysis for Visible-Light Communication System using LED Lights
q = 1.6E-19; % Carga do eletron
k = physconst('Boltzmann');
% c = physconst('LightSpeed');

% banda (100Mb/s)
B = Rb; %Mbps

% Responsividade do Fotodetector
R = 0.54; % A/W

% photocorrent due to ground radiation
I_b = 5400e-6; % uA

% temperatura ambiente
T = 25; % C
T_k = T+273; % K

% % open-loop voltage gain
Gol = 10;

% fixed capacitance of photodetector per unit area 
Cpd = 112*10^-12/10^-4; %pF/cm^2

% channel noise factor 
Gamma = 1.5;

% FET transconducatance  
gm = 30e-3; % mS

% noise bandwidth factor
I_2 = .562;
I_3 = .0868;

% ruido thermal

sigma2_thermal = (8*pi*k*T_k/Gol)*Cpd*Area*I_2*B^2 + (16*pi^2*k*T_k*Gamma/gm)*Cpd^2*Area^2*I_3*B^3;

% background (fundo)
% pbn = 5.8e-6/1e-4;  %uW/cm^2
% d_l = 30e-9; %nm
% sigma2_background = 2*q*I_2*pbn*Area*d_l*R*B;
% ruído devido a SSI (Solar Spectral Irradiance)
% I_sun = w(lambda)*delta_lambda
% w -> Solar Spectral Irradiance
% d_l -> largura de banda do OBPF que procede o photodetector
% ler arquivo SSI
% pegar o valor para w

data = 'ssi_v02r01_daily_s18820101_e18821231_c20170717.nc';

% matrix de dias por comprimento de onda
SSI = ncread(data,'SSI');

% valores de w são relativos a lambda em nm
lambda = ncread(data,'wavelength');
lambda_led = 902.5; %nm (próximo do comprimento de pico do LED)

i_ld = (lambda == lambda_led);
d_l = 30e-9; % 30 nm de banda
w = mean(SSI(i_ld,:));

sigma2_background = 2*q*B*R*w*d_l;

N0 = sigma2_thermal+sigma2_background;
% considerando que estamos utilizando OOK com Ps0 = 0 e P_Rx = Ps1

Pe_menor_erro = 1e-6;
Pe_maior_erro = 1e-4;

P_Tx_max = 3000e-3; %W
pace = 1e-3;

% vetor de tempo em horas
t = linspace(0,24,1e4); % 0 a 24 horas
%t = 1:1e-3:24;

% gaussiana com ruido de background
sgm = std(t);
pd = makedist('Normal',mean(t),sgm);
y = pdf(pd,t);
sigma2_background = 2*q*B*R*(y/max(y))*w*d_l;

P_Tx = PLed;
P_min = PLed;
HLOS_t = HLOS+HNLOS;
I_bg = 10e-9;
% alocação de memória
P_Rx_control = zeros(length(t),1);
N0_plot = zeros(length(t),1);
sigma2_shot = zeros(length(t),1);
Pe_plot = zeros(length(t),1);   
P_Tx_control = zeros(length(t),1); 

%% Controle de potência
for i=1:length(t)
    sigma2_shot(i) = 2*q*R*P_Tx*HLOS_t*B;
    
    % memória
    P_Tx_ant = P_Tx;
    
    % atualiza o valor de potência
    P_Rx_control(i) = P_Tx*HLOS_t;
    
    % para o caso do OOK-NRZ
    Eb = 2*(P_Rx_control(i)*R)^2;
    
    % atualizando o valor do ruído
    if i==1
        N0 = N0+sigma2_background(i); %+sigma2_shot(i);
    else
        % isto é feito porque somente o valor do novo do ruído deve ser
        % considerado
        
        N0 = N0+sigma2_background(i)-sigma2_background(i-1)+sigma2_shot(i)-sigma2_shot(i-1);
    end
    
    % verifica o valor para Pe e a partir desse valor, ele faz o controle
    % para a potência
    Pe = qfunc(sqrt(Eb/N0));
    Pe_plot(i) = Pe;
    
    %% Controle da potência
    %%% caso esteja menor que a SNR, então aumenta a potência
    if(Pe >= Pe_maior_erro)
        % aumenta potência no TX
        P_Tx = P_Tx + pace;
        % verifica saturação no TX
        if(P_Tx >= P_Tx_max) 
            P_Tx = P_Tx_max;
        end
        %atualiza potência que chega no RX
        P_Rx_control(i) = P_Tx*HLOS_t;

        %%% Limite inferior
    elseif(Pe <= Pe_menor_erro)
        % diminui potência se Probabilidade de erro for menor que o
        % estimado
        P_Tx = P_Tx - pace;
        % no máximo para a potência mínima
        if(P_Tx <= P_min)
            P_Tx = P_min;
        end
        
        % atualiza potência que chega no Rx
        P_Rx_control(i) = P_Tx*HLOS_t;
        %P_Rx_control(i) = P_Tx;
    else
        % mantém a potência potência
        P_Rx_control(i) = P_Tx_ant*HLOS_t;
    end
    N0_plot(i) = N0;
    P_Tx_control(i) = P_Tx;
end

figure(1)
title('Potencia');
plot(t,P_Rx_control);
ylabel('P Rx(W)');
xlabel('tempo(h)');
xlim([0 24])
%ylim([0 P_Tx_max])
legend('Pot Control');
grid on;

figure(2)
plot(t,Pe_plot,t,Pe_menor_erro*ones(length(t),1),t,Pe_maior_erro*ones(length(t),1))
grid on;
ylabel('P_e');
xlabel('tempo(h)');
xlim([0 24])
% ylim([Pe_low 1e-1]);

figure(3)
plot(t,P_Tx_control);
grid on;
ylabel('P Tx(W)');
xlabel('tempo(h)');
% xlim([0 24])  
% ylim([0 P_Tx_max]);

% legend('PTx','Pmax','Pmin');

% figure(3)
% title("/sigma_background");
% plot(t,sigma2_background);
% xlabel('tempo(h)');
% grid on;
% xlim([0 24])
% 
figure(4)
plot(t,N0_plot);
xlabel('tempo(h)');
ylabel('N_0');
grid on;
xlim([0 24])