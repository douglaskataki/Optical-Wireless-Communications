clc;
clear;  
%% Propriedades no Detector
% potencia led
PLed = 20; % mW
% campo de visão do fotodetector
FOV = 60; % graus
% area do fotodetector
Area = 5.8E-6; % m^2
% semi angulo do Tx para meia potência
theta = 60; % graus
% n Ordem dem Emissão Lambertiana 
n = -log10(2)/log10(cosd(theta));
% Ganho optico no concentrador
index = 1.5;
G_Con = index^2/sin(FOV)^2;
TS = 1;
%% Posições
% tamanhos iniciais no local
xl = 5; yl = 5; zl = 3; 
% altura do Rx
h = 2.15;
% Posicao Tx
x_Tx = xl/2; y_Tx = yl/2; z_Tx = zl;
Pos_Tx = [x_Tx y_Tx z_Tx];
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
%% Cálculo da resposta para o Rx


% resposta do canal
H = HLOS(Area,FOV,n,N,Nx,Ny,x,y,N_Tx,Pos_Tx,N_Rx);
% potência de chegada
P_Rx = PLed*G_Con.*H;
P_Rx_dBm = 10*log(P_Rx);

% %% SNR
% % ruído
% % Dados retirados Fundamental Analysis for Visible-Light Communication System using LED Lights
% q = 1.6E-19; % Carga do eletron
% k = physconst('Boltzmann');
% c = physconst('LightSpeed');
% % banda (100Mb/s)
% B = 100E6; %Mbps
% % Responsividade do Fotodetector
% R = 0.54/1e-3; % A/W
% % photocorrent due to ground radiation
% I_b = 5400E-6; % uA
% % temperatura ambiente
% T = 25; % C
% T_k = T+273; % K
% % % open-loop voltage gain
% Gol = 10;
% % fixed capacitance of photodetector per unit area 
% Cpd = 112*10^-9/10^-4; %pF/cm^2
% % channel noise factor 
% Gamma = 1.5;
% % FET transconducatance  
% gm = 30E-3; % mS
% % noise bandwidth factor
% I_2 = .562;
% I_3 = .0868;
% % ruidos shot e thermal
% sigma2_shot = 2*q*R.*P_Rx*B+2*q*I_b*I_2*B; 
% sigma2_thermal = (8*pi*k*T_k/Gol)*Cpd*Ar*I_2*B^2 + (16*pi^2*k*T_k*Gamma/gm)*Cpd^2*Ar^2*I_3*B^3;
% noise2 = sigma2_shot+sigma2_thermal;
% SNR = R^2*P_Rx.^2./noise2;
% SNR_dB = 10*log10(SNR);
%% Plot do diagrama de potência em dBm
figure(1)
meshc(x,y,P_Rx_dBm);
grid on;
title('LOS - LED');
xlabel('x(m)');
ylabel('y(m)');
zlabel('Potência(dBm)')
colorbar;
colormap jet;