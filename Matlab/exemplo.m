clc;
clear;  
%% Propriedades no Detector
% Potencia LED
PLed = 20; % mW
% campo de visao do fotodetector
FOV = 60; % graus
% area do fotodetector
Area = 5.8E-6; % m^2
% semi-angulo do Tx para meia potencia
theta = 60; % graus
% n Ordem dem Emissao Lambertiana 
n = -log10(2)/log10(cosd(theta));
% Ganho optico no concentrador
index = 1.5;
G_Con = index^2/sin(FOV)^2;
TS = 1;
%% Posicoes Tx e Rx
% tamanhos iniciais da sala
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

%% Matrizes de Rotações
% ângulo de rotação, lembrar de colocar em graus
% angulosTx = [0 0 0];
% angulosRx = [0 0 0];

% Tx
% rotx = [1 0 0; 0 cosd(angulosTx(1)) -sind(angulosTx(1)); 0 sind(angulosTx(1)) cosd(angulosTx(1))];
% roty = [cosd(angulosTx(2)) 0 sind(angulosTx(2)); 0 1 0; -sind(angulosTx(2)) 0 cosd(angulosTx(2))];
% rotz = [cos(angulosTx(3)) -sin(angulosTx(3)) 0; sin(angulosTx(3)) cos(angulosTx(3)) 0; 0 0 1];
% N_Tx = rotz*roty*rotx*N_Tx;

% Rx
% rotx = [1 0 0; 0 cosd(angulosRx(1)) -sind(angulosRx(1)); 0 sind(angulosRx(1)) cosd(angulosRx(1))];
% roty = [cosd(angulosRx(2)) 0 sind(angulosRx(2)); 0 1 0; -sind(angulosRx(2)) 0 cosd(angulosRx(2))];
% rotz = [cos(angulosRx(3)) -sin(angulosRx(3)) 0; sin(angulosRx(3)) cos(angulosRx(3)) 0; 0 0 1];
% N_Rx = rotz*roty*rotx*N_Rx;

%% Calculo da Resposta do Canal
H = HLOS(Area,FOV,n,N,Nx,Ny,x,y,N_Tx,Pos_Tx,N_Rx);
% potencia de chegada
P_Rx = PLed*G_Con.*H;
% para dBm
P_Rx_dBm = 10*log10(P_Rx);

%% SNR
% ruido
% Dados retirados Fundamental Analysis for Visible-Light Communication System using LED Lights
q = 1.6E-19; % Carga do eletron
k = physconst('Boltzmann');
c = physconst('LightSpeed');
% banda (100Mb/s)
B = 100E6; %Mbps
% Responsividade do Fotodetector
R = 0.54; % A/W
% photocorrent due to ground radiation
I_b = 5400E-6; % uA
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
gm = 30E-3; % mS
% noise bandwidth factor
I_2 = .562;
I_3 = .0868;
% ruidos shot e thermal
sigma2_shot = 2*q*R.*P_Rx*B+2*q*I_b*I_2*B; 
sigma2_thermal = (8*pi*k*T_k/Gol)*Cpd*Area*I_2*B^2 + (16*pi^2*k*T_k*Gamma/gm)*Cpd^2*Area^2*I_3*B^3;
% background (fundo)
pbn = 5.8e-6/1e-4;  %uW/cm^2
d_l = 30e-9; %nm
sigma2_background = 2*q*I_2*pbn*Area*d_l*R*B;
noise2 = sigma2_shot+sigma2_thermal+sigma2_background;
% considerando que estamos utilizando OOK com Ps0 = 0 e P_Rx = Ps1
SNR = R^2*P_Rx.^2./noise2;
SNR_dB = 10*log10(SNR);

%% Plot do diagrama de potencia em dBm
subplot(1,2,1)
meshc(x,y,P_Rx_dBm);
grid on;
title('LOS - LED');
xlabel('x(m)');
ylabel('y(m)');
zlabel('Potencia(dBm)');
colorbar;

subplot(1,2,2)
% com reflexao
rho = 0.8;
H = H + HNLOS(xl,yl,zl,rho,Area,FOV,n,Nx,Ny,Nz,x,y,z,N_Tx,Pos_Tx,N_Rx);
P_Rx = PLed*G_Con.*H;
P_Rx_dBm = 10*log10(P_Rx);
meshc(x,y,P_Rx_dBm);
colorbar;
colormap jet;

figure(2)
meshc(x,y,SNR_dB);
title('LOS - LED');
xlabel('x(m)');
ylabel('y(m)');
zlabel('SNR(dBm)')
colorbar;
colormap jet;