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
theta = 60; % graus
% n Ordem dem Emissao Lambertiana
n = -log10(2)/log10(cosd(theta));
% Ganho optico no concentrador
index = 1.5;
% G_Con = index^2/sin(FOV)^2;
TS = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constantes
% c = physconst('LightSpeed')*1e-9; % tempo em nano segundos
% coeficiente de reflexão da parede, neste caso é considerado o mesmo para todas as paredes na sala.
% rho = 0.8;

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
x_Tx = xl/2; y_Tx = yl/2; z_Tx = zl;
Pos_Tx = [x_Tx y_Tx z_Tx];

vy = 10e-2;
t = 0:yl/vy; %segundos
%fixos
x_Rx = xl/2; y_Rx = 0;  z_Rx = 0;
Pos_Rx = [x_Rx y_Rx z_Rx];
% normais
% rotx = [1 0 0; 0 cosd(angulosTx(1)) -sind(angulosTx(1)); 0 sind(angulosTx(1)) cosd(angulosTx(1))];
% roty = [cosd(angulosTx(2)) 0 sind(angulosTx(2)); 0 1 0; -sind(angulosTx(2)) 0 cosd(angulosTx(2))];
% rotz = [cos(angulosTx(3)) -sin(angulosTx(3)) 0; sin(angulosTx(3)) cos(angulosTx(3)) 0; 0 0 1];
% N_Tx = rotz*roty*rotx*N_Tx;
N_Rx = [0 0 1];

% mudar Tx para -20 graus em y
%angle_init = -atand((xl/2)/zl);
%N_Tx = [0 0 -1]*[1 0 0; 0 cosd(angle_init) -sind(angle_init); 0 sind(angle_init) cosd(angle_init)];

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
P_min = 100e-3;
P_Tx_max = 10000e-3; %W

% Condição inicial
P_Tx = PLed;
pace = 100e-3; %100 mW
% rotação inicial de ângulos, considerando que ele rotaciona até o alvo
% FOV_l = -FOV/2+abs(angle_init);
% FOV_h = FOV/2-abs(angle_init);

% valores para SNR
Pe_menor_erro = 1e-6;
Pe_maior_erro = 1e-4;

% potência relativas aos erros
P_Rx_maior_erro = sqrt((N0*(qfuncinv(Pe_maior_erro))^2)/2)*(1/R);
P_Rx_menor_erro = sqrt((N0*(qfuncinv(Pe_menor_erro))^2)/2)*(1/R);

% alocação de memória para melhorar a velocidade do programa
P_Rx = zeros(1,length(t));
Pe_plot = zeros(1,length(t));
HLOS = zeros(1,length(t));
HNLOS = zeros(1,length(t));
check = zeros(1,length(t));
P_Tx_control = zeros(1,length(t));
costheta = zeros(1,length(t));

% para o reflexo de 1st ordem nas paredes
M = 10;
Nx = xl*M; Ny = yl*M; Nz = round(zl*M);
x = linspace(0,xl,Nx);
y = linspace(0,yl,Ny);
z = linspace(0,zl,Nz);
N = max([Nx Ny Nz]);

%% Constantes
% coeficiente de reflexão da parede, neste caso é considerado o mesmo para todas as paredes na sala.
rho = 0.8;
N_Tx_0 = [0 0 -1];

%lembrar que depende da atualização de potência
for i=1:length(t)
    
    vTxRx = Pos_Tx-Pos_Rx;
    d = sqrt(sum((vTxRx).^2));
    N_Tx = -vTxRx/d;
    %N_Tx = N_Tx_0;
    % ângulo entre Tx e Rx
    cosphi = dot(vTxRx,N_Rx)/d;
    phi_los = acosd(cosphi);
    % fator para saber o novo phi_los rotacionado
        
    if(phi_los<FOV)
        costheta = abs(dot(vTxRx,N_Tx))/d;
        HLOS(i) = (n+1)*Area*costheta^n/(2*pi*d^2);
    else
        HLOS(i) = 0;
    end
    %%%%%%% HLOS %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % atualiza a posição de RX
    y_Tx = y_Tx + vy; % vx metros em um segundo
    Pos_Tx = [x_Rx y_Rx z_Rx];
end


% figure(1)
% %subplot(312)
% p = plot(t,P_Rx,t,P_Rx_menor_erro*ones(1,length(t)),t,P_Rx_maior_erro*ones(1,length(t)));
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% grid on;
% %title('b) Received Power');
% ylabel('P_{Rx}(W)')
% xlabel('time(s)');
% legend('P_{Rx}','P_{low}','P_{high}');
%
% %subplot(313)
% figure(2)
% % p = plot(t,Pe_plot,t,Pe_menor_erro*ones(length(t),1),t,Pe_maior_erro*ones(length(t),1));
% p = plot(t,log10(Pe_plot),'r');
% p.LineWidth = 2;
% % p(2).LineWidth = 2;
% % p(3).LineWidth = 2;
% %title('c) bit error rate');
% ylabel('log(BER)');
% xlabel('time(s)');
% % ylim([Pe_menor_erro 1.5*Pe_maior_erro]);
% grid on;
% % legend('P_e','P_{el}','P_{eh}')
%
% %subplot(311)
% figure(3)
% p = plot(t,P_Tx_control,'g');
% p(1).LineWidth = 2;
% grid on;
% %title('a) Transmitted power');
% xlabel('time(s)');
% ylabel('P_{Tx}(W)');
