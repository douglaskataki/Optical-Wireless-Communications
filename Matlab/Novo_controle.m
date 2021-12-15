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
ang = linspace(-atand((yl/2)/zl),atand((yl/2)/zl),length(t));
clicks = zeros(1,length(t));
clicks_up = zeros(1,length(t));
clicks_down = zeros(1,length(t));
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

count_1 = 0;
count_2 = 0;
count_3 = 0;
%lembrar que depende da atualização de potência
for i=1:length(t)
    % memória de potência
    P_Tx_control(i) = P_Tx;
    %    P_Tx_ant = P_Tx;
    
    % atualiza ângulo
    %N_Tx = [1 0 0; 0 cosd(ang(i)) -sind(ang(i)); 0 sind(ang(i)) cosd(ang(i))]*[0 0 -1]';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% HLOS %%%%%%
    vTxRx = Pos_Tx-Pos_Rx;
    d = sqrt(sum((vTxRx).^2));
    N_Tx = -vTxRx/d;
    % ângulo entre Tx e Rx
    cosphi = dot(vTxRx,N_Rx)/d;
    phi_los = acosd(cosphi);
    % fator para saber o novo phi_los rotacionado
    
    if(phi_los<FOV)
        costheta(i) = abs(dot(vTxRx,N_Tx))/d;
        HLOS(i) = (n+1)*Area*costheta(i)^n/(2*pi*d^2);
    else
        HLOS(i) = 0;
    end
    %%%%%%% HLOS %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% HNLOS %%%%%%%%%%
    if 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        HNLOS(i) = h1+h2+h3+h4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    P_Rx(i) = P_Tx_control(i)*(HLOS(i)+HNLOS(i));
    
    % atualiza a posição de RX
    y_Rx = y_Rx + vy; % vx metros em um segundo
    Pos_Rx = [x_Rx y_Rx z_Rx];
    
    % atualiza normal Tx
    %ang(i) = -atand((xl/2-vx*i)/zl);
    
    % ajusta FOV
    %FOV_l = FOV_l + ang;
    %FOV_h = FOV_h + ang;
    
    % considerando OOK, ajusta potência se for depender dos valores da SNR
    Eb = 2*(P_Rx(i)*R)^2;
    Pe = qfunc(sqrt(Eb/N0));
    
    %% pace variável
    if(Pe >= Pe_maior_erro)
        % aumenta potência no TX para o valor de P Tx para o alcançar o
        % maior erro possível
        P_Tx = P_Rx_maior_erro/(HLOS(i)+HNLOS(i));
        % verifica saturação no TX
        if(P_Tx >= P_Tx_max)
            P_Tx = P_Tx_max;
        end
        %atualiza potência que chega no Rx
        count_1 = count_1+1;
        %%% Limite inferior
    elseif(Pe <= Pe_menor_erro)
        % diminui potência se Probabilidade de erro for menor que o
        % estimado
        P_Tx = P_Rx_menor_erro/(HLOS(i)+HNLOS(i));
        % no máximo para a potência mínima
        if(P_Tx <= P_min)
            P_Tx = P_min;
        end
        
        % atualiza potência que chega no Rx
        count_2 = count_2+1;
    else
        %caso contrário, mantém a potência potência
        P_Tx = P_Tx_control(i);
        
        count_3 = count_3+1;
    end
    P_Rx(i) = P_Tx*(HLOS(i)+HNLOS(i));
    Eb = 2*(P_Rx(i)*R)^2;
    Pe_plot(i) = qfunc(sqrt(Eb/N0));
    % calcula o HLOS(i) + HNLOS(?) -> volta para o começo do loop
end

%%%% sem variacao normal
P_Tx = PLed;
y_Rx = 0;
P_Tx_control_2 = zeros(length(t),1);
P_Rx_2 = zeros(length(t),1);
Pe_plot_2 = zeros(length(t),1);

N_Tx = [0 0 -1];
for i=1:length(t)
    % memória de potência
    P_Tx_control_2(i) = P_Tx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% HLOS %%%%%%
    vTxRx = Pos_Tx-Pos_Rx;
    d = sqrt(sum((vTxRx).^2));
    % ângulo entre Tx e Rx
    cosphi = dot(vTxRx,N_Rx)/d;
    phi_los = acosd(cosphi);
    % fator para saber o novo phi_los rotacionado
    
    if(phi_los<FOV)
        costheta(i) = abs(dot(vTxRx,N_Tx))/d;
        HLOS(i) = (n+1)*Area*costheta(i)^n/(2*pi*d^2);
    else
        HLOS(i) = 0;
    end
    %%%%%%% HLOS %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%% HNLOS %%%%%%%%%%
    if 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        HNLOS(i) = h1+h2+h3+h4;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    P_Rx_2(i) = P_Tx_control_2(i)*(HLOS(i)+HNLOS(i));
    
    % atualiza a posição de RX
    y_Rx = y_Rx + vy; % vx metros em um segundo
    Pos_Rx = [x_Rx y_Rx z_Rx];
    
    % atualiza normal Tx
    %ang(i) = -atand((xl/2-vx*i)/zl);
    
    % ajusta FOV
    %FOV_l = FOV_l + ang;
    %FOV_h = FOV_h + ang;
    
    % considerando OOK, ajusta potência se for depender dos valores da SNR
    Eb = 2*(P_Rx_2(i)*R)^2;
    Pe = qfunc(sqrt(Eb/N0));
    
    %% pace variável
    if(Pe >= Pe_maior_erro)
        if i==1
            clicks(i) = ceil((P_Rx_maior_erro/(HLOS(i)+HNLOS(i))-P_Tx_control_2(i))/pace);
        else
            clicks(i) = ceil((P_Rx_maior_erro/(HLOS(i)+HNLOS(i))-P_Tx_control_2(i))/pace) - clicks(i-1);
        end
        % quantidade de clicks para ajuste da potência
        %         clicks_up(i) = ceil((P_Rx_maior_erro/(HLOS(i)+HNLOS(i))-P_Tx_control(i))/pace);
        
        % aumenta potência no TX para o valor de P Tx para o alcançar o
        % maior erro possível
        P_Tx = P_Rx_maior_erro/(HLOS(i)+HNLOS(i));
        % verifica saturação no TX
        if(P_Tx >= P_Tx_max)
            P_Tx = P_Tx_max;
        end
        %atualiza potência que chega no Rx
        P_Rx_2(i) = P_Tx*(HLOS(i)+HNLOS(i));
        
        %%% Limite inferior
    elseif(Pe <= Pe_menor_erro)
        if i==1
            clicks(i) = ceil((P_Rx_menor_erro/(HLOS(i)+HNLOS(i))- P_Tx_control_2(i))/pace);
        else
            clicks(i) = ceil((P_Rx_menor_erro/(HLOS(i)+HNLOS(i))- P_Tx_control_2(i))/pace) - clicks(i);
        end
        % quantidade de clicks para ajuste da potência
        %         clicks_down(i) = ceil((P_Tx_control(i)-P_Rx_menor_erro/(HLOS(i)+HNLOS(i)))/pace);
        % diminui potência se Probabilidade de erro for menor que o
        % estimado
        P_Tx = P_Rx_menor_erro/(HLOS(i)+HNLOS(i));
        % no máximo para a potência mínima
        if(P_Tx <= P_min)
            P_Tx = P_min;
        end
        
        % atualiza potência que chega no Rx
        P_Rx_2(i) = P_Tx*(HLOS(i)+HNLOS(i));
        %P_Rx_control(i) = P_Tx;
    else
        % mantém a potência potência
        P_Tx = P_Tx_control_2(i);
        P_Rx_2(i) = P_Tx*(HLOS(i)+HNLOS(i));
    end
    Eb = 2*(P_Rx_2(i)*R)^2;
    Pe_plot_2(i) = qfunc(sqrt(Eb/N0));
    % calcula o HLOS(i) + HNLOS(?) -> volta para o começo do loop
end

figure(1)
%subplot(312)
p = plot(t,P_Rx,t,P_Rx_menor_erro*ones(1,length(t)),t,P_Rx_maior_erro*ones(1,length(t)));
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
grid on;
%title('b) Received Power');
ylabel('P_{Rx}(W)')
xlabel('time(s)');
legend('P_{Rx}','P_{low}','P_{high}');

%subplot(313)
figure(2)
% p = plot(t,Pe_plot,t,Pe_menor_erro*ones(length(t),1),t,Pe_maior_erro*ones(length(t),1));
p = plot(t,log10(Pe_plot),'r');
p.LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
%title('c) bit error rate');
ylabel('log(BER)');
xlabel('time(s)');
% ylim([Pe_menor_erro 1.5*Pe_maior_erro]);
grid on;
% legend('P_e','P_{el}','P_{eh}')

%subplot(311)
figure(3)
p = plot(t,P_Tx_control,'g');
p(1).LineWidth = 2;
grid on;
%title('a) Transmitted power');
xlabel('time(s)');
ylabel('P_{Tx}(W)');

% p = plot(t,P_Rx_2,t,P_Rx_menor_erro*ones(1,length(t)),t,P_Rx_maior_erro*ones(1,length(t)));
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% grid on;
% title('Power on RX after regulation on TX, non rotation');
% ylabel('Power(W)')
% xlabel('time(s)');
% legend('Power Rx','Power lowest error','Power biggest error');

% subplot(224)
% p = plot(t,P_Tx_control_2,'g');
% p(1).LineWidth = 2;
% grid on;
% title('Power on TX before regulation, non rotation');
% xlabel('time(s)');
% ylabel('Pot Tx(W)');

% p = plot(t,Pe_plot,t,Pe_menor_erro*ones(length(t),1),t,Pe_maior_erro*ones(length(t),1));
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% title('BER after regulation on TX, non rotation');
% ylabel('P_e');
% xlabel('time(s)');
% ylim([Pe_menor_erro 1.5*Pe_maior_erro]);
% grid on;
% legend('P_e','P_e lowest error','P_e biggest error')

% figure(2)
% subplot(121)
% p = plot(t,HLOS+HNLOS,'k');
% p(1).LineWidth = 2;
% grid on;
% title('Channel response with time');
% xlabel('time(s)');
% ylabel('HLOS+HNLOS');
%
% subplot(122)
% p = plot(t,P_Tx_control_2,'g');
% p(1).LineWidth = 2;
% grid on;
% title('Power on TX before regulation, non rotation');
% xlabel('time(s)');
% ylabel('Pot Tx(W)');

% figure(2)
% p = plot(t,P_Tx_control_2'-P_Tx_control);
% p(1).LineWidth = 2;
% grid on;
% title('Difference between Powers');
% xlabel('time(s)');
% ylabel('Pot Tx(W)');
% y4 = P_Tx_control_2'-P_Tx_control;
% t4 = t;
% % plot(t,clicks_up,t,clicks_down)
% plot(t,clicks,'c');
% grid on;
% xlabel('time(s)');
% ylabel('clicks');
% % legend('Clicks aumentando','Clicks diminuindo')