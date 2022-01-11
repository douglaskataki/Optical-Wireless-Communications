clc;
clear;
close all;
%% Propriedades no Detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potencia LED
PLed = 100e-3; % mW
% campo de visao do fotodetector
FOV = 60; % graus
% area do fotodetector
Area = 5.8E-6; % m^2
% semi-angulo do Tx para meia potencia
theta = 60; % graus
% n Ordem dem Emissao Lambertiana
n = -log10(2)/log10(cosd(theta));
% Ganho optico no concentrador
index = 1;
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
x_Tx = 3*xl/4; y_Tx = yl/3; z_Tx = zl;
shift = yl/4;

Pos_Tx = [x_Tx y_Tx-shift z_Tx; x_Tx 2*y_Tx-shift z_Tx]; %(x_Tx - sqrt(3)*yl/6) yl/2-shift z_Tx];
%Pos_Tx = [x_Tx y_Tx-shift z_Tx];
v = 10e-2;
radius = 1;

n_rot = 3; % numero de rotações

L = 2*pi*radius*n_rot;
t = 0:floor(L/v); %segundos

%fixos
x_Rx = xl/2; y_Rx = yl/2;  z_Rx = 0;
Pos_Rx = [x_Rx y_Rx z_Rx];
N_Rx = [0 0 1];

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
%% Condições inicial
% considerando que estamos utilizando OOK com Ps0 = 0 e P_Rx = Ps1
P_Tx_min = 100e-3;
P_Tx_max = 5000e-3; %W

% valores para SNR
Pe_menor_erro = 1e-6;
Pe_maior_erro = 1e-4;

% potência relativas aos erros
P_Rx_maior_erro = sqrt((N0*(qfuncinv(Pe_maior_erro))^2)/2)*(1/R);
P_Rx_menor_erro = sqrt((N0*(qfuncinv(Pe_menor_erro))^2)/2)*(1/R);

% Alocação de memória para melhorar a velocidade do programa
num_led = ndims(Pos_Tx);

P_Rx = zeros(1,length(t));
Pe_plot = zeros(1,length(t));
HLOS = zeros(num_led,length(t));
HNLOS = zeros(num_led,length(t));
P_Tx = zeros(num_led,length(t));
P_Rx_antes = zeros(1,length(t));
P_Tx_control = zeros(num_led,length(t));
Pe = zeros(num_led,length(t));

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

% alpha
% alpha = linspace(0.8,0.2,length(t));
ang = linspace(0,2*pi*n_rot,length(t));

%lembrar que depende da atualização de potência
%P_Tx(:,1) = [PLed; PLed; PLed];
P_Tx(:,1) = [PLed PLed];
shiftP = 100e-3;
count_1 = 0;
count_2 = 0;
count_3 = 0;

% Referencial
%x<- y-> z^|
%1(0,0)->2(0,yl)
%  |         |
%4(xl,0)<-3(xl,yl)

ledS = [0 0 zl;0 yl zl;xl yl zl;xl 0 zl];
Pot_Rx_Sensor = zeros(4,1);
for i=1:length(t)
    % memoria de potência
    if i==1
        P_Tx_control(:,i) = P_Tx(:,i);
    else
        %recebe a potência anterior
        P_Tx_control(:,i) = P_Tx(:,i-1);
    end
    % calculo da posição estimada do Rx
    for k=1:4
        vect_dist = Pos_Rx-ledS(k,:);
        dist = sqrt(sum(vect_dist).^2);
        N_Tx = [0 0 -1];
        costheta = dot(vect_dist,N_Tx);
        Pot_Rx_Sensor(k) = (n+1)*Area*costheta^n/(2*pi*dist^2)*1e-3;
    end
        Pos_Rx_est = [sum(Pot_Rx_Sensor.*ledS(:,1)) sum(Pot_Rx_Sensor.*ledS(:,2)) 0]/sum(Pot_Rx_Sensor);
    for j=1:num_led
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% HLOS-beginning %%%%%%
        vTxRx = Pos_Tx(j,:)-Pos_Rx;
        d = sqrt(sum((vTxRx).^2));
        N_Tx = -vTxRx/d;
        %N_Tx = N_Tx_0;
        % ângulo entre Tx e Rx
        cosphi = dot(vTxRx,N_Rx)/d;
        phi_los = acosd(cosphi);
        % fator para saber o novo phi_los rotacionado
        
        if(phi_los<FOV)
            costheta = abs(dot(vTxRx,N_Tx))/d;
            HLOS(j,i) = (n+1)*Area*costheta^n/(2*pi*d^2);
        else
            HLOS(j,i) = 0;
        end
        %%%%%%% HLOS-ending %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% HNLOS-beginning %%%%%%%%%%
        
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
                vTxWP = Pos_Tx(j,:) - WP;
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
                vTxWP = Pos_Tx(j,:) - WP;
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
                vTxWP = Pos_Tx(j,:) - WP;
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
                vTxWP = Pos_Tx(j,:) - WP;
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
        HNLOS(j,i) = h1+h2+h3+h4;
        %%%%%%% HNLOS-ending %%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%% calcula valor de P_Rx para que chega no local
    for j=1:num_led
        P_Rx_antes(i) = P_Rx_antes(i) +  P_Tx_control(j,i)*(HLOS(j,i)+HNLOS(j,i));
    end
   
    % considerando OOK, ajusta potência se for depender dos valores da SNR
    Eb_antes = 2*(P_Rx_antes(i)*R)^2;
    Pe(i) = qfunc(sqrt(Eb_antes/N0));
    y_Rx = xl/4+radius*sin(ang(i)); % vx metros em um segundo
    x_Rx = yl/4+radius*cos(ang(i));
    Pos_Rx = [x_Rx y_Rx z_Rx];
    
    %% Controle de acordo com o Pe
    for k = 1:num_led
        if(Pe(i) >= Pe_maior_erro)
            P_Tx(k,i) = (P_Rx_maior_erro/num_led)/(HLOS(k,i)+HNLOS(k,i));
            
            %verifica saturação
            if(P_Tx(k,i) >= P_Tx_max)
                P_Tx(k,i) = P_Tx_max;
            end
            count_1 = count_1+1;
        elseif(Pe(i) <= Pe_menor_erro)
            P_Tx(k,i) = (P_Rx_menor_erro/num_led)/(HLOS(k,i)+HNLOS(k,i));
            
            %verifica saturação
            if(P_Tx(k,i) <= P_Tx_min)
                P_Tx(k,i) = P_Tx_min;
            end
            count_2 = count_2+1;
        else
            %manter a potencia do P_Tx constante ou aumentar ela de um
            %shift
            P_Tx(k,i) = P_Tx_control(k,i);
            count_3 = count_3+1;
        end
        
        P_Rx(i) = P_Rx(i) + P_Tx(k,i)*(HLOS(k,i)+HNLOS(k,i));      
    end
    
    %atualiza o valor de Pe
    Eb = 2*(P_Rx(i)*R)^2;
    Pe_plot(i) = qfunc(sqrt(Eb/N0));
    
%     % atualiza a posição de Rx
%     y_Rx = xl/4+radius*sin(ang(i)); % vx metros em um segundo
%     x_Rx = yl/4+radius*cos(ang(i));
%     Pos_Rx = [x_Rx y_Rx z_Rx];
    
end


figure(1)
% p = plot(t,P_Tx(1,:),t,P_Tx(2,:),t,P_Tx(3,:),t,P_Tx(1,:)+P_Tx(2,:)+P_Tx(3,:));
% p(1).LineWidth = 2;
% p(2).LineWidth = 2;
% p(3).LineWidth = 2;
% p(4).LineWidth = 2;

p =  plot(t,P_Tx(1,:),t,P_Tx(2,:),t,P_Tx(1,:)+P_Tx(2,:));
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
grid on;
title('a) Transmitted power');
xlabel('time(s)');
ylabel('P_{Tx}(W)');
legend('P_{TX1}','P_{TX1}','sumP_{TX}');
%,'P_{TX3}','sumP_{TX}');


figure(2)
subplot(211)
p = plot(t,P_Rx,t,P_Rx_menor_erro*ones(1,length(t)),t,P_Rx_maior_erro*ones(1,length(t)));
p(1).LineWidth = 2;
p(2).LineWidth = 2;
p(3).LineWidth = 2;
grid on;
title('b)Received Power');
xlabel('time(s)');
ylabel('P_{Rx}(W)');
legend('P_{Rx}','P_{low}','P_{high}');

subplot(212)
% p = plot(t,Pe_plot,t,Pe_menor_erro*ones(length(t),1),t,Pe_maior_erro*ones(length(t),1));
plot(t,log10(Pe_plot));
title('c) bit error rate');
ylabel('log(BER)');
xlabel('time(s)');
% ylim([Pe_menor_erro 1.5*Pe_maior_erro]);
grid on;
% legend('P_e','P_{el}','P_{eh}')