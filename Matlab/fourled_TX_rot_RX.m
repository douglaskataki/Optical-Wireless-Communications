clear;
clc;
%% Dados iniciais
% potencia led
PLed = 20; % mW
% campo de visão do fotodetector
FOV = 30; % graus
% area do fotodetector
Area = 5.8E-6; % m^2
% semi angulo do Tx para meia potência
theta_meio = 30; % graus
%% Posições
% tamanhos iniciais no local
xl = 5; yl = 5; zl = 3; 

% altura do Tx
h = 2.15;

% Posicao Tx
x_Tx = xl/4; 
y_Tx = yl/4; 
z_Tx = h;

% Posicao Rx
x_Rx = xl/4; 
y_Rx = yl/4; 
z_Rx = 0;

% vetores
Pos_Tx = [x_Tx y_Tx z_Tx];
Pos_Rx = [x_Rx y_Rx z_Rx];
dim_sala = [xl yl zl];

%% Normal
N_Rx = [0 0 -1].';
N_Rx = -N_Rx;
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

%% Demais dados
% ordem do lambertiano
n = -log(2)/log10(cosd(theta_meio));
% vetor distancia entre Tx e Rx
vTxRx = Pos_Tx - Pos_Rx;
% distancia entre Tx e Rx
d = sqrt(sum((vTxRx).^2));
% angulo de incidencia direto
cosphi = abs(dot(vTxRx,N_Rx))/d;
phi_los = acosd(cosphi);
% ganho no concentrador, considerando que ele tem uma lente com indice de
% refração index
index = 1.5;
G_Con = index^2/sind(FOV)^2;
Ts = 1;

% coeficiente de reflexão
rho = 0.8;

% divisão para o grid da parede
M = 10;
Nx = xl*M; Ny = yl*M; Nz = round(zl*M);
x=linspace(0,xl,Nx);
y=linspace(0,yl,Ny);    
z=linspace(0,zl,Nz);

% normaliza o locais para rotações em torno dos eixos x e y

x_theta = linspace(-90,90,Nx);
y_theta = linspace(-90,90,Ny);

% alocação de memória para as funções de transferência
N = max([Nx Ny Nz]);
h = zeros(N);
h1 = zeros(N);
h2 = zeros(N);
h3 = zeros(N);
h4 = zeros(N);
%% Cálculo do HLOS +  HNLOS
for nled=1:4
    for ii=1:Nx
        rotx = [1 0 0; 0 cosd(x_theta(ii)) -sind(x_theta(ii)); 0 sind(x_theta(ii)) cosd(x_theta(ii))];
        N_Rx = rotx*N_Rx;
        N_Rx1 = N_Rx;
        for jj=1:Ny
            roty = [cosd(y_theta(jj)) 0 sind(y_theta(jj)); 0 1 0; -sind(y_theta(jj)) 0 cosd(y_theta(jj))];
            N_Rx = roty*N_Rx;
            % HLOS
            if(phi_los>=0)&&(phi_los<=FOV)
                costheta = abs(dot(vTxRx,N_Rx))/d;
                h(ii,jj) = (n+1)*Ts*Area*costheta^n/(2*pi*d^2);
            else
                h(ii,jj) = 0;
            end
            % reflexão nas paredes (HNLOS)
            h1(ii,jj)=0;
            h2(ii,jj)=0;
            h3(ii,jj)=0;
            h4(ii,jj)=0;

            %%%%%%%%%%%%%%% wall (Z,Y) %%%%%%%%%%%%%%%%%%
            % fator de área
            dA = zl*yl/(Ny*Nz);
            % normal da parede
            n1 = [1 0 0];
            for kk=1:Ny
                for ll=1:Nz
                    % ponto na parede Z,Y - Wall Point
                    WP = [0 y(kk) z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    % distancia do TX para a parede Z,Y(1)
                    D1 = sqrt(sum((vTxWP).^2));
                    % distancia TX e Rx
                    % angulos Tx e incidencia na parede
                    cos_phi = abs(dot(N_Rx,vTxWP))/D1;
                    cos_alpha = abs(dot(vTxWP,n1))/D1;
                    % distance from WP to receiver
                    vRxWP = Pos_Rx-WP;
                    D2 = sqrt(sum((vRxWP).^2));
                    % angulos Rx e reflexão
                    cos_psi = abs(dot(vRxWP,N_Rx))/D2;
                    cos_beta = abs(dot(vRxWP,n1))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h1(ii,jj)=h1(ii,jj)+(n+1)*Area*rho*dA*...
                        cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end
                end
            end

            %%%%%%%%%%%%%%% Parede (Z,X) %%%%%%%%%%%%%%%%%%
            % fator de área
            dA = zl*xl/(Nx*Nz);
            % normal da parede
            n2 = [0 1 0];
            for kk=1:Ny
                for ll=1:Nz
                    % ponto na parede Z,X - Wall Point
                    WP = [x(kk) 0 z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    % distancia do TX para a parede Z,
                    D1 = sqrt(sum((vTxWP).^2));
                    % distancia TX e Rx
                    % angulos Tx e incidencia na parede
                    cos_phi = abs(dot(N_Rx,vTxWP))/D1;
                    cos_alpha = abs(dot(vTxWP,n2))/D1;
                    % distance from WP to receiver
                    vRxWP = Pos_Rx-WP;
                    D2 = sqrt(sum((vRxWP).^2));
                    % angulos Rx e reflexão
                    cos_psi = abs(dot(vRxWP,N_Rx))/D2;
                    cos_beta = abs(dot(vRxWP,n2))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h2(ii,jj)=h2(ii,jj)+(n+1)*Area*rho*dA*...
                        cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end
                end
            end

            %%%%%%%%%%%%%%% wall (Z,Y) %%%%%%%%%%%%%%%%%%
            % fator de área
            dA = zl*yl/(Ny*Nz);
            % normal da parede
            n3 = [-1 0 0];
            for kk=1:Ny
                for ll=1:Nz
                    % ponto na parede Z,Y(1) - Wall Point
                    WP = [xl y(kk) z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    % distancia do TX para a parede Z,Y(1)
                    D1 = sqrt(sum((vTxWP).^2));
                    % distancia TX e Rx
                    % angulos Tx e incidencia na parede
                    cos_phi = abs(dot(N_Rx,vTxWP))/D1;
                    cos_alpha = abs(dot(vTxWP,n3))/D1;
                    % distance from WP to receiver
                    vRxWP = Pos_Rx-WP;
                    D2 = sqrt(sum((vRxWP).^2));
                    % angulos Rx e reflexão
                    cos_psi = abs(dot(vRxWP,N_Rx))/D2;
                    cos_beta = abs(dot(vRxWP,n3))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h3(ii,jj)=h3(ii,jj)+(n+1)*Area*rho*dA*...
                        cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end
                end
            end

            %%%%%%%%%%%%%%% Parede (Z,X) %%%%%%%%%%%%%%%%%%
            % fator de área
            dA = zl*xl/(Nx*Nz);
            % normal da parede
            n4 = [0 -1 0];
            for kk=1:Ny
                for ll=1:Nz
                    % ponto na parede Z,X - Wall Point
                    WP = [x(kk) yl z(ll)];
                    vTxWP = Pos_Tx - WP; 
                    % distancia do TX para a parede Z,
                    D1 = sqrt(sum((vTxWP).^2));
                    % distancia TX e Rx
                    % angulos Tx e incidencia na parede
                    cos_phi = abs(dot(N_Rx,vTxWP))/D1;
                    cos_alpha = abs(dot(vTxWP,n4))/D1;
                    % distance from WP to receiver
                    vRxWP = Pos_Rx-WP;
                    D2 = sqrt(sum((vRxWP).^2));
                    % angulos Rx e reflexão
                    cos_psi = abs(dot(vRxWP,N_Rx))/D2;
                    cos_beta = abs(dot(vRxWP,n4))/D2;
                    if abs(acosd(cos_psi))<=FOV
                        h4(ii,jj)=h4(ii,jj)+(n+1)*Area*rho*dA*...
                        cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                    end
                end
            end
            N_Rx = N_Rx1;
        end
        N_Rx = [0 0 -1].';
    end
    if(nled==2) 
        Pos_Tx = [xl/2 yl/4 zl];
    end
    if(nled==3)
        Pos_Tx = [xl/4 yl/2 zl];
    end
    if(nled==4)
        Pos_Tx = [xl/2 yl/2 zl];
    end
end

P_Rx = (h1+h2+h3+h4+h)*Ts*G_Con*PLed;
P_Rx_dBm = 10*log10(P_Rx);

%% Plot do diagrama de potência em dBm
figure(1)
meshc(x_theta,y_theta,P_Rx_dBm);
% axis([min(x_theta) max(x_theta) min(y_theta) max(y_theta) -200 max(max(P_Rx_dBm_rot))]);
% %meshc(x,y,SNR_dB);
grid on;
title('LOS+NLOS+Rot_Rx(x,y) - Potência');
xlabel('\alpha(graus)');
ylabel('\beta(graus)');
% % zlabel('Potência(dBm)')
zlabel('Potência(dBm)');
colorbar;
colormap jet;

%% SNR
% ruído
% Dados retirados Fundamental Analysis for Visible-Light Communication System using LED Lights
q = 1.6E-19; % Carga do eletron
k = physconst('Boltzmann');
c = physconst('LightSpeed');
% banda (100Mb/s)
B = 100E6; % Mbps
% Responsividade do Fotodetector
R = 0.54/1e-3; % A/W -> correção para mW
% photocorrent due to ground radiation
I_b = 5100E-6; % uA
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
noise2 = sigma2_shot+sigma2_thermal;
SNR = R^2*P_Rx.^2./noise2;
SNR_dB = 10*log10(SNR);

%% Plot para o SNR
figure(2)
meshc(x_theta,y_theta,SNR_dB);
% axis([0 max(x) 0 max(y) min(min(P_Rx_dBm)) max(max(P_Rx_dBm))]);
grid on;
title('LOS+NLOS+Rot_Rx(x,y)-SNR');
xlabel('\alpha(graus)');
ylabel('\beta(graus)');
zlabel('SNR(dB)');
colorbar;
colormap jet;