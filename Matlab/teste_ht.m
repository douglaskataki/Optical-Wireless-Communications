clc;
clear;  
%% Propriedades no Detector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potencia LED
PLed = 20e-3; % W
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
% Posições Tx e Rx
x_Tx = xl/2; y_Tx = yl/2; z_Tx = zl;
Pos_Tx = [x_Tx y_Tx z_Tx];
x_Rx = xl/5; y_Rx = yl/5; z_Rx = 0;
Pos_Rx = [x_Rx y_Rx z_Rx];
% normaliza o locais para os Rx
M = 10;
Nx = xl*M; Ny = yl*M; Nz = round(zl*M);
x = linspace(0,xl,Nx);
y = linspace(0,yl,Ny);    
z = linspace(0,zl,Nz);
% N = max([Nx Ny Nz]);
% base para os Rx
N_Tx = [0 0 -1];
N_Rx = -N_Tx;
%% Taxa de transmissão
Rb = 8e3; % 8kbps

%% resposta ao impulso do canal
% em ns (1e-9)
delta_t = Rb;
tf = 1/Rb*1e9;
t_vector = linspace(0,tf,delta_t);
ht = zeros(1, length(t_vector));

%% Tx e Rx (LOS)
D0 = sqrt(sum(Pos_Tx-Pos_Rx).^2);
cosphi = z_Tx/D0;
tau = D0/c;
ind = find(round(tau/delta_t)==t_vector);
if abs(acosd(cosphi))<= FOV
    ht(ind) = ht(ind) + (n+1)*Area*cosphi^(n+1)/(2*pi*D0^2);
end

%% Tx e Rx (NLOS)
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
        % tempo para percorrer a primeira reflexão
        tau1 = (D1+D2)/c;
        ind = find(round(tau1/delta_t)==t_vector);
        if abs(acosd(cos_psi))<=FOV
           ht(ind) = ht(ind)+(n+1)*Area*rho*dA*...
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
        % tempo
        tau2 = (D1+D2)/c;
        ind = find(round(tau2/delta_t)==t_vector);
        if abs(acosd(cos_psi))<=FOV
           ht(ind) = ht(ind)+(n+1)*Area*rho*dA*...
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
        % tempo
        tau3 = (D1+D2)/c;
        ind = find(round(tau3/delta_t)==t_vector);
        if abs(acosd(cos_psi))<=FOV
           ht(ind) = ht(ind)+(n+1)*Area*rho*dA*...
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
        % tempo
        tau4 = (D1+D2)/c;
        ind = find(round(tau4/delta_t)==t_vector);
        if abs(acosd(cos_psi))<=FOV
           ht(ind) = ht(ind)+(n+1)*Area*rho*dA*...
           cos_phi^n*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
        end
    end
end

%% Pot
% HLOS
H = 0;


vTxRx = Pos_Tx-Pos_Rx;
d = sqrt(sum((vTxRx).^2));
% ângulo entre Tx e Rx
cosphi = abs(dot(vTxRx,N_Rx))/d;
phi_los = acosd(cosphi);
if(phi_los<=FOV)
    costheta = abs(dot(vTxRx,N_Tx))/d;
    % H = (n+1)*Area*costheta^n/(2*pi*d^2);
    H = HLOS(Area,FOV,n,1,1,1,x_Rx,y_Rx,N_Tx,Pos_Tx,N_Rx);
end

% HNLOS
H = H + HNLOS(xl,yl,zl,rho,Area,FOV,n,1,1,1,x_Rx,y_Rx,z_Rx,N_Tx,Pos_Tx,N_Rx);
% Total
P_Rx = PLed*G_Con.*H;

%% !!!!!!!!!!!!!!!!!!!!! multiplicação de escala "errada" !!!!!!!!!!!!!!!!!!!!!
ht = ht/max(ht);

%% Ruído
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dados retirados Fundamental Analysis for Visible-Light Communication System using LED Lights

q = 1.6E-19; % Carga do eletron
k = physconst('Boltzmann');
c = physconst('LightSpeed');

% banda (100Mb/s)
B = 2*Rb; % igual a taxa de TX

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

% ruído de background (fundo)
pbn = 5.8e-6/1e-4; %uW/cm^2
d_l = 30e-9; %nm
sigma2_background = 2*q*I_2*pbn*Area*d_l*R*B;

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
% lambda = ncread(data,'wavelength');
% i = find(lambda == 5225);

% w = mean(SSI(i,:));

% sigma2_ssi = 2*q*B*Rb*w*d_l;  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ruído total
N0 = sigma2_shot + sigma2_thermal + sigma2_background; %+ sigma2_ssi;
Tb = 1/Rb; % tempo de um bit

% Conversão de Lux para Watts
% P = E*A/eff, eff = 90 para led
% eff = 90;

% Responsividade no Rx
Resp_Rx = 1;

% Energia do bit
EbN0 = 1; %dB
SNR = 10^(EbN0/10);
Eb = N0*SNR;
%Ep = 2*Eb;
P_avg = sqrt(N0*Rb*SNR/(2*Resp_Rx^2));


% Pico de Fotocorrente ip = 2*R*Pr*sqrt(Tb);
ip = 2*Resp_Rx*P_avg*1e2; 
Ep = ip^2*Tb; 

%% Sequencia de bits randomica com tamanho fixo N
N = 1e4;

% bits a serem transmitidos
n = randi([0 1],1,N);
L = length(t_vector);
nf = zeros(1,length(n));

%% Bits para simulação
% Multiplica pela amplitude do bit  
% seq = [];
% soma = 0;
% simulação para cada bit
% Para valor 1 => IM/DD, para valor 2 => PAM-4
type = 1;

%% IM/DD
if type == 1
    % codificação Manchester de n
    N_code = Manchester_code(n);
    N_decode = zeros(1,length(N_code));
    
    for j=1:length(N_code)
    if N_code(j)>0
       %pt = [ones(1,ceil(L/2))*ip zeros(1,ceil(L/2))];
       pt = ones(1,round(L/2))*ip;
    else
       %pt = [zeros(1,ceil(L/2)) ones(1,ceil(L/2))*ip];
       pt = zeros(1,round(L/2));
    end
    
    % Simulação da luz
    % Variação de Potências constantes durante a simulação, entre 1000 e
    % 3000 lux
    % E = randi([1000 3000],1);
    % E = 100;
    % 
    % P_Lux = E*Area/eff;
    % converter em corrente
    % ip1 =  2*R_Tx*P_Lux;
    
    % Soma dos valores 
    % pt1 = pt+ip1;
    
    % Passagem no canal, tendo como saida o sinal x(t)
    xt = conv(ht,pt);
    
    % Ruído, n(t)
    sgma = sqrt(N0/2/Tb);
    nt = sgma*randn(1,length(xt))*Rb;
    it = nt+xt;
    % it = awgn(xt,SNR,10*log10(P_avg*R)); -> precisa do pack
    % communications
    
    % Soma ruído e sinal recebido (considerando ruído aditivo)
    % i(t) = x(t) + n(t)
    %it = xt + nt;
    
    % Saida no filtro casado
    % Pulso retangular para a função kx(T-t), onde T é o tempo de amostragem do
    % bit, ou seja,  r(t) = 1/sqrt(Tb)*p(Tb-t).
    % o fator 2 apareceu devido a condificação Manchester
    
    rt = pt*1/sqrt(2*Tb);
    % Convolução
    
    yt = conv(it,rt);
    
    % Reconhecimento do bit de saida
    
    % th = limiar de decisão
    th = 0.5*Ep;
    if yt(round(L/2)) > th
       seq = 1; 
    else
       seq = 0;
    end
    %soma = soma + abs(seq-N_code(j));
    N_decode(j) = seq;
    % yt = [yt0 yt(1:L)];
    % seq = [seq0 seq]; -> deste modo estava demorando muito para fazer a
    % comparação entre os bits recebidos.
    end
nf = Manchester_decode(N_decode); 

%% PAM-4
Ep = Ep/2;
elseif type == 2
    for j=1:2:length(n)-1
        % Modulação para o 4-PAM, considerando valores de 3,1,-1,-3 para os
        % níveis para transmissão
        if n(j)==1
            if n(j+1)==1 % symbol = 11
               Lvl = ip;
            else % symbol = 10
               Lvl = 2/3*ip;
            end
        else
            if n(j+1)==0 % symbol = 00 
               Lvl = 1/3*ip;
            else % symbol = 01
               Lvl = 0;
            end
        end
        pt = ones(1,round(L/2))*Lvl;
        
        xt = conv(ht,pt);
        
        sgma = sqrt(N0*Ep/2);
        nt = sgma*randn(1,length(xt))*Rb;
        
        it = nt+xt;
        
        rt = pt*1/sqrt(2*Tb);
        % Convolução 
        yt = conv(it,rt);
        
        % Reconhecimento do bit de saida
        % th = limiar de decisão
        if yt(round(L/2)) >= 2/3*Ep
            % bit 11
            nf(j) = 1;
            nf(j+1) = 1;
        elseif yt(round(L/2)) < 2/3*Ep && yt(round(L/2))>= 1/3*Ep
            % bit 10
            nf(j) = 1;
            nf(j+1) = 0;
        elseif yt(round(L/2)) < 1/3*Ep && yt(round(L/2))>= 0*Ep
            % bit 00
            nf(j) = 0;
            nf(j+1) = 0;
        else
            % bit 01
            nf(j) = 0;   
            nf(j+1) = 1;
        end
    end
end

%% bit_error
soma = sum(abs(nf-n))/N
% Mostra que conseguiu identificar a sequência
% bits = abs(seq - n);
% error_bit = sum(bits~=0)/N

