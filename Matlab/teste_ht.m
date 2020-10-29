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
N = max([Nx Ny Nz]);
% base para os Rx
N_Tx = [0 0 -1];
N_Rx = -N_Tx;
%% resposta ao impulso do canal
% em ns (1e-9)
delta_t = 0.5;
t_vector = 0:1000;
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

N = 10;
n = randi([0 1],1,N);

% Função periódica
bits = 1;
P = bits;
htilde = ht'*ones(1,P);
htilde = htilde(:)';

Rb = 1e6; %bps   
Tb = 1/Rb; % tempo de um bit   

Ib = 200e-6;
q = 1.6e-19; % C
N0 = 2*q*Ib;
% Responsividade no Tx
R=1;

% Energia do bit
EbN0 = 130; 
SNR = 10^(EbN0/10);
Eb = N0*SNR;
Ep = 2*Eb;
% P_avg = sqrt(N0*Rb*SNR/(2*R^2));
P_avg = 20e-3; % W
% Pico de Fotocorrente ip = sqrt(Eb) = 2*R*Pr*sqrt(Tb);
ip = 2*R*P_avg;

%% Bits para simulação
pt = ones(1,length(t_vector))*ip;   

%% Gráfico resposta ao impulso
subplot(311)
plot(t_vector,ht)
xlabel("ns");
ylabel("h(t)");
title("Resposta ao impulso");
grid on;

%% Entrada no filtro casado
% Passagem no canal, sinal x(t)
xt = conv(ht,pt);
% Ruído, n(t)
sgma = sqrt(N0*Ep/2);
nt = sgma*randn(1,length(xt))*Rb;
% Plot temporário para o ruído
% subplot(413)
% plot(nt);
% xlabel("ns");
% ylabel("n(t)");
% title("Ruído AGWN");
% grid on;
    
% Soma ruído e sinal recebido
% i(t) = x(t) + n(t)
it = xt + nt;
subplot(312)
plot(it);
xlabel("ns");
ylabel("i(t)");
title("Entrada do filtro Rx");
grid on;

%% Saida no filtro casado
% Pulso retangular para a função kx(T-t), onde T é o tempo de amostragem do
% bit, ou seja,  r(t) = 1/sqrt(Tb)*p(Tb-t). 
rt = pt/sqrt(Tb);
% Convolução 
yt = conv(it,rt);
subplot(313)
plot(yt);
xlabel("ns");
ylabel("y(t)");
title("Saida do filtro Rx");
grid on;

%% Decisor
% Verificar qual limiar deveremos usar para identificação
th = 0.5*Ep;
if yt(length(t_vector)) > th
    bit_out = 1
else
    bit_out = 0
end