clc;
clear;
%% Modulação e Demodulação OOK
%% Entrada

% Quantidade de bits na entrada
N = 10; 

% Sequência aleatória de bits
n = zeros(N,1);
for i=1:N
    n(i) = randi([0 1],1,1);
end

%% Codificação
% Utilizando codificação NRZ
S = 1000;
t = 0:1/S:N;
ind = 1;

% 
m = zeros(length(t),1);
for j=1:length(t)
    if t(j) <= ind
        m(j) = n(ind);
    else
        ind = ind+1;
    end
end
subplot(411)
plot(t,m,'m')
xlabel('tempo');
ylabel('Amplitude');
title('sinal codificado com NRZ');

%% Modulação
% Gerar a portadora
c = cos(2*pi*100*t);
subplot(412);
plot(t,c,'r');
xlabel('tempo');
ylabel('Amplitude');
title('Portadora do sinal');

% Sinal OOK
x = m'.*c;
subplot(413);
plot(t,x,'k');
xlabel('tempo');
ylabel('Amplitude');
title('Sinal OOK');

%% AWGN
EbN0 = 5;
EbN0dB = 10^(EbN0/10);
sigma = sqrt(1/(2*EbN0));
r = n + sigma*randn(N,1);
%% Detecção Coerente
y = x+r;
subplot(414);
y1 = y.*c;
plot(t,y1,'k');
xlabel('tempo');
ylabel('Amplitude');
title('Sinal Demodulado');

% Integrador
int_op = [];
for ii=0:S:length(y1)-S
    int_o = (1/S)*trapz(y1(ii+1:ii+S));
    int_op = [int_op int_o];
end

% Decisor
% Limiar para decisão
Th = 0.5;
disp('Detected bits:')
det = (round(int_op,1)>=Th)

% Calculo da BER
ber = sum(n~=det')/N   