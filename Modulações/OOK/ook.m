clc;
clear;
%% Modulação e Demodulação OOK
%% Entrada

% Quantidade de bits na entrada
N = 10; 

% Sequência aleatória de bits
n = randi([0 1],1,N);

%% Codificação
% Utilizando codificação NRZ unipolar
S = 1000; % fs -> frequencia de amostragem
t = 0:1/S:N;

ind = 1;
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
c = cos(2*pi*2*t);
subplot(412);
plot(t,c,'r');
xlabel('tempo');
ylabel('Amplitude');
title('Portadora do sinal');

% Sinal OOK + Portadora
x = m'.*c;
subplot(413);
plot(t,x,'k');
xlabel('tempo');
ylabel('Amplitude');
title('Sinal OOK');

%% Detecção Coerente
y = x;
subplot(414);
y1 = y.*c;
plot(t,y1,'k');
xlabel('tempo');
ylabel('Amplitude');
title('Sinal Demodulado');

% Integrador (retirar a portadora e considerando que a fase de demodulação
% está sincrona com a fase de envio
int_op = [];
for ii=0:S:length(y1)-S
    int_o = (1/S)*trapz(y1(ii+1:ii+S));
    int_op = [int_op int_o];
end

% Decisor
% Limiar para decisão
Th = .5;
disp('Bits detectados:')
det = (round(int_op,1)>=Th)

% Calculo da BER
ber = sum(n~=det)/N   