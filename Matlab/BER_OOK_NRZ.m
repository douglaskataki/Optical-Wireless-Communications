%% Simulação da BER da modulação OOK-NRZ
clear;
clc;
%% Parâmtros iniciais
% Carga do eletron
q = 1.6E-19; 
% Corrente Background + interferência
Ib = 200e-6;
% Densidade espectral de ruído, 2*q*Ib
N0 = 2*q*Ib;
% Responsibidade do Fotodetector TX;
R_Tx = 1;
% Taxa de bit
Rb = 1e6; % 1Mbps
% Duração de um bit
Tb = 1/Rb;
% Número de bits
sig_length = 1e5;
% Número de amostras por símbolo
nsamp = 10;
% Tempo de amostragem
Tsamp = Tb/nsamp;

% Razão sinal ruído em dB
EbN0 = 1:12;
% Razão EbN0
SNR = 10.^(EbN0./10);
% Alocação de memória para algumas variáveis
P_avg = zeros(1,length(SNR));
i_peak = zeros(1,length(SNR));
Ep = zeros(1,length(SNR));
sgma = zeros(1,length(SNR));
ber = zeros(1,length(SNR));

% ------ Simulação da probabilida de erros ----- %
for i=1:length(SNR)
    % Gerando sinal OOK randômico
    OOK = randi([0 1],1,sig_length);
    
    % Potência média óptica transmitida (Pr)
    P_avg(i) = sqrt(N0*Rb*SNR(i)/(2*R_Tx^2));
    % Pico de amplitude elétrica
    i_peak(i)=2*R_Tx*P_avg(i);
    % Pico de energia (Energia por bit é Ep/2)
    Ep(i) = i_peak(i)^2*Tb;
    
    %% Filtro de transmissão p(t)
    pt = ones(1,nsamp)*i_peak(i);
    
    %% Formação do Pulso retangular 
    Tx_signal = rectpulse(OOK,nsamp)*i_peak(i);
  
    %% Entrada no filtro casado (y = x+n)
    % Variancia do ruído
    sgma(i) = sqrt((N0/2)/Tsamp);
    Rx_signal=R_Tx*Tx_signal+sgma(i)*randn(1,length(Tx_signal));
    % Filtro RX casado com o filtro de TX
    rt = pt;
    % Saída do filtro casado
    MF_out = conv(Rx_signal,rt)*Tsamp;
    
    %% Amostragem
    MF_out_downsamp = MF_out(nsamp:nsamp:end);
    % Truncagem
    MF_out_downsamp = MF_out_downsamp(1:sig_length);
    
    %% Decisor
    Rx_th = zeros(1,sig_length);
    % Se for maior que Ep/2, então é o bit 1, caso contrário, será o bit 0.
    Rx_th(find(MF_out_downsamp>Ep(i)/2))=1;
    
    % Cálculo da BER
    [nerr ber(i)] = biterr(OOK,Rx_th);

end

semilogy(EbN0,ber,'b');
hold on;
semilogy(EbN0,qfunc(sqrt(10.^(EbN0/10))),'r-X','linewidth',2);
grid on;
legend('Simulação','Teoria');
xlabel('Eb/No, dB');
ylabel('BER');
title('BER para modulação OOK');