clc;
clear;
%% Modulação e Demodulação OOK(banda base) + AWGN
%% Entrada

% Quantidade de bits na entrada
N = 1000000; 

% Sequência aleatória de bits
x = zeros(N,1);
for i=1:N
    x(i) = randi([0 1],1,1);
end

ber_sim = [];
ber_theory_erfc=[];

for EbN0dB = 0:1:15
    EbN0 = 10^(EbN0dB/10);
    sigma = sqrt(1/(2*EbN0));
    
    % Canal AWGN 
    r = x' + sigma*randn(1,N);
    % Decisor
    m_cap = (r>0.5);
    
    % Cálculo da BER (number of error)
    ber_sim1 = sum(m_cap~=x')/N;
    ber_sim = [ber_sim ber_sim1];
    ber_th_er = 0.5*erfc(sqrt(EbN0/4));
    ber_theory_erfc = [ber_theory_erfc ber_th_er];
end

EbN0dB = 0:1:15;
semilogy(EbN0dB,ber_sim,'r*-',EbN0dB,ber_theory_erfc,'ko-');
xlabel('Eb/N_0 (dB)');
ylabel('BER');
legend('Simulated','Theoretical(erfc)');

grid on;