clear;
clc;
% Potência média do sinal OOK
p_avg = 1;
% Sensitividade do fotodetector
R = .54; % A/W
% Valor da taxa de bit normalizada
Rb_n = 1;
% duração de um bit
Tb = 1/Rb_n;
% Resolução espectral
df = Rb_n/100;
% vetor de frequências
f=0:df:5*Rb_n;
% frequências normalizada
x = f*Tb;
% Valor que iremos normalizar o valor da PSD
% a = 2*R*p_avg;
% Para a PSD
p = (sinc(x)).^2;
% função delta em 0
p(1) = p(1)+(sinc(0)^2)*(1/Tb);
% Gráfico PSD normalizada por f/Rb
plot(f,p,'r');
% Normalização: energia por bit
xlabel("f/R_{b}");