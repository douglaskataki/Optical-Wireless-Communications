clear;
clc;
% Pot�ncia m�dia do sinal OOK
p_avg = 1;
% Sensitividade do fotodetector
R = .54; % A/W
% Valor da taxa de bit normalizada
Rb_n = 1;
% dura��o de um bit
Tb = 1/Rb_n;
% Resolu��o espectral
df = Rb_n/100;
% vetor de frequ�ncias
f=0:df:5*Rb_n;
% frequ�ncias normalizada
x = f*Tb;
% Valor que iremos normalizar o valor da PSD
% a = 2*R*p_avg;
% Para a PSD
p = (sinc(x)).^2;
% fun��o delta em 0
p(1) = p(1)+(sinc(0)^2)*(1/Tb);
% Gr�fico PSD normalizada por f/Rb
plot(f,p,'r');
% Normaliza��o: energia por bit
xlabel("f/R_{b}");