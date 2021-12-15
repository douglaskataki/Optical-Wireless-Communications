clear;
clc;

load('perfisdepot_3_5_7_10.mat');
figure(1)
%subplot(121)
p = plot(t1,y1,t2,y2,t3,y3,t4,y4);
for i=1:4
    p(i).LineWidth = 2;
end
grid on;
%title('a) Power difference vs time');
xlabel('time(s)');
ylabel('\DeltaP_{Tx} (W)');
ylim([0 3.5]);
legend('y_l=3 m','y_l=5 m','y_l=7 m','y_l=10 m','Location','southeast');
% legend('boxoff');

figure(2)
%subplot(122)
x = [0 3 5 7 10];
y = [0 max(y1) max(y2) max(y3) max(y4)];
coef4 = polyfit(x,y,length(y)-1);
coef3 = polyfit(x,y,length(y)-2);
coef2 = polyfit(x,y,length(y)-3);

d = 0:0.01:10;
yf4 = zeros(1,length(d));
yf3 = zeros(1,length(d));
yf2 = zeros(1,length(d));
for i=1:length(d)
    yf4(i) = coef4(1)*d(i)^4+coef4(2)*d(i)^3+coef4(3)*d(i)^2+coef4(4)*d(i)+coef4(5);
    yf3(i) = coef3(1)*d(i)^3+coef3(2)*d(i)^2+coef3(3)*d(i)+coef3(4);
    yf2(i) = coef3(1)*d(i)^3+coef3(2)*d(i)^2+coef3(3);
end

p = plot(d,yf4,'r',x,y,'go');
p(1).LineWidth = 2;
p(2).LineWidth = 2;
%p(3).LineWidth = 2;
grid on;
%title('b) Power difference vs distance');
xlabel('Room length (m)');
ylabel('\DeltaP_{max} (W)');
ylim([0,max(yf3)*1.25])
% legend('y_l=3 m','y_l=5 m','y_l=7 m','y_l=10 m','Location','southeast');
% legend('boxoff');
legend('fitting curve','original data','Location','southeast')