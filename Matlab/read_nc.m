clear;
clc;
data = 'ssi_v02r01_daily_s18820101_e18821231_c20170717.nc';
%ncinfo(data);
%ncdisp(data);
SSI = ncread(data,'SSI');
lambda = ncread(data,'wavelength');

i = find(lambda == 5225);
figure(1)
plot(SSI(i,:));
xlabel('day');
ylabel('W/m^2nm');
grid on;

% figure(2)
% [x,y] = meshgrid(1:365,lambda);
% surf(x,y,SSI);
% xlabel('day');
% ylabel('nm');
% zlabel('W/m^2nm');
% grid on;