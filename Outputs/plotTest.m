% Julie Read, 2015
% This script reads in the results file from the MCPI CUDA code
% and plots the position, velocity, and hamiltonian (energy check)

close all; clc; clear all;

global omega Deg GM Re C S

% Load GM, Re, C, S
load('aeroegm2008.mat'); 

%GM = GM*1e-9;
%Re = Re/1000;
                            
Deg =200;
omega =0%7292115.0e-011;

fileID = fopen('resultMCPISHMSerialHEO2.txt','r');
% h1 = fscanf(fileID,'%e',[7 100])';
h1 = fscanf(fileID,'%e',[7 Inf])';

time = h1(:,7);
h1(:,1:6) = h1(:,1:6); % Convert to meters, meters/s^2

fclose(fileID);

%simpleEarth

% Plot position
figure
plot3(h1(:,1),h1(:,2),h1(:,3)) % position

set(findobj('type','axes'),'FontSize',20)
title('Position','FontSize',20)
xlabel('m','FontSize',14)
ylabel('m','FontSize',14)
zlabel('m','FontSize',14)
grid on
saveas(gcf,'orbitMCPISHMSerialHEO2.png')

figure
subplot(3,1,1)
plot(time,h1(:,1)) % position
subplot(3,1,2)
plot(time,h1(:,2)) % position
subplot(3,1,3)
plot(time,h1(:,3)) % position
title('Position','FontSize',20)
xlabel('m','FontSize',14)
ylabel('T','FontSize',14)
grid on

% Plot velocity
figure
plot3(h1(:,4),h1(:,5),h1(:,6),'g') % velocity
%axis equal
%saveas(gcf,'orbitMCPIFEMRedSerial2GEO.png')
set(findobj('type','axes'),'FontSize',20)
title('Velocity','FontSize',20)
xlabel('m/s','FontSize',14)
ylabel('m/s','FontSize',14)
zlabel('m/s','FontSize',14)
grid on

disp('Relative error between starting norm(r0) and final segment norm(rf)')
r0 = h1(1,1:3);
v0 = h1(1,4:6);
Rel_Err_Pos = (norm(r0) - norm(h1(end,1:3)))/norm(r0)
disp('Relative error between s5400tarting norm(v0) and final segment norm(vf)')
Rel_Err_Vel = (norm(v0) - norm(h1(end,4:6)))/norm(v0)


% Convert ECI frame estimates to ECEF for energy check - uses flipped result
% since tau uses a +cos term
th_interp = omega*time;   % Increases from 0 to th_max
xECEF(:,1) =  (cos(th_interp)).*h1(:,1) + (sin(th_interp)).*h1(:,2);
xECEF(:,2) = -(sin(th_interp)).*h1(:,1) + (cos(th_interp)).*h1(:,2);
xECEF(:,3) =  h1(:,3);
xECEF(:,4) =  h1(:,4) + omega*h1(:,2);
xECEF(:,5) =  h1(:,5) - omega*h1(:,1);
xECEF(:,6) =  h1(:,6);

%  Plot energy
Jacobi_Integral(time(1:1:end), xECEF(1:1:end,:),'')
saveas(gcf,'hamMCPISHMSerialHEO3.png')