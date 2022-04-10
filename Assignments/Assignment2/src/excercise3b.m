% Assignment 2 - Exercise 3b

clear; clc; close all;
%% Approximation using 2nd-order Central Difference

N = 200;
j = linspace(0,N,N+1); %odd number of grid points
dx = (2*pi - 0)/(N+1);
x = j.*dx;

u = @(x,t) exp(sin(x - 2*pi*t)); %solution 
A = CentralDiff_2nd(N);
du_dx_CD2 = @(u) (A * u')'; %transpose is necessary for RK

numerical_1 = zeros(N+1,3);
analytic = zeros(N+1,3);
n = 0;
for t = [0 100 200]
    n = n+1;

    t0 = 0;
    resolution = 100; %factor to make dt smaller, i.e., to make grid finer
    dt = 1/resolution * (t - t0)/(N+1);

    u0 = u(x,t0); % initial solution at t = 0
    u_n = u0; %initial solution at t = 0
    for iter = dt:dt:t
        u_n = RungeKutta_4(u_n,du_dx_CD2,dt); %advance in time
    end
    tend = t;

    % analytic solution for comparison
    analytic(:,n) = u(x,tend);
    numerical_1(:,n) = u_n;
end


%% Approximation using global Fourier method
% clear; clc;

N = 10;
j = linspace(0,N,N+1); %odd number of grid points
dx = (2*pi - 0)/(N+1);
x2 = j.*dx; 

u = @(x,t) exp(sin(x - 2*pi*t)); %solution 
du_dx_GM = @(u) (D_odd(N) * u')'; %transpose is necessary for RK

numerical_2 = zeros(N+1,3);

n = 0;
for t = [0 100 200]
    n = n+1;

    t0 = 0;
    u0 = u(x2,t0); % initial solution at t = 0
    u_n = u0; %initial solution at t = 0
    for iter = dt:dt:t 
        u_n = RungeKutta_4(u_n,du_dx_GM,dt); 
    end
    tend = iter; %this is the last point of the iteration

    % analytic solution for comparison
    % analytic(:,n) = u(x,tend);
    numerical_2(:,n) = u_n;
end


%% Visualise
fig = figure('position',[400 400 750 1000]);

% top left, t = 0, 2ND ORDER FD
subplot(3,2,1);
plot(x,analytic(:,1),'--',...
    'displayname','Analytic');
hold on;
plot(x,numerical_1(:,1),'displayname','$2^{\mathrm{nd}}$-order FD');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('$u(x,t)$','interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);
legend('interpreter', 'latex','fontsize',14);
title('N = 200, t = 0');
hold on;

% middle left, t = 100, 2ND ORDER FD
subplot(3,2,3);
plot(x,analytic(:,1),'--',...
    'displayname','Analytic');
hold on;
plot(x,numerical_1(:,2),'displayname','$2^{\mathrm{nd}}$-order FD');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('$u(x,t)$','interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);
legend('interpreter', 'latex','fontsize',14);
title('N = 200, t = 100');
hold on;

% bottom left, t = 200, 2ND ORDER FD
subplot(3,2,5);
plot(x,analytic(:,1),'--',...
    'displayname','Analytic');
hold on;
plot(x,numerical_1(:,3),'displayname','$2^{\mathrm{nd}}$-order FD');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('$u(x,t)$','interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);
legend('interpreter', 'latex','fontsize',14);
title('N = 200, t = 200');
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% top right, t = 0, GLOBAL METHOD
subplot(3,2,2)
plot(x,analytic(:,1),'--',...
    'displayname','Analytic');
hold on;
plot(x2,numerical_2(:,1),'displayname','Fourier');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('$u(x,t)$','interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);
legend('interpreter', 'latex','fontsize',14);
title('N = 10, t = 0');
hold on;

% middle right, t = 100, GLOBAL METHOD
subplot(3,2,4)
plot(x,analytic(:,2),'--',...
    'displayname','Analytic');
hold on;
plot(x2,numerical_2(:,2),'displayname','Fourier');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('$u(x,t)$','interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);
legend('interpreter', 'latex','fontsize',14);
title('N = 10, t = 100');
hold on;

% bottom right, t = 200, GLOBAL METHOD
subplot(3,2,6)
plot(x,analytic(:,3),'--',...
    'displayname','Analytic');
hold on;
plot(x2,numerical_2(:,3),'displayname','Fourier');
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('$u(x,t)$','interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);
legend('interpreter', 'latex','fontsize',14);
title('N = 10, t = 200');
hold on;

exportgraphics(fig,'ex_3b.pdf','Resolution',300);
