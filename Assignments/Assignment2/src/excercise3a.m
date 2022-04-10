% Assignment 2 - Exercise 3

clear; clc; close all;

% set parameters for all methods
N = [8 16 32 64 128 256 512 1024 2048];
Linf_norm = zeros(length(N),4);
Linf_norm(:,1) = N'; %set the N values for the x-axis

t0 = 0;
t = pi;
% dt = 1e-2;
dt = 1e-3/3;
% dt = 1e-4;
%% Approximation using 2nd-order Central Difference
N = [8 16 32 64 128 256 512 1024 2048];

n = 0; % for indexing of the Linf_vector
for N = N
    n = n+1;

    j = linspace(0,N,N+1); %odd number of grid points
    dx = (2*pi - 0)/(N+1);
    x = j.*dx;
    u = @(x,t) exp(sin(x - 2*pi*t)); %analytic solution 

    % 2nd order central difference operator for spatial derivative
    A = CentralDiff_2nd(N);
    du_dx_CD2 = @(u) (A * u')'; %transpose is necessary for RungeKutta


    % change parameters here if you want other values than defined above
%     t0 = t0;
%     t = t;
    % resolution = 5; %factor to make dt smaller, i.e., to make grid finer
    % dt = 1/resolution * (t - t0)/(N+1);
    % dt = 1e-3;

    u0 = u(x,t0); % initial solution at t = 0
    u_n = u0;
    for iter = dt:dt:t
        u_n = RungeKutta_4(u_n,du_dx_CD2,dt); %advance in time
    end
    tend = iter;
    % tend = t;


    % analytic solution for comparison
    u_analytic = u(x,tend);
    err = abs(u_n - u_analytic);
    % Linf_CD2 = norm(err,inf);

    Linf_norm(n,2) = norm(err,inf);
end



%% Approximation using 4th-order Central Difference
N = [8 16 32 64 128 256 512 1024 2048];
% Linf_norm = zeros(length(N),2);

n = 0;
for N = N
    n = n+1;

    j = linspace(0,N,N+1); %odd number of grid points
    dx = (2*pi - 0)/(N+1);
    x = j.*dx;
    u = @(x,t) exp(sin(x - 2*pi*t)); %solution 

    % 4th order central difference operator for spatial derivative
    A = CentralDiff_4th(N);
    du_dx_CD4 = @(u) (A * u')'; %transpose is necessary for RK

    % change parameters here if you want other values than defined above
%     t0 = t0;
%     t = t;
    % resolution = 100; %factor to make dt smaller, i.e., to make grid finer
    % dt = 1/resolution * (tend - t0)/(N+1);
    % dt = 1e-3;

    u0 = u(x,t0); % initial solution at t = 0
    u_n = u0;
    for iter = dt:dt:t
        u_n = RungeKutta_4(u_n,du_dx_CD4,dt); %advance in time
    end
    tend = iter;
    % tend = t;


    % analytic solution for comparison
    u_analytic = u(x,tend);
    err = abs(u_n - u_analytic);
    % Linf_CD4 = norm(err,inf);

    Linf_norm(n,3) = norm(err,inf);
end




%% Approximation using global Fourier method
N = [8 16 32 64 128 256 512 1024 2048];
% Linf_norm = zeros(length(N),2);

n = 0;
for N = N
    n = n+1;

    j = linspace(0,N,N+1); %odd number of grid points
    dx = (2*pi - 0)/(N+1);
    x = j.*dx;
    u = @(x,t) exp(sin(x - 2*pi*t)); %solution 

    % numerical derivative 2nd order CD
    du_dx_GM = @(u) (D_odd(N) * u')'; %transpose is necessary for RK


    % change parameters here if you want other values than defined above
%     t0 = t0;
%     t = t;
    % resolution = 10; %factor to make dt smaller, i.e., to make grid finer
    % dt = 1/resolution * (t - t0)/(N+1);
    % dt = 1e-3;

    u0 = u(x,t0); % initial solution at t = 0
    u_n = u0; 
    for iter = dt:dt:t 
    u_n = RungeKutta_4(u_n,du_dx_GM,dt); %advance in time
    end
    tend = iter; %this is the last point of the iteration


    % analytic solution for comparison
    u_analytic = u(x,tend);
    err = abs(u_n - u_analytic);
    Linf_GM = norm(err,inf);


    Linf_norm(n,4) = norm(err,inf);
end


%% visualise the infinity norm over values of N
fig = figure('position',[100 100 650 400]);
loglog(Linf_norm(:,1),Linf_norm(:,2),...
    'displayname','$2^{\mathrm{nd}}$-order CD','linewidth',1.5);
hold on;
loglog(Linf_norm(:,1),Linf_norm(:,3),...
    'displayname','$4^{\mathrm{th}}$-order CD','linewidth',1.5);
hold on;
loglog(Linf_norm(1:3,1),Linf_norm(1:3,4),...
    'displayname','Global method','linewidth',1.5);
hold on;
title(['The pointwise error $L_\infty$-norm at t = ', ...
    + num2str(tend),' for $\Delta t$ = ',+num2str(dt)], 'interpreter','latex','fontsize',14);
ylabel('$L_\infty$','fontsize',16,'interpreter','latex');
xlabel('Gridpoints N','fontsize',14,'interpreter','latex');
xlim([4 10^3.5])
legend('interpreter','latex','fontsize',14);

% savefig(fig,'Linf_norm.pdf');














