clear; clc; close all;

% INITIALISATION
% N = 64; %even number
Nvec = [16 32 48 64 96 128 192 256];
c = 4.0;
v = 0.1;

Linf_norm = zeros(1,length(Nvec));
i = 1;
for N = Nvec
    j = linspace(0,N,N+1); %odd number of grid points
    dx = (2*pi - 0)/(N+1); % --> same as Ts in the example below
    x = j.*dx; % --> same as t = 0:Ts:T-Ts in the example below

    t0 = 0;
    t = pi/4;
    CFL = 1.4;

    % 1.) Define u(x,0)
    u = @(x,t) c-2*v.*(d_phi_dx(x,c,v,t)./phi_func(x-c*t, t+1, v)); 
    u_0 = u(x,0);

    % 2.) Calculate complex-valued coefficients
    % N = N; % number of spectral components --> same as % K = 10
    K = (N/2); % N = 2*K+1
    a_k0 = complex(ones(1,N+1)); %from -N to +N including 0
    modes = -K:K;
    for k = modes
        integrand = @(x) u(x,0).*exp(-1i*k*x);
        a_k0(k+K+1) = (1/(2*pi)) * integral(integrand,0,2*pi);
    end
%     figure;
%     stem(modes,real(a_k0));
%     hold on;
%     stem(modes,imag(a_k0));
%     hold off;
%     legend('real','imag')


    % this is not necessary, just a test here if it's correct. we actually need to
    % compute the derivative using FourierGalerkin and for that, we need the coefficients
    % 3.) Fourier series complex exponentials form
    u_n0 = zeros(1,length(x));
    for k = modes
        u_n0 = u_n0 + (a_k0(k+K+1) .* exp(1i*k*x)); %approx. of the initial solution
    end
    % figure;
    % plot(x,u_0,'displayname','analytic');
    % hold on;
    % plot(x,real(u_n0),'displayname','fourier');
    % title('fourier vs. analytic for initial condition');
    % legend('fontsize',14);

    % err = real((u_n0 - u_0));
    % Linf_norm = norm(err,inf); %error measured on initial solution



    % ADVANCE IN TIME USING RUNGE KUTTA
    dt = CFL .* (max(abs(u(x,t0)).*(N/2) + v.*(N/2)^2))^-1;
    for iter = dt:dt:t

        u1 = a_k0 + dt/2 * F(a_k0,K);
        u2 = a_k0 + dt/2 * F(u1,K);
        u3 = a_k0 + dt * F(u2,K);
        a_knew = (1/3) * (-a_k0 + u1 + 2*u2 + u3 + (dt/2) * F(u3,K));

        a_k0 = a_knew;
    end
    tend = iter;

    u_n = zeros(1,length(x));
    for k = modes
        u_n = u_n + (a_knew(k+K+1) .* exp(1i*k*x));
    end

    err = abs(u_n - u(x,tend));
    Linf_norm(i) = norm(err,inf);
%     Linf_norm = norm(err,inf);
    i = i + 1;
end

figure;
plot(x,real(u_n),'displayname','fourier');
hold on;
plot(x,u(x,tend),'displayname','analytic');
legend('fontsize',14);


%% Visualise convergence rate

figure;
loglog(Nvec,Linf_norm(1,:),'linewidth',1.2);
title('Convergence rate for $t = \pi/4$ (Fourier-Galerkin method)',...
    'interpreter','latex','fontsize',15);
xticks(Nvec)
xlabel('$N$','interpreter','latex','fontsize',16);
ylabel('$L_\infty$-error','interpreter','latex','fontsize',15)
grid on;


%% Visualise the final plot at t = pi/4 for N = 128
figure;
plot(x,u_n,'o','displayname','numerical');
hold on;
plot(x, u(x,tend),'displayname','analytical');
legend('fontsize',14);
title(['N = ',+num2str(N), ', dt = ',+num2str(dt),', $L_\infty$ = ',+num2str(Linf_norm(1,end)),', CFL = ',+num2str(CFL)],...
    'interpreter','latex','fontsize',14);
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel(['$u(x,t = ',+num2str(tend),')$'],'interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);





%% Plot for N = 128 and t = [0 pi/8 pi/6 pi/4]
% clear; clc; close all;

% INITIALISATION
N = 64;
c = 4.0;
v = 0.1;

fig = figure;
i = 1;
for t = [pi/8 pi/6 pi/4]
    j = linspace(0,N,N+1); 
    dx = (2*pi - 0)/(N+1); 
    x = j.*dx;

    t0 = 0;
    % t = pi/4;
    CFL = 1.4;

    % 1.) Define u(x,0)
    u = @(x,t) c-2*v.*(d_phi_dx(x,c,v,t)./phi_func(x-c*t, t+1, v)); 
    u_0 = u(x,0);

    % 2.) Calculate complex-valued coefficients
    K = (N/2);
    a_k0 = complex(ones(1,N+1)); 
    modes = -K:K;
    for k = modes
        integrand = @(x) u(x,0).*exp(-1i*k*x);
        a_k0(k+K+1) = (1/(2*pi)) * integral(integrand,0,2*pi);
    end

    u_n0 = zeros(1,length(x));
    for k = modes
        u_n0 = u_n0 + (a_k0(k+K+1) .* exp(1i*k*x)); 
    end

    % ADVANCE IN TIME USING RUNGE KUTTA
    dt = CFL .* (max(abs(u(x,t0)).*(N/2) + v.*(N/2)^2))^-1;
    for iter = dt:dt:t

        u1 = a_k0 + dt/2 * F(a_k0,K);
        u2 = a_k0 + dt/2 * F(u1,K);
        u3 = a_k0 + dt * F(u2,K);
        a_knew = (1/3) * (-a_k0 + u1 + 2*u2 + u3 + (dt/2) * F(u3,K));

        a_k0 = a_knew;
    end
    tend = iter;

    % get back to space dimension
    u_n = zeros(1,length(x));
    for k = modes
        u_n = u_n + (a_knew(k+K+1) .* exp(1i*k*x));
    end

    err = abs(u_n - u(x,tend));
    Linf_norm = norm(err,inf);


    subplot(2,2,i+1)
    plot(x,u_n,'o','displayname','numerical');
    hold on;
    plot(x, u(x,tend),'displayname','analytical');
    legend('fontsize',14);
    title(['Fourier-Galerkin method for $t = $',+num2str(t)],...
        'interpreter','latex','fontsize',14);
    xlim([0 2*pi]);
    xticks(0:pi/2:2*pi)
    xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
    ylabel(['$u(x,t = ',+num2str(tend),')$'],'interpreter','latex','fontsize',16);
    xlabel('$x$','interpreter','latex','fontsize',16);

    i = i + 1;
end

%% include plot for t = 0
hold on;
subplot(2,2,1)
plot(x,u(x,0),'o','displayname','numerical');
hold on;
plot(x, u(x,0),'displayname','analytical');
legend('fontsize',14);
title(['$t = 0$'],...
    'interpreter','latex','fontsize',14);
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel(['$u(x,t = ',+num2str(tend),')$'],'interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);

exportgraphics(fig,'final_FG.pdf','Resolution',300)



%% Visualise the results
% figure;
% plot(x,u_n,'displayname','numerical');
% hold on;
% plot(x, u_0(x,tend),'displayname','analytical');
% legend('fontsize',14);
% title(['N = ',+num2str(N), ', dt = ',+num2str(dt),', $L_\infty$ = ',+num2str(Linf_norm),', CFL = ',+num2str(CFL)],...
%     'interpreter','latex','fontsize',14);
% xlim([0 2*pi]);
% xticks(0:pi/2:2*pi)
% xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
% ylabel(['$u(x,t = ',+num2str(tend),')$'],'interpreter','latex','fontsize',16);
% xlabel('$x$','interpreter','latex','fontsize',16);






