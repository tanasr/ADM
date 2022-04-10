% clear; clc; close all;

% INITIALISATION
% N = 256; %even number
c = 4.0;
v = 0.1;
Nvec = [16 32 48 64 96 128 192 256];

t0 = 0;
t = pi;
CFL = 0.75; % up until N = 96, CFL = 1 is fine
% for N = 128, CFL = 0.99 is stable but not good, CFL = 0.98 is good
% for N = 192, CFL = 0.86
% for N = 256, CFL = 0.75

Linf_norm = zeros(1,length(Nvec));
i = 1;
for N = Nvec
    j = linspace(0,N,N+1); %odd number of grid points
    dx = (2*pi - 0)/(N+1);
    x = j.*dx;

    %ufunc with analytic derivative of function phi
    u_func = @(x,t) c-2*v.*(d_phi_dx(x,c,v,t)./phi_func(x-c*t, t+1, v)); 
    %ufunc with numerical derivative of function phi
    % u_func = @(x,t) c-2*v.*(dphidx_n./phi_func(x-c*t, t+1));

    % function to compute derivative
    dudx = @(func) transpose(D_odd(N) * func'); %transpose output for RungeKutta algo
    dt = CFL * (max(abs(u_func(x,t0))/(dx) + v/(dx)^2))^-1;
    Fu = @(un) -un.*dudx(un) + v.*dudx(dudx(un));
    u_n = u_func(x,t0); % initial solution at t = 0
    for iter = dt:dt:t

        u1 = u_n + dt/2 * Fu(u_n);
        u2 = u_n + dt/2 * Fu(u1);
        u3 = u_n + dt * Fu(u2);
        un_new = (1/3) * (-u_n + u1 + 2*u2 + u3 + (dt/2) * Fu(u3));

        u_n = un_new;
    end
    tend = iter;


    err = abs(u_n - u_func(x,tend));
    Linf_norm(i) = norm(err,inf);
    i = i + 1;
end
% CFL = CFL - 0.005;
% end



%% Visualise convergence rate

figure;
loglog(Nvec,Linf_norm(1,:),'linewidth',1.2);
title('Convergence rate for $t = \pi/4$',...
    'interpreter','latex','fontsize',15);
xticks(Nvec)
xlabel('$N$','interpreter','latex','fontsize',16);
ylabel('$L_\infty$-error','interpreter','latex','fontsize',15)
grid on;



%% Visualise the final plot at t = pi/4 for N = 128
figure;
plot(x,u_n,'o','displayname','numerical');
hold on;
plot(x, u_func(x,tend),'displayname','analytical');
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
N = 128; %even number
c = 4.0;
v = 0.1;

t0 = 0;
CFL = 0.75; % up until N = 96, CFL = 1 is fine

fig = figure;
i = 1;
for t = [pi/8 pi/6 pi/4]
    j = linspace(0,N,N+1);
    dx = (2*pi - 0)/(N+1);
    x = j.*dx;

    u_func = @(x,t) c-2*v.*(d_phi_dx(x,c,v,t)./phi_func(x-c*t, t+1, v)); 
    dudx = @(func) transpose(D_odd(N) * func'); 
    dt = CFL * (max(abs(u_func(x,t0))/(dx) + v/(dx)^2))^-1;
    Fu = @(un) -un.*dudx(un) + v.*dudx(dudx(un));
    u_n = u_func(x,t0);
    for iter = dt:dt:t

        u1 = u_n + dt/2 * Fu(u_n);
        u2 = u_n + dt/2 * Fu(u1);
        u3 = u_n + dt * Fu(u2);
        un_new = (1/3) * (-u_n + u1 + 2*u2 + u3 + (dt/2) * Fu(u3));

        u_n = un_new;
    end
    tend = iter;

    err = abs(u_n - u_func(x,tend));
    Linf_norm(i) = norm(err,inf);

    subplot(2,2,i+1)
    plot(x,u_n,'o','displayname','numerical');
    hold on;
    plot(x, u_func(x,tend),'displayname','analytical');
    legend('fontsize',14);
    title(['$t = $',+num2str(t)],...
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
plot(x,u_func(x,0),'o','displayname','numerical');
hold on;
plot(x, u_func(x,0),'displayname','analytical');
legend('fontsize',14);
title(['$t = 0$'],...
    'interpreter','latex','fontsize',14);
xlim([0 2*pi]);
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel(['$u(x,t = ',+num2str(tend),')$'],'interpreter','latex','fontsize',16);
xlabel('$x$','interpreter','latex','fontsize',16);

% exportgraphics(fig,'final_sol.pdf','Resolution',300)














