clear; clc; close all;

% INITIALISATION
x = linspace(-1,1,1000);
u = @(x,t) cos(5*pi.*x.*exp(-t)); %objective function
% dudx = @(x,t) -5*pi.*exp(-t).*sin(5*pi.*exp(-t).*x); %analytic derivative
% N = 64;
Nvec = [16 32 64 128];

% N = 128;
t0 = 0;
t = 5;
Linf_norm = zeros(1,length(Nvec));
i = 1;
for N = Nvec
    [D,xcheb] = chebyshev(N);
    % dt = t/N^2; %N^2*1e-7;
    dt = N^2*1e-7;

    Fu = @(func) -xcheb.*(D*func);
    u_n = u(xcheb,t0); % initial solution at t = 0
    for iter = dt:dt:t

        u1 = u_n + dt/2 * Fu(u_n);
        u2 = u_n + dt/2 * Fu(u1);
        u3 = u_n + dt * Fu(u2);
        un_new = (1/3) * (-u_n + u1 + 2*u2 + u3 + (dt/2) * Fu(u3));

        u_n = un_new;
    end
    tend = iter;

    err = abs(u_n - u(xcheb,tend));
%     Linf_norm = norm(err,inf);    
    Linf_norm(i) = norm(err,inf);
    i = i + 1;
end

%%
figure;
loglog(Nvec,Linf_norm(1,:),'linewidth',1.1);
title('Convergence rate at $t = 5$',...
    'interpreter','latex','fontsize',15);
xticks(Nvec);
xlabel('$N$','interpreter','latex','fontsize',16);
ylabel('$L_\infty$-error','interpreter','latex','fontsize',15);
grid on;


%%
figure;
plot(xcheb,u(xcheb,tend),'-','displayname','true cheb grid');
hold on;
plot(xcheb,u_n,'mo','displayname','chebyshev points');
hold on;
plot(x,u(x,tend),'k--','displayname','true fine grid');
title('');
legend('fontsize',14);























