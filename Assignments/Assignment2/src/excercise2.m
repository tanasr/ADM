% Assignment 2 - Excercise 2

clear; 
% clc; 
close all;


%% Function a): f(x) = cos(10x)

N = 2^7;
N = linspace(2,N,N/2); %in steps of 2
pw_err = zeros(length(N),3);

i = 0; %for indexing
k = 0; %for indexing

for N = N
    
    i = i+1;
    k = k+1;
    
    j = linspace(0,N-1,N); %from 0 to N-1 with N grid points: even
    % (2*pi - 0) is right bound - left bound
    dx = (2*pi - 0)/(N); %we need to divide by the same amount of grid points
    x = j.*dx;
    
    % compute derivative
    f = cos(10*x);
    df = -10*sin(10*x); %analytic
    df_n = D_even(N)*f'; %numeric
    
    % compute the error
    err = df_n - df';
    pw_err(i,1) = N;
    pw_err(i,2) = norm(err,inf); %keep track of the error
    pw_err(k,3) = norm(err,2);
end


%% plot the errors of function a)
figure('position', [400, 400, 900, 400]);
fig = tiledlayout(1,2,'TileSpacing','loose','Padding','Compact');

nexttile;
semilogy(pw_err(:,1),pw_err(:,2)); %plot inf norm
title('$L_\infty-\mathrm{norm}$',...
    'interpreter','latex','FontSize',14);
xlabel('$N$','interpreter','latex','FontSize',14);
ylabel('Error');


nexttile;
semilogy(pw_err(:,1),pw_err(:,3)); %plot L2 norm
title('$L_2-\mathrm{norm}$',...
    'interpreter','latex','FontSize',14);
xlabel('$N$','interpreter','latex','FontSize',14);
ylabel('Error');

% exportgraphics(fig,'error_a.pdf','Resolution',300) 

%% visualise the approximation of function a)
figure('position', [400, 400, 950, 350]);
fig2 = tiledlayout(1,2,'TileSpacing','loose','Padding','Compact');


nexttile;
plot(x,df,'-o','DisplayName','$f^\prime$');
hold on;
plot(x,df_n,'-*','DisplayName','$f^\prime_n$');
xlabel('$N$','interpreter','latex','FontSize',14);
title('Derivative $f^\prime(x)$','interpreter','latex','FontSize',14);
legend('fontsize',14,'interpreter','latex');

nexttile;
plot(x,err,'DisplayName','$|f^\prime_n - f^\prime|$');
xlabel('$N$','interpreter','latex','FontSize',14);
title('Error','FontSize',14);
legend('fontsize',14,'interpreter','latex');

% exportgraphics(fig2,'error_a.pdf','Resolution',300) 




%% Function b): f(x) = cos(x/2)
%here I only compute the derivative for one single N value

N = 2^7;
j = linspace(0,N-1,N); %from 0 to N-1 with N grid points: even
dx = (2*pi - 0)/(N); %we need to divide by the same amount of grid points
x = j.*dx;

% compute derivative
f = cos(x/2);
df = -0.5*sin(x/2); %analytic
df_n = D_even(N)*f'; %numeric
err = abs(df_n - df');
inf_norm = norm(err, inf);
L2_norm = norm(err,2);


%% visualise
figure('position', [400, 400, 950, 350]);
fig3 = tiledlayout(1,2,'TileSpacing','loose','Padding','Compact');


nexttile;
plot(x,df,'-o','DisplayName','$f^\prime$');
% hold on;
% plot(x,df_n,'+','DisplayName','$f^\prime_n$');
xlabel('$x$','interpreter','latex','FontSize',14);
title('Analytic derivative $f^\prime(x)$','FontSize',14,...
    'interpreter','latex');
xlim([0 2*pi]);
% legend('fontsize',14,'interpreter','latex');
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})

nexttile;
plot(x,df_n,'*-','DisplayName','$f^\prime_n$');
xlabel('$x$','interpreter','latex','FontSize',14);
title('Approx. derivative $f^\prime(x)_n$',...
    'interpreter','latex','FontSize',14);
xlim([0 2*pi]);
% legend('fontsize',14,'interpreter','latex');
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})

% exportgraphics(fig3,'error_b.pdf','Resolution',300) 


% figure;
% plot(x,err,'DisplayName','$|f^\prime_n - f^\prime|$');
% title('Error');
% legend('fontsize',14,'interpreter','latex');




%% Function c)

% init
N = 2^7;
j = linspace(0,N-1,N); %from 0 to N-1 with N grid points: even
dx = (2*pi - 0)/(N); %we need to divide by the same amount of grid points
x = j.*dx;


% compute derivative
f = x;
% df = 1; %analytic
df = ones(1,length(df_n)); %analytic
df_n = D_even(N)*f'; %numeric
err = abs(df_n - df');


%% visualise
figure;
plot(x,df,'-o','DisplayName','$f^\prime$');
hold on;
plot(x,df_n,'-*','DisplayName','$f^\prime_n$');
title('Derivative $f^\prime(x)$','interpreter','latex');
legend('fontsize',14,'interpreter','latex');

figure;
plot(x,err,'DisplayName','$|f^\prime_n - f^\prime|$');
title('Error');
legend('fontsize',14,'interpreter','latex');