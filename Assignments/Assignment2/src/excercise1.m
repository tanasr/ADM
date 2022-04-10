% Assignment 2 - Excercise 1
% Implementation of the Fourier Differentiation method (even) and comparison with odd
clear; 
close all;


%% Fourier Differentiation Matrix (odd method)
k_values = [2 4 6 8 10 12]; 

for k = k_values
    for N = 10:2:80 %N is an even integer

        j = linspace(0,N-1,N); %from 0 to N-1 with N grid points: even
        % (2*pi - 0) is right bound - left bound
        dx = (2*pi - 0)/(N); %we need to divide by the same amount of grid points
        x = j.*dx;

        %objective function
        u = exp(k*sin(x));

        %compute its derivative
        analytic = k*(exp(k*sin(x))).*cos(x);
        approx = D_even(N)*u';
%         approx = D_even(x,N)*u';

        %compute the deviation form the analytical solution
        err = approx' - analytic;
        rel_err = err./analytic;
        err_norm = norm(err,inf);

        if norm(err,inf) <= 1e-5
            disp(['k = ', num2str(k)])
            disp(['n = ', num2str(N)])
            disp(err_norm)
            break
        end
    end
end


%% Visualise

figure('position',[100 100 800 300])
xlabel('$x$','interpreter','latex','fontsize',14);

yyaxis left;
plot(x,u,'k:','Linewidth',1.2,'DisplayName',"$u(x)$");
hold on;
plot(x, analytic,'-','LineWidth',1.7,'DisplayName',"$u'(x)$");
hold on;
plot(x, approx,'--k','LineWidth',1.4,'DisplayName',"$u'_n(x) = D\cdot u$");
hold on;
ylabel('Function value');
hold on;


%add the error to the plot
yyaxis right;
ylabel('Error');
plot(x,abs(err),'g','DisplayName','$|u^\prime_n - u^\prime|$');
hold on;
plot(x,rel_err,'-o','DisplayName','$\frac{u^\prime_n - u^\prime}{u^\prime}$');
hold on;


title(['Analytic vs. approx. derivative with N = ',num2str(N)],...
    ['obj. func.: $u(x) = \exp(\mathrm{k}\sin(x))$ for $k = $',num2str(k)], ...
    'Interpreter', 'Latex','Fontsize',16); 
leg = legend;
leg.set('Interpreter', 'Latex', 'Fontsize',12); 
hold off;


