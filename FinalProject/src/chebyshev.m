function [D,x] = chebyshev(N)
    % returns D = Differentiation matrix
    % returns x = Chebyshev Gauss–Lobatto quadrature points
    
    D = zeros(N+1,N+1);
    %j goes from 0 to N, including 0, so N+1 points
    j = linspace(0,N,N+1); 
    % Chebyshev Gauss–Lobatto quadrature points (equally distanced on the arc)
    x = -cos(j*pi/N)';
    % at endpoints = 2, at inner points = 1
    c = [2; ones(N-1,1); 2];% .* (1).^(0:N)';
    
    for i = 0:N
        for j = 0:N
            
            if i == 0 && j == 0
                D(i+1,j+1) = -(2*N.^2+1)/6;
            end
            
            if i == N && j == N
                D(i+1,j+1) = (2*N^2+1)/6;
            end
            
            if i ~= j
                D(i+1,j+1) = c(i+1)/c(j+1) .* ((-1)^(i+j+2))/(x(i+1) - x(j+1));
            end
            
            if i == j && i ~= 0 && i ~= N
                D(i+1,j+1) = -x(i+1)./(2*(1-x(i+1).^2));
            end
        end
    end
end