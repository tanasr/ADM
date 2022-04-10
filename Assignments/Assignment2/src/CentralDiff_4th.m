function [A] = CentralDiff_4th(N)
% set up matrix A for 4th order central difference (1st deriv.)
% this way we ensure symmetry of the matrix, which can be handy
% we include u_0, u_1, up to u_N-1 points, which is = N points
% I have u_0, then grid points u_1, up to u_N, where u_N = u_0
% and therefore, we only need u_0, so we have N+1 and not N+2 points

A = zeros(N+1);
stencil = [1 -8 0 +8 -1]; %five point stencil

A(1,1:3) = stencil(3:end); %first row left
A(1,end-1:end) = stencil(1:2); %first row wrap around
A(2,1:4) = stencil(2:end); %second row left
A(2,end) = stencil(1); %second row wrap around
for j = 3:N-1 %only inner grid points
    A(j, [j-2 j-1 j j+1 j+2]) = stencil; % third row until third last
end
A(end-1,end-3:end) = stencil(1:4); % second last row
A(end-1,1) = stencil(end); %second last row wrap around
A(end,end-2:end) = stencil(1:3); % last row
A(end,1:2) = stencil(end-1:end); %last row wrap around

dx = (2*pi - 0)/(N+1);
A = (1/(12*dx))*A;
end
