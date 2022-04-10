function [A] = CentralDiff_2nd(N)
% set up matrix A for 2nd order central difference (1st deriv)
% this way we ensure symmetry of the matrix, which can be handy
% we include u_0, u_1, up to u_N-1 points, which is = N points
% I have u_0, then grid points u_1, up to u_N, where u_N = u_0
% and therefore, we only need u_0, so we have N+1 and not N+2 points

A = zeros(N+1);
stencil = [-1 0 1]; %three point stencil

A(1,1:2) = stencil(2:3);
for j = 2:N %only inner grid points
    A(j, [j-1 j j+1]) = stencil;
end
A(end,end-1:end) = stencil(1:2);

% include "boundaries", i.e., periodicity
A(1,end) = stencil(1);
A(end,1) = stencil(end);

dx = (2*pi - 0)/(N+1);
A = (1/(2*dx))*A;
end
