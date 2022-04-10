function [D] = D_even(n)
D = zeros(n,n);

%we can use the vector j and perform iterations over columns only
j = linspace(0,n-1,n); %from 0 to N-1 with N grid points: even
    for i = 0:n-1
        D(:,i+1) = ((-1).^(i+j)./2).*cot(((j-i).*pi)/(n));
%         D(:,i+1) = ((-1).^(i+j)./2).*cot((x(j+1)-x(i+1))/2);
        % to use this, add vector x as function input

        %set the diagonal entries equal 0
        D(i+1,i+1) = 0;
    end
end

% 
% function [D] = D(n)
% D = zeros(n+1,n+1);
% for i = 1:n+1
%     for j = 1:n+1
%         if i ~= j
%             D(j,i) = ((-1)^(j+i)/2).*(sin(((j-i)*pi)/(n+1)))^(-1);
%         end
%     end
% end
% end