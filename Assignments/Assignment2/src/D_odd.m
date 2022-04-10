function [D] = D_odd(n)
D = zeros(n+1,n+1);

%we can use the vector j and perform iterations over columns only
j = linspace(0,n,n+1);
    for i = 0:n
        D(:,i+1) = ((-1).^(j+i)./2).*(sin(((j-i).*pi)/(n+1))).^(-1);
        % D(:,i+1) = ((-1).^(j+i)./2).*(sin((x(j+1)-x(i+1))/(2)).^(-1));
        % to use this, add vector x as function input

        %set the diagonal entries equal 0
        D(i+1,i+1) = 0;
    end
end

% 
% function [D] = D_odd(n)
% D = zeros(n+1,n+1);
% for i = 1:n+1
%     for j = 1:n+1
%         if i ~= j
%             D(j,i) = ((-1)^(j+i)/2).*(sin(((j-i)*pi)/(n+1)))^(-1);
%         end
%     end
% end
% end