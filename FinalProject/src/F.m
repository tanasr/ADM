function f = F(func,K)
    freq = -K:K;
    N = 2*K;
    f = zeros(1,length(freq));
    for n = freq
        for k = freq
            if abs(n-k) <= N/2
            f(n+K+1) = f(n+K+1) - (1i*k).*func(n-k+K+1).*func(k+K+1);   
            end
        end
        f(n+K+1) = f(n+k+1) - 0.1*n^2.*func(n+K+1);
    end
end














