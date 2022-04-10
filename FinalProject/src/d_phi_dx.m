function f = d_phi_dx(x, c, v, t)

    f = 0;
    for k = -100:100 %it's from -inf to inf
        f = f + ((exp((-(x-c*t-(2*k+1).*pi).^2)./(4*v*(t+1)))).*...
            (-1*(x-c*t-(2*k+1)*pi)))./(2*v*(t+1));
    end
    
end