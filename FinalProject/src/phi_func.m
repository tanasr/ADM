function phi = phi_func(a,b,v)%phi_func(a,b,v,j)

    phi = 0;
    for k = -100:100 %it's from -inf to inf, but for small time range it's sufficient
        phi = phi + exp((-(a-(2*k+1).*pi).^2)/(4*v*b));
    end

end