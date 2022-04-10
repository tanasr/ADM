function u_n = RungeKutta_4(un,dudx,dt)
    
    
    Fu = @(method,un) -2*pi*method(un);
    
    u1 = un + dt/2 * Fu(dudx,un);
    u2 = un + dt/2 * Fu(dudx,u1);
    u3 = un + dt * Fu(dudx,u2);
    un_new = (1/3) * (-un + u1 + 2*u2 + u3 + (dt/2) * Fu(dudx,u3));

    % for each time step, there is a new un
    u_n = un_new; %update

end