function dxdt = deriv(t,state,data)

    n = size(state,1)/2;
    
    x   = state(1:n);
    dx  = state(n+1:end);
    ddx = -data.g*(sin(x)./data.l);
    
    dxdt = zeros (2*n,1);
    dxdt(1:n,1) = dx;
    dxdt(n+1:end,1) = ddx;  

end