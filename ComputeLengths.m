function lengths = ComputeLengths(nPend,nOs)

    g=9.8;

    % First pendulum (shortest) Nos oscillations per minute, second Nos+1, third Nos+2, .....
    freqs = (nOs:1:(nOs+nPend-1))/nOs;
  
    % Length freq= 1/(2*pi*sqrt(g/L)) => L = g/(freq*2*pi)^2
    lengths = g./(2*pi*freqs).^2;

end