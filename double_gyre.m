function dx = double_gyre(tau,x)
idx1 = x(1);
idx2 = x(2);
dx = nan(size(x));

dx(1) = -pi*0.1*sin(pi*(0.25*sin(2*pi*tau/10)*idx1.^2+idx1-2*0.25*sin(2*pi*tau/10)*idx1)).*cos(pi*idx2);
dx(2) = pi*0.1*cos(pi*(0.25*sin(2*pi*tau/10)*idx1.^2+idx1-2*0.25*sin(2*pi*tau/10)*idx1)).*sin(pi*idx2).*(2*0.25*sin(2*pi*tau/10)*idx1+1-2*0.25*sin(2*pi*tau/10));