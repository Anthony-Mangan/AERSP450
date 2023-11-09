function [dx] = TwoBP(~,x,mu)
    dx = zeros(size(x,1),1);
    dx(1:3) = x(4:6);
    dx(4:6) = -mu/norm(x(1:3))^3*x(1:3);
end