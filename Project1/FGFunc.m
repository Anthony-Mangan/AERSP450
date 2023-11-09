function [r1,v1] = FGFunc(r0,v0,dt,mu)
    [a,e,~,~,~,f0] = RV2OE(r0,v0,mu);
    n = sqrt(mu/a^3);
    E0 = 2*atan(sqrt((1-e)/(1+e))*tan(f0/2));
    M0 = E0 - e*sin(E0);
    M1 = M0 + n*dt;
    E1 = NR(e,M1);
    F = 1 - a/norm(r0)*(1-cos(E1-E0));
    G = dt - 1/n*((E1-E0) - sin(E1-E0));
    r1 = F*r0 + G*v0;
    Fdot = -sqrt(mu*a)/(norm(r1)*norm(r0))*sin(E1-E0);
    Gdot = 1 - a/norm(r1)*(1-cos(E1-E0));
    v1 = Fdot*r0 + Gdot*v0;
end