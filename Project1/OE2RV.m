function [r,v,rvecORB] = OE2RV(a,e,I,RAAN,AOP,f,mu)
    p = a*(1-e^2);
    h = sqrt(mu*p);
    rnorm = p/(1 + e*cos(f));
    rvecORB = [rnorm,0,0]';
    vvecORB = [mu/h * e*sin(f), mu/h*(1 + e*cos(f)),0]';
    rvORB = [rvecORB, vvecORB];
    rvECI = (DCM(AOP + f,3)*DCM(I,1)*DCM(RAAN,3))'*rvORB; 
    r = rvECI(:,1);
    v = rvECI(:,2);
end

