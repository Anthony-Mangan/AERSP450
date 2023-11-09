function [a,e,I,RAAN,AOP,f] = RV2OE(rvec,vvec,mu)
    h = norm(cross(rvec,vvec));
    hhat = cross(rvec,vvec)/h;
    xhat = [1,0,0];
    yhat = [0,1,0];
    zhat = [0,0,1];
    nhat = cross(zhat,hhat)/norm(cross(zhat,hhat));

    r = norm(rvec);
    rdot = dot(rvec,vvec)/r;
    v = norm(vvec);
    evec = cross(vvec,hhat*h)/mu - rvec/r;
    a = 1/(2/r - v^2/mu);
    e = norm(evec);
    I = acos(dot(hhat,zhat));
    cOmega = dot(nhat,xhat);
    sOmega = dot(nhat,yhat);
    if cOmega >= 0
        RAAN = atan(sOmega/cOmega);
    else
        RAAN = atan(sOmega/cOmega) + pi;
    end
    if dot(evec,zhat) < 0
        AOP = 2*pi-acos(dot(evec,nhat)/e);
    else
        AOP = acos(dot(evec,nhat)/e);
    end
    if rdot > 0
        f = acos(1/e*(h^2/(mu*r) -1));
    elseif rdot < 0
        f = 2*pi-acos(1/e*(h^2/(mu*r) -1));
    else
        if round(r,8) == round(a*(1-e),8)
            f = 0;
        else 
            f = pi/2;
        end
    end





end