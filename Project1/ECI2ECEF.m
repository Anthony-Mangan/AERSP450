function[r_ecef,v_ecef] = ECI2ECEF(r_eci,v_eci,omega, t, GMST)
    rvECI = [r_eci,v_eci];
    rvECEF = DCM(GMST + omega*t,3)*rvECI; 
    r_ecef = rvECEF(:,1);
    v_ecef = rvECEF(:,2);
end
% mu is not used in the function?