mu = 3.9860044e5;
a = (mu*(6*3600)^2/(4*pi^2))^(1/3);
e = 0.3;
I = 0*pi/180;
AOP = 0*pi/180;
RAAN = 0*pi/180;
f = 0*pi/180;
omega = 2*pi/(86400);
[r,v,ashdf] = OE2RV(a,e,I,RAAN,AOP,f,mu);
[r_ecef,v_ecef] = ECI2ECEF(r,v,omega, 0, 175.4942*pi/180)
norm(r) * norm(v)