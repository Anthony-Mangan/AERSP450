%% Project 1
% Anthony Mangan, Arnold Vento, Sadhana Dasari
clear 
clc
%% Question 1
r0 = [-7669.55591738356;4428.02030100983;7698.43847858985];
v0 = [-3.32261842988279;-5.75494384030483;4.1037677802895e-08];
rvcombo = [r0;v0];
dt = 10;
tspan = 0:dt:6*3600*24;
mu = 3.986E5;
options = odeset("RelTol",1E-12,"AbsTol",1E-12);
[T,rvorbit] = ode45(@(t,x) TwoBP(t,x,mu),tspan,rvcombo,options);
figure(1)
plot3(rvorbit(:,1),rvorbit(:,2),rvorbit(:,3))
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
%% Question 2
rECEF = zeros(length(T),3);
omega = 2*pi/86400;
t = 0;
GMST = 228*pi/180;
for i = 1:size(rvorbit,1)
    [r_ecef,~] = ECI2ECEF(rvorbit(i,1:3)',rvorbit(i,4:6)',omega, t, GMST);
    rECEF(i,:) = r_ecef';
end
figure(2)
subplot(2,2,1)
plot3(rECEF(:,1),rECEF(:,2),rECEF(:,3))
title("3D")
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
subplot(2,2,2)
plot(rECEF(:,1),rECEF(:,2))
title("X-Y")
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
subplot(2,2,3)
plot(rECEF(:,2),rECEF(:,3))
title("Y-Z")
grid on
xlabel('Y')
ylabel('Z')
zlabel('Z')
subplot(2,2,4)
plot(rECEF(:,1),rECEF(:,3))
title("X-Z")
grid on
xlabel('X')
ylabel('Z')
zlabel('Z')
%% Question 3
figure(3)
load('topo.mat','topo');
topoplot = [topo(:,181:360),topo(:,1:180)];
contour(-180:179,-90:89,topoplot,[0,0],'black','linewidth',1);
grid on
grid minor
axis equal
hold on
[long,lat] = LongLat(rECEF(:,1),rECEF(:,2),rECEF(:,3));
plot(long,lat,'m.','LineWidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
set(gca,'FontSize',18)
%% Question 4
r0 = [-7669.55591738356;4428.02030100983;7698.43847858985];
v0 = [-3.32261842988279;-5.75494384030483;4.1037677802895e-08];
tspan = 6*3600*24;
dt = 10;
tt = 0:dt:tspan;
mu = 3.986E5;
for i = 1:length(tt)
    [rfgf,vfgf] = FGFunc(r0,v0,dt,mu);
    rFGF(:,i) = real(rfgf);
    vFGF(:,i) = real(vfgf);
    r0 = rfgf;
    v0 = vfgf;
end
figure(4)
plot3(rFGF(1,:),rFGF(2,:),rFGF(3,:))
grid on
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
%% Question 5
rvFGF = [rFGF;vFGF]';
error = (rvFGF - rvorbit)./rvFGF;
figure(5)
subplot(2,3,1)
plot(tt,(rvorbit(:,1) - rvFGF(:,1))./norm(rvorbit(:,1) - rvFGF(:,1)))
xlabel('Time')
ylabel('X Position Error')
subplot(2,3,2)
plot(tt,(rvorbit(:,2) - rvFGF(:,2))./norm(rvorbit(:,2) - rvFGF(:,2)))
xlabel('Time')
ylabel('Y Position Error')
subplot(2,3,3)
plot(tt,(rvorbit(:,3) - rvFGF(:,3))./norm(rvorbit(:,3) - rvFGF(:,3)))
xlabel('Time')
ylabel('Z Position Error')
subplot(2,3,4)
plot(tt,(rvorbit(:,4) - rvFGF(:,4))./norm(rvorbit(:,4) - rvFGF(:,4)))
xlabel('Time')
ylabel('X Velocity Error')
subplot(2,3,5)
plot(tt,(rvorbit(:,5) - rvFGF(:,5))./norm(rvorbit(:,5) - rvFGF(:,5)))
xlabel('Time')
ylabel('Y Velocity Error')
subplot(2,3,6)
plot(tt,(rvorbit(:,6) - rvFGF(:,6))./norm(rvorbit(:,6) - rvFGF(:,6)))
xlabel('Time')
ylabel('Z Velocity Error')
sgtitle('Numerical and Analytic Propagation Error')
%% Question 6
[a,e,I,RAAN,AOP,f] = RV2OE(rvorbit(2,1:3)',rvorbit(2,4:6)',mu)
%% Functions
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
        AOP = -acos(dot(evec,nhat)/e);
    else
        AOP = acos(dot(evec,nhat)/e);
    end
    if rdot > 0
        f = acos(1/e*(h^2/(mu*r) -1));
    elseif rdot < 0
        f = -acos(1/e*(h^2/(mu*r) -1));
    else
        if round(r,8) == round(a*(1-e),8)
            f = 0;
        else 
            f = pi/2;
        end
    end
end

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

function[r_ecef,v_ecef] = ECI2ECEF(r_eci,v_eci,omega, t, GMST)
    rvECI = [r_eci,v_eci];
    rvECEF = DCM(GMST + omega*t,3)*rvECI; 
    r_ecef = rvECEF(:,1);
    v_ecef = rvECEF(:,2);
end
% mu is not used in the function?

function[output] = DCM(angle,axis)
    switch axis
        case 1
            output = [1,         0,             0;
                      0,    cos(angle),    sin(angle);
                      0,    -sin(angle),   cos(angle)];
        case 2
            output =  [cos(angle),    0,       -sin(angle);
                           0,         1,           0;
                       sin(angle),    0,        cos(angle)];
         case 3
                   
            output =  [cos(angle),     sin(angle),       0;
                       -sin(angle),    cos(angle),       0;
                            0,              0,           1];
    end
     
end
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
function [E1] = NR(e,M1)
error = 1;
absTol = 1E-8;
Eold = M1;
    while error > absTol
        y = M1 + e*sin(Eold) - Eold;
        z = e*cos(Eold) - 1;
        Enew = Eold - y/z;
        error = abs(Enew-Eold);
        Eold = Enew;
    end
E1 = Enew;
end
function [dx] = TwoBP(~,x,mu)
    dx = zeros(size(x,1),1);
    dx(1:3) = x(4:6);
    dx(4:6) = -mu/norm(x(1:3))^3*x(1:3);
end
function [long,lat] = LongLat(x,y,z)
    long = zeros(size(x));
    lat = long;
    for i = 1:length(x)
        long(i) = atan2(y(i),x(i))*180/pi;
        lat(i) = atan2(z(i),sqrt(x(i)^2 + y(i)^2))*180/pi;
    end
end




