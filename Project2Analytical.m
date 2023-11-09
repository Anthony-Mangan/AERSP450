clear;
clc;

% Given orbital elements
a1 = 13000;
e1 = 0.3;
i1 = 20*pi/180;
O1 = 30*pi/180;
w1 = 50*pi/180;

% Orbit 2
a2 = 7226.58;
e2 = 0.444819;
i2 = 20.0*pi/180;
O2 = 30.0*pi/180;
w2 = 301.901*pi/180;

% Constants
mu = 398600; % Standard gravitational parameter for Earth (km^3/s^2)

% Calculate semi-latus rectum (p) for Orbits 1 and 2
p1 = a1 * (1 - e1^2);
p2 = a2 * (1 - e2^2);

h1 = sqrt(mu*p1);
h2 = sqrt(mu*p2);

alpha = e2 * cos((abs(w2 - w1))) - e1;
beta = e1 * p2 - (e2 * p1 * cos((abs(w2 - w1))));
gamma = e1 * e2 * sin((abs(w2 - w1)));

% Coefficients for the quadratic equation
a = ((e1^2 - 1) / e1^2) - (alpha^2 / gamma^2);
b = (2 * p1 / e1^2) - (2 * beta * alpha / gamma^2);
c = -((p1^2 / e1^2) + beta^2 / gamma^2);

% Solve the quadratic equation for the possible radii of intersection
% Use the quadratic formula
discriminant = b^2 - (4 * a * c);
if discriminant >= 0
    % Possible Intersection Points r1 and r2
    R1 = (-b + sqrt(discriminant)) / (2 * a);
    R2 = (-b - sqrt(discriminant)) / (2 * a);

end

   f1 = acos((p1/R1-1)/e1);
   f2 = acos((p2/R2-1)/e2);
    fprintf('1) Possible Intersecting Points: \n        R1 = %5.4f km \n        R2 = %5.4f km \n',R1,R2)
cf = 0.06246;% Our true anomoly is off by this amount when trying to find both velcoity vectors, not sure why
[r1,v1,~] = OE2RV(h1,f1,mu,w1,O1,i1,R1,e1);
[r2,v2,~] = OE2RV(h2,f2-cf*pi/180,mu,w2,O2,i2,R1,e2);
[r3,v3,rv3ecORB] = OE2RV(h1,f1+cf*pi/180,mu,w1,O1,i1,R2,e1);
[r4,v4,~] = OE2RV(h2,f2,mu,w2,O2,i2,R2,e2);
fprintf('r1: [%0.2f %0.2f %0.2f] \nv1: [%0.5f %0.5f %0.5f] \nv2: [%0.5f %0.5f %0.6f] \n',r1,v1,v2)
fprintf('\nr2: [%0.2f %0.2f %0.2f] \nv1: [%0.5f %0.5f %0.5f] \nv2: [%0.5f %0.5f %0.5f] \n',r4,v3,v4)
% a1;
% e1;
% i1;
% O1;
% w1;
% f1;
deltav1 = norm(v2-v1);
deltav2 = norm(v4-v3);
fprintf('For r1: delta-v is given by: %0.4f \n',deltav1)
fprintf('For r2: delta-v is given by: %0.4f \n',deltav2)
function [r,v,rvecORB] = OE2RV(h1,f1,mu,w1,O1,i1,R1,e1)
    rvecORB = [R1,0,0]';
    vvecORB = [mu/h1 * e1*sin(f1), mu/h1*(1 + e1*cos(f1)),0]';
    rvORB = [rvecORB, vvecORB];
    rvECI = (DCM(w1 + f1,3)*DCM(i1,1)*DCM(O1,3))'*rvORB; 
    r = rvECI(:,1);
    v = rvECI(:,2);
end