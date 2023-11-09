rvec1 = [-6988.98,6531.78,3330.75];
rvec2 = [-6998.36,6525.51,3330.49];
vvec11 =[-5.76386,-3.84599,-0.16335];
vvec12 =[-4.03629,-2.69071,-0.113587];
vvec21 =[-5.75949,-3.85007,-0.16543];
vvec22 = [-4.03005,-2.69654,-0.116559];
mu = 3.986e5;
[a,e,I,RAAN,AOP,f] = RV2OE(rvec1,vvec11,mu);
[rrRr,vvVv,rvecORBas] = OE2RV(13000,0.3,20*pi/180,30*pi/180,50*pi/180,0.976841623137754,mu);
[rrRr2,vvVv2,rvecORBas2] = OE2RV(7226.58,0.444819,20*pi/180,30*pi/180,301.901*pi/180,2.864615316541500,mu);
I = I*180/pi;
RAAN = RAAN*180/pi;
AOP = AOP*180/pi;
f = f*180/pi;
fprintf(' %0.3f \n %0.6f \n %0.3f \n %0.3f \n %0.3f \n %0.3f \n \n',a,e,I,RAAN,AOP,f)
% p1 = 2*pi*sqrt(13000^3/mu);
% p2 = 2*pi*sqrt(a2^3/mu);
% tspan = [0 p1];
% tspan2 = [0 p2];
% rvcombo = [rrRr;vvVv];
% rvcombo2 = [rrRr2;vvVv2];
% options = odeset("RelTol",1E-12,"AbsTol",1E-12);
% [T,rvorbit] = ode45(@(t,x) TwoBP(t,x,mu),tspan,rvcombo,options);
% [T2,rvorbit2] = ode45(@(t,x) TwoBP(t,x,mu),tspan2,rvcombo2,options);
% figure(1)
% plot(rvorbit(:,1),rvorbit(:,2))
% hold on
% plot(rvorbit(1,1),rvorbit(1,2),'b.','MarkerSize',17)
% plot(rvorbit2(1,1),rvorbit2(1,2),'r.','MarkerSize',17)
% plot(rvorbit2(:,1),rvorbit2(:,2))
% grid on
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% function [dx] = TwoBP(~,x,mu)
%     dx = zeros(size(x,1),1);
%     dx(1:3) = x(4:6);
%     dx(4:6) = -mu/norm(x(1:3))^3*x(1:3);
% end