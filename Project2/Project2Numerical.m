clc;
clear;

alpha1 = 48.699;
alpha2 = 48.684;
alpha3 = 48.67;

delta1 = -51.221;
delta2 = -51.054;
delta3 = -50.887;

R1 = [-1522.8,-5212.8,3386.4];
R2 = [-1521,-5213.3,3386.4];
R3 = [-1519.1,-5213.9,3386.4];

tau1 = -4;
tau3 = 5;
tau31 = tau3 - tau1;

L1 = [cos(alpha1)*cos(delta1);cos(delta1)*sin(alpha1);sin(delta1)];
L2 = [cos(alpha2)*cos(delta2);cos(delta2)*sin(alpha2);sin(delta2)];
L3 = [cos(alpha3)*cos(delta3);cos(delta3)*sin(alpha3);sin(delta3)];

D0 = dot(L1,cross(L2,L3));
D11 = dot(R1,cross(L2,L3));
D21 = dot(R2,cross(L2,L3));
D31 = dot(R3,cross(L2,L3));
D12 = dot(R1,cross(L1,L3));
D22 = dot(R2,cross(L1,L3));
D32 = dot(R3,cross(L1,L3));
D13 = dot(R1,cross(L1,L2));
D23 = dot(R2,cross(L1,L2));
D33 = dot(R3,cross(L1,L2));

A = 1/D0 * (-tau3/tau13 *D12 + D22 + tau1/tau13 * D32);
B = 1/(6*D0) * (-(tau13^2 - tau3^2) * tau3/tau13 * D12 + (tau13^2 - tau1^2) * tau1/tau13 * D32);

a = -A^2 - 2*A * dot(L2,R2) - norm(R2); 


