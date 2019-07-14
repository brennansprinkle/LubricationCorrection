clc
close all

A = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/many_bodyMCMC/examples/Two_Sphere_Lub/run_TS_Lub_soft.Two_Sphere.config');
A(1:3:end,:) = [];

B = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/many_bodyMCMC/examples/Two_Sphere_Lub/run_TS_Lub_grav_2.Two_Sphere.config');
B(1:3:end,:) = [];

C = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/many_bodyMCMC/examples/Two_Sphere_Lub/run_TS_Lub_soft_debye_0.01a.Two_Sphere.config');
C(1:3:end,:) = [];

s1 = A(1:2:end,1:3);
s2 = A(2:2:end,1:3);

r = sqrt(sum((s1-s2).^2,2));
r(1:1000) = [];
[h,b] = hist(r,100);
h = h./trapz(b,h);
plot(b,h)

hold all
s1 = B(1:2:end,1:3);
s2 = B(2:2:end,1:3);

r = sqrt(sum((s1-s2).^2,2));
r(1:1000) = [];
[h,b] = hist(r,100);
h = h./trapz(b,h);
plot(b,h)

hold all
s1 = C(1:2:end,1:3);
s2 = C(2:2:end,1:3);

r = sqrt(sum((s1-s2).^2,2));
r(1:1000) = [];
[h,b] = hist(r,100);
h = h./trapz(b,h);
plot(b,h)