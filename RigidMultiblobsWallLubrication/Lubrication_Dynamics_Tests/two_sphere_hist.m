clc
close all

cd '/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests'
 
A = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/Lub_Spring_Test_MCMC_no_rejection_dt_0.01.Two_Sphere_MCMC_1.config');
A(1:3:end,:) = [];

% B = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/many_bodyMCMC/examples/Two_Sphere_Lub/run_TS_Lub_grav_2.Two_Sphere.config');
% B(1:3:end,:) = [];
% 
% C = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/Mob_Spring_Test_MCMC_no_rejection.Two_Sphere_MCMC_1.config');
% C(1:3:end,:) = [];

s1 = A(1:2:end,1:3);
s2 = A(2:2:end,1:3);

r = sqrt(sum((s1-s2).^2,2));
r(1:10000) = [];
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
r(1:10000) = [];
[h,b] = hist(r,100);
h = h./trapz(b,h);
plot(b,h)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


B = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/many_bodyMCMC/examples/Two_Sphere_Lub/run_TS_Lub_grav_2.Two_Sphere.config');
B(1:3:end,:) = [];

% A = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/Lub_Spring_Test_MCMC_with_rejection.Two_Sphere_MCMC_1.config');
% A(1:3:end,:) = [];


figure(2)
s1 = A(1:2:end,1:3);
s2 = A(2:2:end,1:3);

height = [s1(1000:end,3); s2(1000:end,3)];
[h,b] = hist(height,100);
h = h./trapz(b,h);
plot(b,h)

hold all
s1 = B(1:2:end,1:3);
s2 = B(2:2:end,1:3);

height = [s1(1000:end,3); s2(1000:end,3)];
[h,b] = hist(height,100);
h = h./trapz(b,h);
plot(b,h)

hold all
s1 = C(1:2:end,1:3);
s2 = C(2:2:end,1:3);

height = [s1(1000:end,3); s2(1000:end,3)];
[h,b] = hist(height,100);
h = h./trapz(b,h);
plot(b,h)
hold all
