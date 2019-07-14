clc
%close all

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

% D = dlmread('/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/many_bodyMCMC/examples/Two_Sphere_Lub/run_TS_Lub_grav_2.Two_Sphere.config');
% D(1:3:end,:) = [];

b = linspace(0,6,600);
db = (b(2)-b(1));
b = b + 0.5*db;

bh = linspace(0,5,250);
dbh = (bh(2)-bh(1));
bh = bh + 0.5*dbh;

figure(1)
s1 = D(1:2:end,1:3);
s2 = D(2:2:end,1:3);

r = sqrt(sum((s1-s2).^2,2));
r(1:1000) = [];
[h,b] = hist(r,b);
h = h./trapz(b,h);
plot(b,h,'k','DisplayName','MCMC: hard sphere')
hold all
ylim([0 1.2*max(h)])

figure(2)
h1 = D(1:2:end,3);
h2 = D(2:2:end,3);

height = [h1(20000:end); h1(20000:end)];
[h,bh] = hist(height,bh);
h = h./trapz(bh,h);
plot(bh,h,'k','DisplayName','MCMC: hard sphere')
hold all
ylim([0 1.2*max(h)])

drawnow

delta = {'1e2','1e4'};
dts = {{'0.01','0.02','0.04','0.1'},{'0.005','0.01','0.02'}};

% delta = {'1e2'};
% dts = {{'0.01','0.02','0.04'}};

% delta = {'1e4'};
% dts = {{'0.0025','0.005','0.01','0.02'}};

delta = {'1e1'};
dts = {{'0.16','0.08','0.04','0.02'}}; %,'0.01'


for k = 1 %1:length(delta)
    delta{k}
    for p = 4 %1:length(dts{k})
        dts{k}{p}
        %f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
        %          'Lub_Spring_Test_collision_dt_' dts{k}{p} '.Two_Sphere_Overlap.config']; %
%         f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                  'Lub_Spring_Test_delta_overlap_tanh_' delta{k} '_dt_' dts{k}{p} '.Two_Sphere_Overlap.config']; %
%         f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                  'Lub_Spring_Test_delta_overlap_tanh_cut_1e1_dt_' dts{k}{p} '.Two_Sphere_Overlap.config']; %
%         f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_' delta{k} '_overlap_no_tanh_cut_1e2_dt_' dts{k}{p} '.Two_Sphere_Overlap.config']; %
%         f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_' delta{k} '_overlap_no_tanh_cut_1e1_dt_' dts{k}{p} '.Two_Sphere_Overlap.config']; %
%         f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_' delta{k} '_overlap_no_tanh_cut_1e1_dt_' dts{k}{p} '.Two_Sphere_Overlap.config']; %
%         f_name = ['/fhd/bsprinkle/Lubrication_Two_Sphere/eps_clamp_delta_1e2_dt_0.02.Reject_Test_2.config']; %
          f_name = ['/fhd/bsprinkle/Lubrication_Two_Sphere/Trapz_eps_clamp_delta_1e2_dt_0.02.Reject_Test_2.config']; %
        A = dlmread(f_name);
        A(1:3:end,:) = [];
        
        figure(1)
        s1 = A(1:2:end,1:3);
        s2 = A(2:2:end,1:3);

        skip = floor(100/str2num(dts{k}{p}));
        
        r = sqrt(sum((s1-s2).^2,2));
        r(1:skip) = [];
        [h,b] = hist(r,b);
        h = h./trapz(b,h);
        disp_name = ['$$\delta = \frac{1}{' delta{k} '}$$, dt = ' dts{k}{p}];
        plot(b,h,'DisplayName',disp_name)

        hold all
        figure(2)
        h1 = A(1:2:end,3);
        h2 = A(2:2:end,3);

        height = [h1(skip:end); h1(skip:end)];
        [h,bh] = hist(height,bh);
        h = h./trapz(bh,h);
        plot(bh,h,'DisplayName',disp_name)
        hold all

        drawnow
        legend show
        leg = legend;
        set(leg,'fontsize',25)
    end
end






