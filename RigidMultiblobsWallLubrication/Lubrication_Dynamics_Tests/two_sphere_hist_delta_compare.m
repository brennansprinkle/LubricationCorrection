% clc
% close all

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

delta = {'1e1','1e2','1e4'};


for k = 1%:length(delta)
%     switch k
%         case 1
%            f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_' delta{k} '_overlap_no_tanh_cut_1e1_dt_0.02.Two_Sphere_Overlap.config']; % 
%         case 2
%            f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_overlap_tanh_cut_1e1_dt_0.02.Two_Sphere_Overlap.config']; % 
%         case 3
%            f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_' delta{k} '_overlap_no_tanh_cut_1e1_dt_0.02.Two_Sphere_Overlap.config']; %
%     end


        f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
                  'Lub_Spring_Test_delta_1e2_overlap_split_RFD_fix_fmin_mob_dt_0.02.Reject_Test_2.config']; %
%         f_name = ['/home/bs162/LubricationCorrection/RigidMultiblobsWallLubrication/Lubrication_Dynamics_Tests/data/' ...
%                   'Lub_Spring_Test_delta_1e2_overlap_dense_center_RFD_fix_fmin_mob_dt_0.02.Reject_Test_2.config']; %

        A = dlmread(f_name);
        A(1:3:end,:) = [];
        
        figure(1)
        s1 = A(1:2:end,1:3);
        s2 = A(2:2:end,1:3);

        skip = floor(100/str2num('0.02'));
        
        r = sqrt(sum((s1-s2).^2,2));
        r(1:skip) = [];
        [h,b] = hist(r,b);
        h = h./trapz(b,h);
        disp_name = ['$$\delta = \frac{1}{' delta{k} '}$$, dt = 0.02'];
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

% 
% n_bods = 2;
% N = length(A)/n_bods;
% n_steps = 1000;
% n_start = 18300;
% xl = min(A(n_start:(n_start+n_steps),1));
% xr = max(A(n_start:(n_start+n_steps),1));
% yl = min(A(n_start:(n_start+n_steps),2));
% yr = max(A(n_start:(n_start+n_steps),2));
% zl = 0;
% zr = max(A(n_start:(n_start+n_steps),3));
% for i = n_start:(n_start+n_steps)
%     clf
%     x = A((i-1)*n_bods+1:i*n_bods,1);
%     y = A((i-1)*n_bods+1:i*n_bods,2);
%     z = A((i-1)*n_bods+1:i*n_bods,3);
%     for j = 1:length(x)
%         h = surface(x(j)+a*sx,y(j)+a*sy,z(j)+a*sz,'facecolor','r','edgecolor','none');
%         set(h,'FaceLighting','gouraud', ...
%         'AmbientStrength',0.5, ...
%         'DiffuseStrength',0.3, ... 
%         'Clipping','off',...
%         'BackFaceLighting','lit', ...
%         'SpecularStrength',0.3, ...
%         'SpecularColorReflectance',0.7, ...
%         'SpecularExponent',1)
%         daspect([1 1 1])
%         view([29 46])
%         xlim([xl xr])
%         ylim([yl yr])
%         zlim([0 zr])
%         hold all
%     end
%     l1 = light('Position',[15 15 100], ...
%     'Style','local', ...
%     'Color',1*[1 1 1]); 
%     drawnow
%     
%     hold off
% end



