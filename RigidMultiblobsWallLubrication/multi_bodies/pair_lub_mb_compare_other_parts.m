clc
close all

set(0,'defaulttextfontsize',25);
set(0,'defaultaxesfontsize',25);
set(0,'defaultaxeslinewidth',3);
set(0, 'DefaultLineLineWidth',3);
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex');

xyz = {'x','y','z'};
dc = 0
for dist = {'0.025','0.25','1.0','4.0','6.0'} % '
    dc = dc + 1;
    figure(dc)
    clf
    v_lub = dlmread(['tetra_particle_2_velocities_Lub_dist_' dist{1} '.dat']);
    v_mb = dlmread(['./examples/Tetra_Multiblobs/tetra_particle_2_velocities_2562_dist_' dist{1} '.dat']);
    v_mb_12 = dlmread(['./examples/Tetra_Multiblobs/tetra_particle_2_velocities_12_dist_' dist{1} '.dat']);
    
    h_lub = v_lub(:,1);
    v_lub(:,1) = [];
    
    h_mb = v_mb(:,1);
    v_mb(:,1) = [];
    
    h_mb_12 = v_mb_12(:,1);
    v_mb_12(:,1) = [];
    
    
    v_mb = v_mb(:,1:18);
    v_mb_12 = v_mb_12(:,1:18);
    v_lub = v_lub(:,1:18);
    
    p = 0;
    for k = 1:length(v_mb(1,:))
        m = max(abs(v_mb(end,:)));
        if(abs(v_mb(end,k)) > 1e-3*m)
            p = p + 1;
            k
            if(dc == 1)
            subplot(4,4,p)
            else
            subplot(4,4,p)   
            end
            plot(h_mb-1,v_mb(:,k),'-k')
            hold all
            plot(h_mb_12-1,v_mb_12(:,k),':r')
            hold all
            plot(h_lub-1,v_lub(:,k),':b')
            hold off
            comp = mod(k-1,6)+1;
            uw = comp<4;
            comp_x = mod(k-1,3)+1;
            part = floor((k-1)/6.0)+1;
            if(uw)
                title(['$$U^' xyz{comp_x} '_' num2str(part) '$$'])
            else
                title(['$$\omega^' xyz{comp_x} '_' num2str(part) '$$'])
            end
            xlim([0.0 0.51])
            xticks(0.0:0.25:0.5)
            if(p == 4 && dc > 2)
                leg = legend('2562 blobs','12 blobs','Lubrication corrected');
                %set(leg,'fontsize',40,'position',[0.580593659713630   0.049601115909170   0.287982726253463   0.204994791132106])
                set(leg,'fontsize',40,'position',[0.710281159713629   0.014702330822018   0.287982726253464   0.222763224990072])
            end
            if(p == 14)
                xlabel('$$\delta$$')
            end
            
        end
    end
tit = suplabel(['$$\epsilon = $$' dist{1}],'t',[.075 .075 .85 .85]);
set(tit,'fontsize',40)

set(gcf,'units','normalized','outerposition',[0 0 1 1])
drawnow
% print('-dpng','-r90',['/home/bs162/LubricationCode/Tetra_plots_pull_2_in_eps_' dist{1} '.png'])
end
