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

    clf
    v_lub = dlmread(['tetra_particle_2_velocities_Lub_epsilon.dat']);
    v_mb = dlmread(['./examples/Tetra_Multiblobs/tetra_particle_2_velocities_2562_epsilon.dat']);
    v_mb_12 = dlmread(['./examples/Tetra_Multiblobs/tetra_particle_2_velocities_12_epsilon.dat']);
    
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
    for k = 1:length(v_lub(1,:))
        m = max(abs(v_lub(end,:)));
        if(abs(v_lub(end,k)) > 1e-4*m)
            p = p + 1;
            k
            figure(1)
            subplot(3,4,p)
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
            xlim([0.0 2])
            xticks(0.0:0.5:2)
            if(p == 4)
                leg = legend('2562 blobs','12 blobs','LubCorr');
                %set(leg,'fontsize',40,'position',[0.580593659713630   0.049601115909170   0.287982726253463   0.204994791132106])
                set(leg,'fontsize',25,'position',[0.710281159713629   0.014702330822018   0.287982726253464   0.222763224990072])
            end
            if(p == 14)
                xlabel('$$\delta$$')
            end
            
            figure(2)
            subplot(3,4,p)
            plot(h_mb-1,abs(v_mb(:,k)-v_lub(:,k)),'-b')
            hold all
            plot(h_mb_12-1,abs(v_mb(:,k)-v_mb_12(:,k)),'-r')
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
            xlim([0.0 2])
            xticks(0.0:0.5:2)
            if(p == 4)
                leg = legend('error in LubCorr','error in 12 blobs');
                %set(leg,'fontsize',40,'position',[0.580593659713630   0.049601115909170   0.287982726253463   0.204994791132106])
                set(leg,'fontsize',25,'position',[0.710281159713629   0.014702330822018   0.287982726253464   0.222763224990072])
            end
            if(p == 14)
                xlabel('$$\delta$$')
            end
            
            
        end
    end


set(gcf,'units','normalized','outerposition',[0 0 1 1])
drawnow

figure(3)
[max_lub, I_lub] = max(abs(v_mb-v_lub),[],2);
[max_mb, I_mb] = max(abs(v_mb-v_mb_12),[],2);

lub_vals = unique(I_lub);
mb_vals = unique(I_mb);

lub_col = [0 0 1; 0 0.8 1];
mb_col = [1 0 0; 1 0.8 0];



subplot(2,1,1)
p = 0;
for k = 1:length(lub_vals)
    p = p + 1;
    v = lub_vals(k);
    
    comp = mod(v-1,6)+1;
    uw = comp<4;
    comp_x = mod(v-1,3)+1;
    part = floor((v-1)/6.0)+1;
    if(uw)
        str_leg = ['$$U^' xyz{comp_x} '_' num2str(part) '$$'];
    else
        str_leg = ['$$\omega^' xyz{comp_x} '_' num2str(part) '$$'];
    end
    
    
    cut = (I_lub == v+0*I_lub);
    x = h_mb-1;
    y = max_lub;
    y(cut) = NaN;
    plot(x,y,'color',lub_col(k,:))
    legn{p} = ['max error in LubCorr, from ' str_leg]
    hold all
end
hold all
for k = 1:length(mb_vals)
    p = p + 1;
    v = mb_vals(k);
    
    comp = mod(v-1,6)+1;
    uw = comp<4;
    comp_x = mod(v-1,3)+1;
    part = floor((v-1)/6.0)+1;
    if(uw)
        str_leg = ['$$U^' xyz{comp_x} '_' num2str(part) '$$'];
    else
        str_leg = ['$$\omega^' xyz{comp_x} '_' num2str(part) '$$'];
    end
    
    cut = (I_mb == v+0*I_mb);
    legn{p} = ['max error in 12 blobs, from ' str_leg]
    x = h_mb-1;
    y = max_mb;
    y(cut) = NaN;
    plot(x,y,'color',mb_col(k,:))
    hold all
end
hold off
leg = legend(legn);
set(leg,'fontsize',40)
xlabel('$$\epsilon$$')
ylabel('Error')

subplot(2,1,2)
plot(h_mb-1,sqrt(sum((v_mb-v_lub).^2,2)),'-b')
hold all
plot(h_mb_12-1,sqrt(sum((v_mb-v_mb_12).^2,2)),'-r')
hold off
leg = legend('$$L^2$$ error in LubCorr','$$L^2$$ error in 12 blobs');
set(leg,'fontsize',40)
xlabel('$$\epsilon$$')
ylabel('Error')