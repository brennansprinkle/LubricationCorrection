clc
clear all
%close all
clf

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

cols = lines(7);
cols(1,:) = 0*cols(1,:);

cd '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies'

for k = [1]
    if k == 1
        A = dlmread('./Pair_Mob_Data/FarFEO_mob_coefs.txt'); %farFW %FarFEO
    elseif k == 2
        A = dlmread('./examples/Tetra_Multiblobs/Pair_Mob_Data/FarFEO_mob_coefs_2562.txt');
    elseif k == 3
        A = dlmread('./examples/Tetra_Multiblobs/Pair_Mob_Data/FarFEO_mob_coefs_642.txt');
    elseif k == 4
        A = dlmread('./examples/Tetra_Multiblobs/Pair_Mob_Data/FarFEO_mob_coefs_162.txt');
    end
c_11 = A(1:2:end,:);
c_12 = A(2:2:end,:);

if(k == 1)
x = c_11(:,1);
end
c_11(:,1) = [];
c_12(:,1) = [];

[m,n] = size(c_11);

for p = 1:n
    if(any(p == [2 3 5]))
        XY = 'Y';
    else
        XY = 'X';
    end
    
    if(any(p == [1 2]))
        ABC = 'a';
    elseif(p == 3)
        ABC = 'b';
    else
        ABC = 'c';
    end
    
    
    subplot(1,5,p)
    loglog(x,c_11(:,p),'linewidth',3,'color',cols(k,:))
    hold all
    set(gca,'linewidth',3,'fontsize',20)
    title([XY ABC],'fontsize',30)
    xlabel('$$\epsilon$$')

%     subplot(2,5,p+5)
%     plot(x,c_12(:,p),'linewidth',3,'color',cols(k,:))
%     hold all
%     %title(['cross: ' XY ABC],'fontsize',30)
%     xlabel('$$\epsilon$$')
%     set(gca,'linewidth',3,'fontsize',30)
end
end

A = dlmread('Resistance_Coefs/mob_scalars_wall_MB_2562.txt'); %2562
x = A(1:end,1);
A(:,1) = [];
self = A(1:end,:);

for p = 1:n
    if(any(p == [2 3 5]))
        XY = 'Y';
    else
        XY = 'X';
    end
    
    if(any(p == [1 2]))
        ABC = 'a';
    elseif(p == 3)
        ABC = 'b';
    else
        ABC = 'c';
    end
    
    subplot(1,5,p)
    hold all
    loglog(x,self(:,p),':','linewidth',3,'color','r')
    hold all
    %title(['self: ' XY ABC],'fontsize',30)
    xlabel('$$\epsilon$$')
    set(gca,'linewidth',3,'fontsize',25)
end

subplot(1,5,1)
xlabel('$$\epsilon$$: wall distance')
subplot(1,5,1)
legend('Lub','2562')
