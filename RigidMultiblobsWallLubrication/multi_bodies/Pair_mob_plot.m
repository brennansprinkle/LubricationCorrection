clc
clear all
close all

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

cols = lines(7);
cols(1,:) = 0*cols(1,:);

cd '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies'

for k = [1 2 3 4]
    if k == 1
        A = dlmread('./Pair_Mob_Data/mob_coefs.txt');
    elseif k == 2
        A = dlmread('./examples/Tetra_Multiblobs/Pair_Mob_Data/mob_coefs_642.txt');
    elseif k == 3
        A = dlmread('./examples/Tetra_Multiblobs/Pair_Mob_Data/mob_coefs_162.txt');
    elseif k == 4
        A = dlmread('./examples/Tetra_Multiblobs/Pair_Mob_Data/mob_coefs_12.txt');
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
    
    if(k > 1)
        sym = ':';
    else
        sym = '-';
    end
    
    subplot(2,5,p)
    plot(x,c_11(:,p),sym,'linewidth',3,'color',cols(k,:))
    hold all
    title(['self: ' XY ABC],'fontsize',30)
    xlabel('$$\epsilon$$')
    set(gca,'linewidth',3,'fontsize',30)
    subplot(2,5,p+5)
    plot(x,c_12(:,p),sym,'linewidth',3,'color',cols(k,:))
    hold all
    title(['cross: ' XY ABC],'fontsize',30)
    xlabel('$$\epsilon$$')
    set(gca,'linewidth',3,'fontsize',30)
end
end
subplot(2,5,1)
xlabel({'$$\epsilon$$: 2 particle dist.'; 'AND wall dist.'},'fontsize',20)
subplot(2,5,1)
legend('Lub','642','162','12')