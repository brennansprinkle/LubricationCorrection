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
    if k == 2
        A = dlmread(['./Pair_Mob_Data/mob_eigs.txt']); %farFW %FarFEO
    elseif k == 1
        A = dlmread(['./Pair_Mob_Data/mob_eigs_2562.txt']);
    elseif k == 3
        A = dlmread(['./Pair_Mob_Data/mob_eigs_642.txt']);
    elseif k == 4
        A = dlmread(['./Pair_Mob_Data/mob_eigs_12.txt']);
    end

if(k == 1)
x = A(:,1);
end
A(:,1) = [];

[m,n] = size(A);

for p = 1:n
    if(k > 1)
        sym = ':';
    else
        sym = '-';
    end
    subplot(3,4,p)
    plot(x,A(:,p),sym,'linewidth',3,'color',cols(k,:))
    hold all
    axis tight
    set(gca,'linewidth',3,'fontsize',20)
    xlabel('$$\epsilon$$')
end

figure(2)
if(k == 1)
    ref = A;
else
for p = 1:n
    if(k > 1)
        sym = ':';
    else
        sym = '-';
    end
    subplot(3,4,p)
    plot(x,A(:,p)-ref(:,p),sym,'linewidth',3,'color',cols(k,:))
    hold all
    axis tight
    set(gca,'linewidth',3,'fontsize',20)
    xlabel('$$\epsilon$$')
end
end
end
figure(1)
subplot(3,4,1)
xlabel({'$$\epsilon$$: wall distance,';' and particle distance'})
subplot(3,4,1)
legend('2562','Lub','642','12')
figure(2)
subplot(3,4,1)
xlabel({'$$\epsilon$$: wall distance,';' and particle distance'})
subplot(3,4,1)
legend('Lub','642','12')

% A = dlmread('Resistance_Coefs/mob_scalars_MB_2562.txt'); %2562
% x = A(1:2:end,1)-2;
% A(:,1) = [];
% self = A(1:2:end,:);
% self(:,3) = -self(:,3);
% cross = A(2:2:end,:);
% cross(:,3) = -cross(:,3);
% 
% for p = 1:n
%     if(any(p == [2 3 5]))
%         XY = 'Y';
%     else
%         XY = 'X';
%     end
%     
%     if(any(p == [1 2]))
%         ABC = 'a';
%     elseif(p == 3)
%         ABC = 'b';
%     else
%         ABC = 'c';
%     end
%     
%     subplot(2,5,p)
%     plot(x,self(:,p),'linewidth',3,'color','b')
%     hold all
%     %title(['self: ' XY ABC],'fontsize',30)
%     xlabel('$$\epsilon$$')
%     set(gca,'linewidth',3,'fontsize',30)
%     subplot(2,5,p+5)
%     plot(x,cross(:,p),'linewidth',3,'color','b')
%     hold all
%     %title(['cross: ' XY ABC],'fontsize',30)
%     xlabel('$$\epsilon$$')
%     set(gca,'linewidth',3,'fontsize',30)
% end
