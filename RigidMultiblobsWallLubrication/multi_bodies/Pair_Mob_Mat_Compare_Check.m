clc
clear all
close all

set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

cols = lines(7);
cols(1,:) = 0*cols(1,:);

A = dlmread(['./Pair_Mob_Data/mob_matrix_2562.txt']); 
x = A(:,1);
A(:,1) = [];

syms = {'.','.','.','.'};

cd '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/multi_bodies'
for k = [1 2 3]
    if k == 1
        B = dlmread(['./Pair_Mob_Data/mob_matrix.txt']);
    elseif k == 2
        B = dlmread(['./Pair_Mob_Data/mob_matrix_2562.txt']);
    elseif k == 3
        B = dlmread(['./Pair_Mob_Data/mob_matrix_PC_test.txt']);
        x = B(:,1);
    end


B(:,1) = [];

[m,n] = size(B);

for p = 1:m
    p
    M_mb = reshape(B(p,:),12,12);
    h(k) = plot(x(p),norm(M_mb),syms{k},'markersize',10,'color',cols(k,:));
    hold all
    axis tight
    set(gca,'linewidth',3,'fontsize',20)
    xlabel('$$\epsilon$$')
end


end
figure(1)
xlabel({'$$\epsilon$$: wall distance,';' and particle distance'})
title('Error in mobility: $$||M_{Lub,642,12} - M_{2562} ||_{2}$$')
leg = legend(h,'Lub','642','Lub extended')
set(leg,'fontsize',40)
xlim([0 6])

