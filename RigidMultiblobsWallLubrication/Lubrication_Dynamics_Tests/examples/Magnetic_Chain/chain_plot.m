clc
clf

n_bods = 60;

a = 0.5;
[sx, sy, sz] = sphere(20);
figure(1)

%A = dlmread('./data/test.chain.config');
%A = dlmread('./data/test_long.chain_long.config');
% A = dlmread('/fhd/bsprinkle/Mag_Chain/test_theta_hat_eq_2.chain_long_eq.config');
% A = dlmread('/fhd/bsprinkle/Mag_Chain/test_theta_hat_rot_more_bend.chain_long_eq.config');
% A = dlmread('/fhd/bsprinkle/Mag_Chain/test_lasso_no_lub_fs.chain_long_fs.config');
% A = dlmread('/fhd/bsprinkle/Mag_Chain/test_lasso_weak.chain_long_eq.config'); % chain.config
% A = dlmread('/fhd/bsprinkle/Mag_Chain/test_helix_weak.chain_long_eq.config');
% A = dlmread('/fhd/bsprinkle/Mag_Chain/test_lasso_rep_bend.chain_long_eq.config'); %chain_long_eq.config %_bend


A = dlmread('/fhd/bsprinkle/Mag_Chain/lasso_toy.chain_long_eq.config');
%A = dlmread('/fhd/bsprinkle/Mag_Chain/helix_anis_strong.chain_long_eq.config');

N = length(A)/n_bods;
A(1:(n_bods+1):end,:) = [];
dt = 20*0.001;
skip = 4*4; %4*20;

Nhist = 100;
cols = jet(Nhist);
k = 0;

[X, Y] = meshgrid(a*[-30:0.5:160],a*[-30:0.5:160]);

for i = 1:skip:(length(A)/n_bods)
    clf
    k = k+1;
    x = A((i-1)*n_bods+1:i*n_bods,1);
    y = A((i-1)*n_bods+1:i*n_bods,2);
    z = A((i-1)*n_bods+1:i*n_bods,3);
    
    for j = 1:length(x)
        
        h = surface(x(j)+a*sx,y(j)+a*sy,z(j)+a*sz,'facecolor','r','edgecolor','none');
        set(h,'FaceLighting','gouraud', ...
        'AmbientStrength',0.5, ...
        'DiffuseStrength',0.3, ... 
        'Clipping','off',...
        'BackFaceLighting','lit', ...
        'SpecularStrength',0.3, ...
        'SpecularColorReflectance',0.7, ...
        'SpecularExponent',1)
        daspect([1 1 1])
        view([53 49])% view([90 90]) %
        xlim(a*[-30 160])
        ylim(a*[-30 160])
        zlim(a*[0 20])
        hold all
    end
    surface(X,Y,0*X,'facecolor','k','edgecolor','none')
    l1 = light('Position',[15 15 max(z)+100], ...
    'Style','local', ...
    'Color',1*[1 1 1]); 
    title(['t = ' num2str((i-1)*dt)])
    drawnow
    
    hold off
    %print('-dpng',['chain_pngs/lasso_anis_' num2str(k) '.png'],'-r100')
end

% configs_file = './chain_long_fs.clones';
% dlmwrite(configs_file,length(x),'delimiter','\t','precision',5)
% dlmwrite(configs_file,[x y z+100000 0*z+1 0*z 0*z 0*z],'-append','delimiter','\t','precision',12)

% close all
% for t = 0:0.01:10*pi
%     plot3([0 1.43], [0 2.94*sin(2*pi*t)], [0 2.94*cos(2*pi*t)],'-o')
%     axis([-3 10 -5 5 -5 5])
%     daspect([1 1 1])
%     view([53 49])
%     drawnow
%     pause(0.1)
% end