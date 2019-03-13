k = 0;
set(0,'defaulttextInterpreter','latex')
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
cols = [jet(30)];
for dist = [1.01 1.2 2]
if(dist ~= 2)
    rads = 100:100:1000;
else
    rads = 225:225:2250;
end
for rad = rads
    k = k + 1;
    its = dlmread(['tri_lattice_rad_squared_' num2str(rad) '_dist_' num2str(dist) '.gmres.txt']);
    parts = dlmread(['tri_lattice_rad_squared_' num2str(rad) '_dist_' num2str(dist) '.clones']);
    semilogy(its,'-o','color',cols(k,:),'DisplayName',['num bods = ' num2str(parts(1,1)) ' dist = ' num2str(dist)])
    hold all
end

end
legend show
leg = legend;
set(leg,'fontsize',12)