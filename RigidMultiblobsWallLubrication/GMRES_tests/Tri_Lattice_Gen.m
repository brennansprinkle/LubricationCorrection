cd '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/GMRES_tests/GMRES_test_configs'
close all

np = 1001;

ep = 7;
r=sqrt(3)/2;
[x,y] = meshgrid(0:1:np);
n = size(x,1);
x = r * x;
y = y + repmat([0 0.5],[n,n/2]);

x = ep*x(:);
y = ep*y(:);
x = x - max(x)/2;
y = y - max(y)/2;

xref = x;
yref = y;

for cut = 400:400:4000

x = xref;
y = yref; 

R2 = x.^2 + y.^2;
x(R2 >= cut) = [];
y(R2 >= cut) = [];

clf
scatter(x,y,'markerfacecolor','k')
daspect([1 1 1])
drawnow
length(x)

f_name = ['randn_tri_lattice_rad_squared_' num2str(cut) '_dist_' num2str(ep-1) '.clones'];
dlmwrite(f_name,length(x),'delimiter','\t','precision',5)
dlmwrite(f_name,[x y abs(randn(size(x)))+1 0*x+1 0*x 0*x 0*x],'-append','delimiter','\t','precision',12) %ep-1
end

% xy = [x y];
% min(pdist(xy))