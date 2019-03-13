cd '/home/bs162/LubricationCode/RigidMultiblobsWallLubrication/GMRES_tests/GMRES_test_configs'
close all

for nbods = [2500 5000 10000 15000 20000] %2500 5000 10000 15000 20000

A = dlmread(['np' num2str(nbods) 'h3.particles']);
B = A;
B(1:(nbods+1):end,:) = [];
xyz = B(end-nbods+1:end,:);

x = xyz(:,1);
y = xyz(:,2);
zz = xyz(:,3);

xy = [x y];
sum(pdist(xy)<2.00011)

%f_name = '/home/bs162/Sedimentation/RigidMultiblobsWall/multi_bodies/examples/sedimenting_spheres/no_steric.clones';
f_name = ['np' num2str(nbods) 'h3.clones'];
dlmwrite(f_name,nbods,'delimiter','\t','precision',5)
dlmwrite(f_name,[x y zz 0*x+1 0*x 0*x 0*x],'-append','delimiter','\t','precision',12)
end
