clc
close all

a = 0.5;
N = 39;
x = (0:2*a:2*a*N)';
y = 0.005*(x-15).^2; %0.25*randn(length(x),1);
z = a+0.2+0*x; %0.25*randn(length(x),1);

plot(x,y)

configs_file = './chain_40.clones';
dlmwrite(configs_file,length(x),'delimiter','\t','precision',5)
dlmwrite(configs_file,[x y z 0*z+1 0*z 0*z 0*z],'-append','delimiter','\t','precision',12)