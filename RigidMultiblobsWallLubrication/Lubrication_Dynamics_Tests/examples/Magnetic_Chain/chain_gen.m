clc
close all

a = 0.5;
N = 59;
x = (0:2*a:2*a*N)';
y = 0*x; %0.25*randn(length(x),1);
z = 2*a+0*x; %0.25*randn(length(x),1);
z(z<1) = z(z<1)+1;

configs_file = './chain_y.clones';
dlmwrite(configs_file,length(x),'delimiter','\t','precision',5)
dlmwrite(configs_file,[x y z 0*z+1 0*z 0*z 0*z],'-append','delimiter','\t','precision',12)