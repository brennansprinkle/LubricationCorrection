clc
close all

% a = 0.5;
% N = 59;
% x = (0:2*a:2*a*N)';
% y = 0*x; %0.25*randn(length(x),1);
% z = a+0.2+0*x; %0.25*randn(length(x),1);
% 
% plot(x,y)
% 
% configs_file = './chain_60.clones';
% dlmwrite(configs_file,length(x),'delimiter','\t','precision',5)
% dlmwrite(configs_file,[x y z 0*z+1 0*z 0*z 0*z],'-append','delimiter','\t','precision',12)

a = 0.5;
%N = 28; 
N = 25;
x = (2.1*a:2.1*a:2.1*a*N)';
r = 3*2.1*a;
% r = 1*2.1*a;
theta = linspace(pi/2,3*pi/2,10);
% theta = linspace(pi/2,3*pi/2,4);
cx = r*cos(theta');
cy = r*sin(theta');
xp = [flipud(x); cx;  x];
yp = [0*x-r; flipud(cy); 0*x+r];
zp = a+0.2+0*xp;

for k = 1:length(xp)
    if(k < length(xp))
        dist = sqrt((xp(k+1)-xp(k)).^2 + (yp(k+1)-yp(k)).^2);
        disp(dist)
    end
    circle2(xp(k),yp(k),a);
    hold all
end

plot(xp,yp,'o')
daspect([1 1 1])

configs_file = './hairpin_60.clones';
dlmwrite(configs_file,length(xp),'delimiter','\t','precision',5)
dlmwrite(configs_file,[xp yp zp 0*zp+1 0*zp 0*zp 0*zp],'-append','delimiter','\t','precision',12)

