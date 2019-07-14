clc
close all

bump = @(x) exp(1).*exp(1./(x.^2-1));
cut = @(r,ep) 1.0*(r<=0) + bump((1/ep)*r).*((r > 0) & (r < ep));

delta = 1e-3;
debye = 1e-2;

r = -debye:0.001*debye:debye;

f = @(x) (1./(x+eps)).*(x>0.1*debye) + (1./(0.1*debye+eps)).*(x<0.1*debye);

g = f(0.101*debye).*cut(r-0.1*debye,2*delta) + f(r).*(1-cut(r-0.1*debye,2*delta));
plot(r,g)
hold all
plot(r,f(r))