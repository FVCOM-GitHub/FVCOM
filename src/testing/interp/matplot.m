
%clear all
close all

n = 121
nsq = n^2;

y = load('fort.12');
z = load('fort.13');
x = load('fort.11');

yz = reshape(y,n,n);
yz =yz';

xz = reshape(x,n,n);
xz=xz';

zz = reshape(z,n,n);
zz=zz';
pcolor(xz,yz,zz)
colorbar
shading interp

hold on

xb = load('fort.14');
yb = load('fort.15');

xb(5) = xb(1);
yb(5) = yb(1);

plot(xb,yb,'k')




