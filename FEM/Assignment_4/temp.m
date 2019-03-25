clear all
close all
clear
clc

alpha = 0.001;
delx = 0.01;
delt= 0.04;

tstart = 0;
tend = 300;

xstart = 0;
xend = 1;

timevalues = [tstart:delt:tend]';
spacenodes = [xstart:delx:xend]';
[rowt colt] = size(timevalues);
[rowx colx] = size(spacenodes);

TemperatureMatrix = 20*ones(rowt,rowx);

TemperatureMatrix( : , 1 ) = 100;
TemperatureMatrix( : , rowx ) = 40;

Courant = alpha*delt/(delx^2);
T = TemperatureMatrix;

for n = 1:rowt-1
    for j = 2:rowx-1
      T(n+1,j) = T(n,j) + Courant*(T(n,j-1)-2*T(n,j)+T(n,j+1));
    end
end

figure(1);hold;

for s=0:5:round(tend)
   row = round(s*(rowt-1)/tend)+1;
   plot(spacenodes, T(row,:));
end
