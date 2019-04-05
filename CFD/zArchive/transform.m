clear all
clc

syms y x
syms alpha
syms beta
syms h
syms neta
h=1;
alpha=0.5;
beta=1.3;
f = alpha + (1-alpha)*(log((beta+y*((1+2*alpha)/h)-2*alpha )/(beta-y*((1+2*alpha)/h)+2*alpha )))/log((beta+1)/(beta-1));
g = finverse(f,y);
y = subs(g,y,neta);
dy = diff(y,neta);

netax = 0:0.2:1.0;
temp = [0 0 0 0 0 0 0];
a = eval(subs(dy,neta,netax));
[row col]=size(a);
gridMatrix = ones(row,col+1);

for i = 1:col+1
    if i==1
        gridMatrix(1)=0;
    else
        gridMatrix(i)=gridMatrix(i-1)+a(i-1);
    end
end

disp(gridMatrix);
figure
scatter(gridMatrix,temp)
