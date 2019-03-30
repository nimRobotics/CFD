from numpy import *
import scipy.linalg
# https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
# http://mathesaurus.sourceforge.net/matlab-numpy.html

a =array([[1.,2.,3.], [4.,5.,6.]])
print(a[0,2]*0.4)
a[:]=55
print(a)

% [row col]=size(a);
% gridMatrix = ones(row,col+1);
% for i = 1:col+1
%     if i==1
%         gridMatrix(1)=0;
%     else
%         gridMatrix(i)=gridMatrix(i-1)+a(i-1);
%     end
% end
