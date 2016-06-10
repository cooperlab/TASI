function [xlarge, ylarge] = enlargepoint(x, y, padsize)
% this function is to include neighbour points of a point (x,y)
% input:
% x: a M*1 or 1*M vector, x coordinate of all point want to enlarge
%
% y: a M*1 or 1*M vector, y coordinate of all point want to enlarge
%
% padsize: padsize is a vector of nonnegative integers that specifies both 
% the amount of padding to add and the dimension along which to add it. 
% For example, a padsize value of [2 3] means add 2 elements of padding 
% along the first dimension (vertical, column) and 3 elements of padding
% along the second dimension (horizontal, row), so the enlarged matrix
% is a 5*7 matrix (2*2+1 3*2+1).
%
% output:
% xlarge: a M*(2*padsize(1)+1)*(2*padsize(2)+1) vector of all the enlarged
% x coordinate
% ylarge: a M*(2*padsize(1)+1)*(2*padsize(2)+1) vector of all the enlarged
% y coordinate
%
% Written by Yue Hou 2016 <lotushouyue@hotmail.com>

[xx, yy] = ndgrid(-padsize(1):padsize(1), -padsize(2):padsize(2));
for i = 1:length(x)
    xpad = padarray(x(i), padsize, x(i));
    ypad = padarray(y(i), padsize, y(i));
    xi = xpad+xx;
    yi = ypad+yy;
    xi = xi(:);
    yi = yi(:);
    if i == 1
        xlarge = xi; ylarge = yi;
    else
        xlarge = cat(1,xlarge,xi); ylarge = cat(1,ylarge,yi);
    end
end