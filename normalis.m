function [xN,a,b]= normalis(xx,rb,ra)
a = min(xx(:));
b = max(xx(:));
xN = (((ra-rb) * (xx - a)) / (b - a)) + rb;
 %x=(y)*(b - a)/5+a;
