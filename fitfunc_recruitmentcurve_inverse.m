function x = fitfunc_recruitmentcurve_inverse(y,a,b,c)
x = (-1/b) .* log(a./(y) - 1) + c;
% x = (-1/b) .* log(a./(y-d) - 1) + c;
end