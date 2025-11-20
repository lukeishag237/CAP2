function y = fitfunc_recruitmentcurve(x,a,b,c)
y = a./(1+exp(-b.*(x - c)));
% y = a./(1+exp(-b.*(x - c))) + d;

% y = log10(a./(1+exp(-b.*(x - c))) + d);
% y = a./(1+exp(-b.*(x))) + d;
end