function [x ft] = EProjSimplex_new(v,c)

%
%% Problem
%
%  min  1/2 || x - v||^2
%  s.t. x>=0, 1'x=c
%

if nargin < 2
    c = 1;
end;
k=1;
ft=1;
n = length(v);

v0 = v-mean(v) + c/n;
%vmax = max(v0);
vmin = min(v0);
if vmin < 0
    f = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        v1 = v0 - lambda_m;
        posidx = v1>0;
        npos = sum(posidx);
        g = -npos;
        f = sum(v1(posidx)) - k;
        lambda_m = lambda_m - f/g;
        ft=ft+1;
        if ft > 100
            x = max(v1,0);
            break;
        end;
    end;
    x = max(v1,0);

else
    x = v0;
end;