% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2024-12-19 19:52
% In in Z function
function In=funIn(n)
if(mod(n,2)==1)
    In=0;
else
    % In=2.0/sqrt(pi)*gamma(0.5*(n+1));
    In=1.0/sqrt(pi)*gamma(0.5*(n+1));
end
end