% Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, 2018-12-02 01:08
% Adaptive Simpson quadrature formula to calculate the v_perp integral
% for ring beam dispersion relation.
% To do: Speed up using Gauss integral
% 24-12-19 19:06 modified to arbitrary distribution cases
function Cn=funCn(n,a,d,m,p)
na=length(a); Cn=0.*a;
% ym=5.0;
ym=10.0;
reltol=1e-15; abstol=1e-10;
ymin=max(0,d-ym);
% ymin=0;
ymax=ym+d;
if(m>=0)
    for ja=1:na
        f=@(y) 0.25.*(besselj(n-1,a(ja)*y)-...
            besselj(n+1,a(ja)*y)).^2.*exp(-(y-d).^2).*(y-d).^m.*y.^p;
        % Cn(ja)=integral(f,ymin,ymax);
            Cn(ja)=integral(f,ymin,ymax,'RelTol',reltol,'AbsTol',abstol);
        %     Cn(ja)=integral(f,ymin,ymax,'RelTol',reltol);
    end
end