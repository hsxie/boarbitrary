% 2022-01-14 16:27 Hua-sheng XIE, huashengxie@gmail.com
% 25-01-02 08:46 update
function G=funGnmp(a,n,m,p,icase)
% p=0,1,2
% imethod=1,2
% global zmx;
nj=length(a); G=0.*a;
zmx=500e0;

for j=1:nj
    % a(j)=max(a(j),1e-3);
    % zmx=max(1e3,1e3/a(j)); % to update
    if(icase==1)
        f=@(x) x.^p.*((1+1i*x)./(1-1i*x)).^m./(1+x.^2).^2.*(besselj(n,a(j).*x)).^2;
        G(j)=integral(f,0,zmx);
    elseif(icase==2)
        f=@(x) x.^p.*((1+1i*x)./(1-1i*x)).^m./(1+x.^2).^2.*besselj(n,a(j).*x).*(...
            besselj(n-1,a(j).*x)-besselj(n+1,a(j).*x))*0.5;
        G(j)=integral(f,0,zmx);
    elseif(icase==3)
        f=@(x) x.^p.*((1+1i*x)./(1-1i*x)).^m./(1+x.^2).^2.*(...
            besselj(n-1,a(j).*x)-besselj(n+1,a(j).*x)).^2*0.25;
        G(j)=integral(f,0,zmx);
    end
end
end