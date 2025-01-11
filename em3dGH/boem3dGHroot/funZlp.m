% 2024-12-31 08:24, Huasheng XIE, huashengxie@gmail.com, ENN
% calculate function Zlp for GPDF basis
% 14:10 to check
function Z=funZlp(z,l,p)
Z=0.*z;
if(p==0)
    if(l==0)
        Z=-1i./(1-1i*z);
    elseif(l>0)
        Z=-2i./(1+z.^2).*((1+1i*z)./(1-1i*z)).^l;
    end
elseif(p==1)
    Z=z.*funZlp(z,l,0)-(l==0);
elseif(p==2)
    % Z=z.*funZlp(z,l,1)-sign(l)*(-1)^l*1i;
    Z=z.*funZlp(z,l,1)+sign(l)*(-1)^l*1i;
end
% Z=-Z;
end