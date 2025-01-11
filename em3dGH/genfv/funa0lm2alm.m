% Hua-sheng XIE, huashengxie@gmail.com, 2024-12-13 14:31
% Hermite coefficient a0lm -> alm
% 24-12-31 09:04 update for GPDF-Hermite
function alm=funa0lm2alm(a0lm)
% clear;clc; close all;
% % 
% a0lm=zeros(5,6);
% a0lm(1,1)=1; a0lm(1,2)=2; a0lm(2,2)=-2.4;
% a0lm(4,3)=1; a0lm(3,5)=-1.4;

[lmax,mmax]=size(a0lm);
nmax=max(lmax,mmax);

% dz=.0; dx=0.0; Lz=4.0; Lx=2.0;
% 
% frhol=@(z,l) Lz^2./(Lz^2+z.^2).*((Lz+1i*z)./(Lz-1i*z)).^l;
% fum=@(x,m) 1/sqrt(2^m*factorial(m)*sqrt(pi))*hermiteH0(m, ...
%     sqrt(2)*(x-dx)/Lx).*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);
% 
% fgzl=@(z,l) Lz^2./(Lz^2+z.^2).*((Lz+1i*z)./(Lz-1i*z)).^l;
% fgxm=@(x,m) ((x-dx)/Lx).^m.*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);

cHn0=zeros(nmax,nmax);
cHn0(1,1)=1; % H_0(x)=1
if(nmax>=2)
    cHn0(2,2)=2; % H_1(x) = 2x
end
% a_{n+1,k}=-a_{n,k+1} (k=0), a_{n+1,k}=2*a_{n,k-1}-(k+1)*a_{n,k+1} (k>0)
% H_n(x) = 2x H_{n-1}(x) - 2(n-1)H_{n-2}(x)
for n=2:(nmax-1)
    cHn0(n+1,0+1)=-cHn0(n, 2); % 2x * H_{n-1}(x)
    for k=1:1:n
        cHn0(n+1,k+1)=cHn0(n+1,k+1)+2*cHn0(n, k)-2*(n-1)*cHn0(n-1,k+1);
    end
end

cHn=0.*cHn0;
for n=1:nmax
    cHn(n,:)=cHn0(n,:)/sqrt(2^(n-1)*factorial((n-1))*sqrt(pi));
    for k = 1:1:n
        cHn(n,k)=cHn0(n,k)*1/sqrt(2^((n-1)-(k-1))*factorial((n-1))*sqrt(pi));
    end
end

Nx=mmax; Nz=lmax;
alm=0.*a0lm;
for jz=1:1:Nz
    for jx=1:1:Nx
        % l=jz-1; m=jx-1;
        % alm=alm+a0lm(jz,jx)*(cHn(jz,1:lmax).'*cHn(jx,1:mmax)); % cHn()
        tmp=zeros(1,lmax); tmp(1,jz)=1;
        alm=alm+a0lm(jz,jx)*(tmp.'*cHn(jx,1:mmax)); % cHn()
    end
end

% [xx,zz]=ndgrid(0:0.1:5,-5:0.1:5); fxz0=0.*xx; fxz=0.*xx;
% for jz=1:1:Nz
%     for jx=1:1:Nx
%         l=jz-1; m=jx-1;
%         fxz0=fxz0+a0lm(jz,jx).*frhol(zz,l).*fum(xx,m);
%         fxz=fxz+alm(jz,jx).*fgzl(zz,l).*fgxm(xx,m);
%     end
% end
% 
% subplot(131);
% surf(zz,xx,real(fxz0));
% subplot(132);
% surf(zz,xx,real(fxz));
% subplot(133);
% surf(zz,xx,real(fxz0-fxz));
end