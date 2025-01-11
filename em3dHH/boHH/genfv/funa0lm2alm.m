% Hua-sheng XIE, huashengxie@gmail.com, 2024-12-13 14:31
% Hermite coefficient a0lm -> alm
function alm=funa0lm2alm(a0lm)
% clear;clc; close all;
% % 
% a0lm=zeros(5,6);
% a0lm(1,1)=1; a0lm(1,2)=2; a0lm(2,2)=-2.4;
% a0lm(4,3)=1; a0lm(3,5)=-1.4;

[lmax,mmax]=size(a0lm);
nmax=max(lmax,mmax);

% dz=.0; dx=0.0; Lz=1.0; Lx=2.0;
% 
% frhol=@(z,l) 1/sqrt(2^l*factorial(l)*sqrt(pi))*hermiteH0(l, ...
%     sqrt(2)*(z-dz)/Lz).*exp(-(sqrt(2)*(z-dz)/Lz).^2/2);
% fum=@(x,m) 1/sqrt(2^m*factorial(m)*sqrt(pi))*hermiteH0(m, ...
%     sqrt(2)*(x-dx)/Lx).*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);
% 
% fgzl=@(z,l) ((z-dz)/Lz).^l.*exp(-(sqrt(2)*(z-dz)/Lz).^2/2);
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
        alm=alm+a0lm(jz,jx)*(cHn(jz,1:lmax).'*cHn(jx,1:mmax)); % cHn()
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
% surf(zz,xx,fxz0);
% subplot(132);
% surf(zz,xx,fxz);
% subplot(133);
% surf(zz,xx,fxz0-fxz);
end