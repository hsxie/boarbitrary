% Hua-sheng XIE, huashengxie@gmail.com, 2024-12-06 20:35
% f(z,x)=sum{a_{lm}*rho_l(z)*u_m(x)} of Hermite rho_l & u_m
% expansion distribution function f_0(vpara,vperp)
% can run
% 24-12-13 16:08 update the vx<0 region
% 24-12-21 14:50 update to obtain a0lm -> alm
close all; clear; clc;

sid=1; % species number index
load(['fvdata',num2str(sid),'.mat']);
fv=fvdat.fv;
vz=fvdat.vz;
vx=fvdat.vx;
dvz=fvdat.dvz;
dvx=fvdat.dvx;

% dz=0.0*fvdat.vdz; dx=0.0*fvdat.vdx; Lz=1.8*fvdat.vtz; Lx=2*fvdat.vtx;
dz=0.0*fvdat.vdz; dx=0.0*fvdat.vdx; Lz=1.0*fvdat.vtz; Lx=1.0*fvdat.vtx;
% dz=0.0*fvdat.vdz; dx=0.0*fvdat.vdx; Lz=0.8*fvdat.vdx; Lx=0.8*fvdat.vdx;

Nx=2^4; Nz=2^4+0;

frhol=@(z,l) 1/sqrt(2^l*factorial(l)*sqrt(pi))*hermiteH0(l, ...
    sqrt(2)*(z-dz)/Lz).*exp(-(sqrt(2)*(z-dz)/Lz).^2/2);
fum=@(x,m) 1/sqrt(2^m*factorial(m)*sqrt(pi))*hermiteH0(m, ...
    sqrt(2)*(x-dx)/Lx).*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);

fgzl=@(z,l)((z-dz)/Lz).^l.*exp(-(sqrt(2)*(z-dz)/Lz).^2/2);
fgxm=@(x,m)((x-dx)/Lx).^m.*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);

runtime=cputime;

a0lm=zeros(Nz+1,Nx+1);

% for jz=0:1:Nz
%     for jx=0:2:Nx
%         l=jz; m=jx;
%         alm(jz+1,jx+1)=sum(fv(:,1).*frhol(vz(:,1),l).*fum(vx(:,1),m)+ ...
%             2*sum(fv(:,2:end).*frhol(vz(:,2:end),l).*fum(vx(:,2:end),m),2))*dvx*dvz*2/(Lx*Lz);
%     end
% end
for jz=0:1:Nz
    for jx=0:1:Nx
        l=jz; m=jx;
        % 24-12-13 16:08 update the vx<0 region
        a0lm(jz+1,jx+1)=sum(fv(:,1).*frhol(vz(:,1),l).*fum(vx(:,1),m)+ ...
            sum(fv(:,2:end).*frhol(vz(:,2:end),l).*fum(vx(:,2:end),m),2)+ ...
            sum(fv(:,2:end).*frhol(vz(:,2:end),l).*fum(-vx(:,2:end),m),2) ...
            )*dvx*dvz*2/(Lx*Lz);
    end
end

%%
fxz=0.*vx; zz=vz; xx=vx;
% [xx,zz]=ndgrid(0:0.1:5,-5:0.1:5); fxz=0.*xx;
for jz=0:1:Nz
    for jx=0:1:Nx
        l=jz; m=jx;
        fxz=fxz+a0lm(jz+1,jx+1).*frhol(zz,l).*fum(xx,m);
    end
end

bs=dx/Lx;
As=exp(-bs^2)+sqrt(pi)*bs*erfc(-bs);
cs0=1/(sqrt(pi^3)*Lz*Lx^2*As);

alm=funa0lm2alm(a0lm/cs0); % 24-12-21 14:58
runtime=cputime-runtime;

fvc.vdz=dz; fvc.vdr=dx; fvc.vtz=Lz; fvc.vtp=Lx; fvc.alm=alm; fvc.a0lm=a0lm;
save(['fvceff',num2str(sid),'.mat'],'fvc');
%%
close all;
subplot(221);
contour(vz,vx,real(fv),(0.01:0.01:1.05)*max(max(real(fv)))); colorbar;
% xlim([-5,5]);ylim([0,5]);
xlabel('v_z'); ylabel('v_x');
title(['fv, Lz=',num2str(Lz,3),', Lx=',num2str(Lx,3),', Nz=',num2str(Nz),', Nx=',num2str(Nx)]);
subplot(222);
contour(zz,xx,real(fxz),(0.01:0.01:1.05)*max(max(real(fv)))); colorbar;
% xlim([-5,5]);ylim([0,5]);
xlabel('v_z'); ylabel('v_x'); title('fv_{fit}');
subplot(223);
surf(vz,vx,real(fv));
% xlim([-5,5]);ylim([0,5]);
xlabel('v_z'); ylabel('v_x'); zlabel('fv');
subplot(224);
surf(zz,xx,real(fxz));
% xlim([-5,5]);ylim([0,5]);
xlabel('v_z'); ylabel('v_x'); zlabel('fv_{fit}');
print('-dpng',['expand_fv2d_Lz=',num2str(Lz,3),',Lx=',num2str(Lx,3),...
    ',Nz=',num2str(Nz),',Nx=',num2str(Nx),'.png']);
%%
figure;
surf(vz,vx,real(fv)-real(fxz)); hold on;
%%
figure;
plot(vz(:,1),real(fv(:,1)),vz(:,1),real(fxz(:,1)),':','LineWidth',2);
% surf(zz,xx,real(fxz));
%%
sumfv=0;
for jz=0:1:Nz
    for jx=0:1:Nx
        l=jz; m=jx;
        f=@(y) exp(-(y-dx/Lx).^2).*(y-dx/Lx).^m.*y;
        Am=integral(f,0,10);
        sumfv=sumfv+2/As*alm(jz+1,jx+1)*funIn(l)*Am;
    end
end

% vtps=sqrt(2*kB*Tps./ms); % perp thermal velocity, note the sqrt(2)
% wcs=B0*qs./ms; % cyclotron frequency
% rhocs=sqrt(kB*Tps./ms)./wcs; % cyclotron radius, 2018-06-13 21:47
% % as=kx*rhocs*sqrt(2);
% as=0.2;
% 
% s11=0; s22=0; s33=0; % 24-12-21 18:45
% Nmax=10;
% for n=-Nmax:1:Nmax
% for jz=0:1:Nz
%     for jx=0:1:Nx
%         l=jz; m=jx;
%         s11=s11+2/As/as^2*alm(jz+1,jx+1)*funIn(l)*n^2*( ...
%             2*funAn(n,as,bs,m+1,0)-m*funAn(n,as,bs,m-1,0));
%         s22=s22+2/As*alm(jz+1,jx+1)*funIn(l)*( ...
%             2*funCn(n,as,bs,m+1,2)-m*funCn(n,as,bs,m-1,2));
%         s33=s33+2/As*alm(jz+1,jx+1)*( bs*(2*funIn(l+1)-l*funIn(l-1))+...
%             (2*funIn(l+2)-l*funIn(l)))*funAn(n,as,bs,m,1);
%     end
% end
% end

save('fvfitdat.mat');
