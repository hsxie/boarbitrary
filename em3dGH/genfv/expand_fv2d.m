% Hua-sheng XIE, huashengxie@gmail.com, 2024-12-02 13:33
% f(z,x)=Wz(z)*sum{a_{ml}*rho_m(x)*h_l(z)} of GPDF rho_m & Hermite h_l to
% expansion distribution function f_0(vpara,vperp)
% can run
% 24-12-31 10:18 update
close all; clear; clc;

sid=1; % species number index
load(['fvdata',num2str(sid),'.mat']);
fv=fvdat.fv;
vz=fvdat.vz;
vx=fvdat.vx;
dvz=fvdat.dvz;
dvx=fvdat.dvx;

% dz=0.0*fvdat.vdz; dx=0.0*fvdat.vdx; Lz=1.8*fvdat.vtz; Lx=2*fvdat.vtx;
dz=0.0*fvdat.vdz; dx=0.0*fvdat.vdx; Lz=1.15*fvdat.vtz; Lx=1.15*fvdat.vtx; % dx=0.1*Lz;

Nx=2^5; Nz=2^4;


% frhol=@(z,l) 1/sqrt(2^l*factorial(l)*sqrt(pi))*hermiteH0(l, ...
%     sqrt(2)*(z-dz)/Lz).*exp(-(sqrt(2)*(z-dz)/Lz).^2/2);
fum=@(x,m) 1/sqrt(2^m*factorial(m)*sqrt(pi))*hermiteH0(m, ...
    sqrt(2)*(x-dx)/Lx).*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);
% fgzl=@(z,l)((z-dz)/Lz).^l.*exp(-(sqrt(2)*(z-dz)/Lz).^2/2);
fgxm=@(x,m)((x-dx)/Lx).^m.*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);

fWz=@(z) Lz^2./(Lz^2+z.^2);
frhol=@(z,l) ((Lz+1i*z)./(Lz-1i*z)).^l;
% fhm=@(x,m) 1/sqrt(2^m*factorial(m)*sqrt(pi))*hermiteH0(m,x/Lx).*exp(-(x/Lx).^2/2);
% fgm=@(x,m) ((x-dx)/Lx).^m.*exp(-(sqrt(2)*(x-dx)/Lx).^2/2); %.*exp(-(x/Lx).^2/2);
% fum=@(x,m) 1/sqrt(2^m*factorial(m)*sqrt(pi))*hermiteH0(m, ...
%     sqrt(2)*(x-dx)/Lx).*exp(-(sqrt(2)*(x-dx)/Lx).^2/2);

runtime=cputime;

a0lm=zeros(2*Nz+1,Nx+1);

for jz=-Nz:1:Nz
    for jx=0:1:Nx
        l=jz; m=jx;
        a0lm(jz+Nz+1,jx+1)=sum(fv(:,1).*conj(frhol(vz(:,1),l)).*fum(vx(:,1),m)+ ...
            sum(fv(:,2:end).*conj(frhol(vz(:,2:end),l)).*fum(vx(:,2:end),m),2)+ ...
            sum(fv(:,2:end).*conj(frhol(vz(:,2:end),l)).*fum(-vx(:,2:end),m),2) ...
            )*dvx*dvz*sqrt(2)/(Lx*Lz)/pi;
    end
end

%%
bs=dx/Lx;
As=exp(-bs^2)+sqrt(pi)*bs*erfc(-bs);
% cs0=1/(sqrt(pi^3)*Lz*Lx^2*As);
cs0=1/(pi^2*Lz*Lx^2*As);

alm=funa0lm2alm(a0lm/cs0); % 24-12-21 14:58

fxz=0.*vx; zz=vz; xx=vx; fxz1=0.*fxz;
% [xx,zz]=ndgrid(0:0.1:5,-5:0.1:5); fxz=0.*xx;
for jz=-Nz:1:Nz
    for jx=0:1:Nx
        l=jz; m=jx;
        fxz=fxz+a0lm(jz+Nz+1,jx+1)*fWz(zz).*frhol(zz,l).*fum(xx,m);
        fxz1=fxz1+cs0*alm(jz+Nz+1,jx+1)*fWz(zz).*frhol(zz,l).*fgxm(xx,m);
    end
end

runtime=cputime-runtime;

fvc.vdz=dz; fvc.vdr=dx; fvc.vtz=Lz; fvc.vtp=Lx; fvc.alm=alm; fvc.a0lm=a0lm;
save(['fvceff',num2str(sid),'.mat'],'fvc');
runtime=cputime-runtime;

save('fvfitdat.mat');
%%
close all;
subplot(221);
contour(vz,vx,real(fv),(0.01:0.01:1.05)*max(max(real(fv)))); colorbar;
% xlim([-5,5]);ylim([0,5]);
xlabel('v_z'); ylabel('v_x');
title(['fv, Lz=',num2str(Lz),', Lx=',num2str(Lx),', Nz=',num2str(Nz),', Nx=',num2str(Nx)]);
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
print('-dpng',['expand_fv2d_Lz=',num2str(Lz),',Lx=',num2str(Lx),...
    ',Nz=',num2str(Nz),',Nx=',num2str(Nx),'.png']);
%%
figure;
surf(vz,vx,real(fv)-real(fxz)); hold on;
% surf(vz,vx,real(fxz1)-real(fxz)); hold on;
% surf(zz,xx,real(fxz));
%%
figure;
subplot(121);
plot(vz(:,1),real(fv(:,1)),vz(:,1),real(fxz1(:,1)),':','LineWidth',2);
subplot(122);
plot(vz(:,1),real(real(fxz1(:,1))./fv(:,1)),':','LineWidth',2);
ylim([-2,4]);