% Hua-sheng XIE, huashengxie@gmail.com, 2025-01-02 10:34
% f(z,x)=Wz(z)*Wx(x)^2*sum{a_{lm}*rho_l(z)*u_m(x)} of GPDF rho_l & u_m to
% expansion distribution function f_0(vpara,vperp)
% 11:36 use fft
% 20:24 test ok, agree with direct method
close all; clear; clc;

sid=1; % species number index
load(['fvdata',num2str(sid),'.mat']);
fv=fvdat.fv;
vz=fvdat.vz;
vx=fvdat.vx;
dvz=fvdat.dvz;
dvx=fvdat.dvx;

dz=0.0*fvdat.vdz; dx=0.0*fvdat.vdx; Lz=1.15*fvdat.vtz; Lx=1.15*fvdat.vtx;

Nx=2^5; Nz=2^4;

fWz=@(z) Lz^2./(Lz^2+z.^2);
frhol=@(z,l) ((Lz+1i*z)./(Lz-1i*z)).^l;
fWx=@(x) Lx^2./(Lx^2+x.^2);
fum=@(x,m) ((Lx+1i*x)./(Lx-1i*x)).^m;

runtime=cputime;
imethod=2;
if(imethod==1) % use FFT

    t=(-(Nz-1):1:(Nz-1))*pi/Nz; dt=t(2)-t(1); nt=length(t);
    p=(0:1:(Nx-1))*pi/Nx; dp=p(2)-p(1); np=length(p);
    [tt,pp]=ndgrid(t,p);
    % [tt,pp]=meshgrid(t,p);
    
    vzp=Lz*tan(tt/2); vxp=Lx*tan(pp/2);
    gzx=zeros(2*Nz,2*Nx);
    g0=fv./(fWx(vx).^2.*fWz(vz));
    gzx(2:end,(Nx+1):end)=griddata(vz,vx,g0,vzp,vxp);
    gzx(isnan(gzx))=0;
    gzx(2:end,2:Nx)=fliplr(gzx(2:end,(Nx+2):end));

    % G = fftshift(fft2(gzx));
    G = fftshift(fft2(fftshift(gzx)));
    a0lm=zeros(2*Nz+1,2*Nx+1);
    % a0lm(2:end,2:end) = G / (4 * pi^2*(2*Nx)^2*(2*Nz)^2);
    % a0lm(1:(2*Nz),1:(2*Nx)) = G / (4 * pi^2*(2*pi));
    a0lm(1:(2*Nz),1:(2*Nx)) = G / ((2*Nx)*(2*Nz));

else % direct calculate
a0lm=zeros(2*Nz+1,2*Nx+1);

for jz=-Nz:1:Nz
    for jx=-Nx:1:Nx
        l=jz; m=jx;
        a0lm(jz+Nz+1,jx+Nx+1)=sum(fv(:,1)./fWx(vx(:,1)).*conj(frhol(vz(:,1),l)).*conj(fum(vx(:,1),m))+ ...
            sum(fv(:,2:end)./fWx(vx(:,2:end)).*conj(frhol(vz(:,2:end),l)).*conj(fum(vx(:,2:end),m)),2)+ ...
            sum(fv(:,2:end)./fWx(vx(:,2:end)).*conj(frhol(vz(:,2:end),l)).*conj(fum(-vx(:,2:end),m)),2) ...
            )*dvx*dvz/(Lx*Lz*pi^2);
    end
end
end
%%
% close all;
% %     gzx(2:end,2:Nx)=fliplr(gzx(2:end,(Nx+2):end));
% %     surf(gzx);
% subplot(221);surf(real(a0lm));
% subplot(222);surf(real(a0lm1));
% subplot(223);surf(imag(a0lm));
% subplot(224);surf(imag(a0lm1));
%%
cs0=1/(pi^2*Lz*Lx^2);

alm=a0lm/cs0; % 24-12-21 14:58

fxz=0.*vx; zz=vz; xx=vx; fxz1=0.*fxz;
% [xx,zz]=ndgrid(0:0.1:5,-5:0.1:5); fxz=0.*xx;
for jz=-Nz:1:Nz
    for jx=-Nx:1:Nx
        l=jz; m=jx;
        fxz=fxz+a0lm(jz+Nz+1,jx+Nx+1)*fWz(zz).*fWx(xx).^2.*frhol(zz,l).*fum(xx,m);
        % fxz1=fxz1+cs0*alm(jz+Nz+1,jx+Nx+1)*fWz(zz).*fWx(xx).^2.*frhol(zz,l).*fum(xx,m);
    end
end

runtime=cputime-runtime;

fvc.vdz=dz; fvc.vdr=dx; fvc.vtz=Lz; fvc.vtp=Lx; fvc.alm=alm; fvc.a0lm=a0lm;
save(['fvceff',num2str(sid),'.mat'],'fvc');
runtime=cputime-runtime;
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
plot(vz(:,1),real(fv(:,1)),vz(:,1),real(fxz(:,1)),':','LineWidth',2);
subplot(122);
plot(vz(:,1),real(real(fxz(:,1))./fv(:,1)),':','LineWidth',2);
ylim([-2,4]);
