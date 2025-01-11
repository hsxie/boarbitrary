% 2024-12-10 14:32 Hua-sheng XIE, huashengxie@gmail.com, ENN
% Generate distribution function f(vpara,vperp) on grid
% 24-12-22 14:51 add shell and slowing down
% 25-01-03 10:30 shell case in Min15

close all; clear; clc;

% constants
% % Default SI unit
% c2=(2.99792458E8)^2; % speed of ligth c^2
% epsilon0=8.854187817E-12;
% mu0=1/(c2*epsilon0);
mu0=1.2566e-06;
c2=(1.5173e7)^2; % 25-01-03 11:41
epsilon0=1/(c2*mu0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

munit=mp; % set the default normalization for mass

% ms0=5.447e-4; % m_unit
ms0=1.0; % m_unit
Tzs0=2*5e2; % eV
Tps0=2*5e2; % eV
vdsz0=0.0; % c
vdsr0=0.13333; % c
vA=sqrt(c2)/15;

ms=ms0*munit; % * unit mass, -> kg, 18-12-23 22:59
Tzs=Tzs0*qe/kB; % T//, eV -> K (eV -> J * J -> K)
Tps=Tps0*qe/kB; % Tpr, eV -> K
vdsz=vdsz0*sqrt(c2); % vds_z, speed of light c -> m/s
vdsr=vdsr0*sqrt(c2); % vds_r, speed of light c -> m/s

icase=4;

kappa=5.5;
if(icase==2)
    vtzs=sqrt(2*(1-1.5/kappa)*kB*Tzs./ms); % para thermal velocity, note the sqrt(2)
    vtps=sqrt(2*(1-1.5/kappa)*kB*Tps./ms); % perp thermal velocity, note the sqrt(2)
elseif(icase==3)
    vtzs=sqrt(2*(1-0.5/kappa)*kB*Tzs./ms); % para thermal velocity, note the sqrt(2)
    vtps=sqrt(2*(1-1/kappa)*kB*Tps./ms); % perp thermal velocity, note the sqrt(2)
else
    vtzs=sqrt(2*kB*Tzs./ms); % para thermal velocity, note the sqrt(2)
    vtps=sqrt(2*kB*Tps./ms); % perp thermal velocity, note the sqrt(2)
end

sid=1; % species number index

vtx=vtps; vtz=vtzs; vdz=vdsz; vdx=vdsr;

% dvx=0.05*vtx; dvz=0.05*vtz;
dvx=0.02*vdx; dvz=0.02*vdx;
% [vz,vx]=ndgrid((-8*vtz:dvz:8*vtz)+vdz,0:dvx:8*vtx); % (vpara,vperp) grid
% [vz,vx]=ndgrid((-10*vtz:dvz:10*vtz)+vdz,0:dvx:10*vtx); % (vpara,vperp) grid
[vz,vx]=ndgrid((-2.5*vdx:dvz:2.5*vdx)+vdz,0:dvx:2.5*vdx); % (vpara,vperp) grid
% [vz,vx]=ndgrid((-4*vtz:dvz:4*vtz)+0*vdz,0:dvx:4*vtx); % (vpara,vperp) grid

if(icase==1) % drift bi-Maxwellian ring beam
    bs=vdsr/vtps;
    As=exp(-bs^2)+sqrt(pi)*bs*erfc(-bs);
    coef=1/(sqrt(pi^3)*vtz*vtx^2*As);
    fv=coef*exp(-(vz-vdz).^2/vtz^2-(vx-vdx).^2/vtx^2);
elseif(icase==2) % bi-kappa
    % kappa=4.0;
    coef=1/(sqrt(pi^3*kappa^3)*vtz*vtx^2)*exp(gammaln(kappa+1)-gammaln(kappa-0.5));
    fv=coef*(1+(vz-vdz).^2/(kappa*vtz^2)+vx.^2/(kappa*vtx^2)).^(-kappa-1);
elseif(icase==3) % product bi-kappa
    kappaz=kappa; kappax=kappa;
    coef=1/(sqrt(pi^3*kappaz)*vtz*vtx^2)*exp(gammaln(kappaz+1)-gammaln(kappaz+0.5));
    fv=coef*(1+(vz-vdz).^2/(kappaz*vtz^2)).^(-kappaz-1).*(1+vx.^2/(kappax*vtx^2)).^(-kappax-1);
    % coef=1/(sqrt(pi^3*kappaz)*vtz*vtx^2)*(1-1/kappax)*exp(gammaln(kappaz)-gammaln(kappaz-0.5));
    % fv=coef*(1+(vz-vdz).^2/(kappaz*vtz^2)).^(-kappaz).*(1+vx.^2/(kappax*vtx^2)).^(-kappax);
elseif(icase==6) % kappa-Maxwellian, Bai2024PoP
    % kappa=4.0;
    coef=1/(sqrt(pi^3*kappa)*vtz*vtx^2)*exp(gammaln(kappa)-gammaln(kappa-0.5));
    fv=coef*(1+(vz-vdz).^2/(kappa*vtz^2)).^(-kappa).*exp(-vx.^2/vtx^2);
elseif(icase==4) % shell, 24-12-22 14:51
    bs=vdx/vtx;
    As=2/sqrt(pi)*bs*exp(-bs^2)+(2*bs^2+1)*erfc(-bs);
    coef=1/(sqrt(pi^3)*vtx^3*As);
    fv=coef*exp(-(sqrt(vz.^2+vx.^2)-vdx).^2./vtx^2);
    % fv=exp(-(sqrt(vz.^2+vx.^2)-vdz).^2./vtz^2);
elseif(icase==5) % slowing down, 24-12-22 14:
    coef=3/(4*pi*log(1+vdx^3/vtx^3));
    fv=coef*(sqrt(vz.^2+vx.^2)<=vdx)./(sqrt(vz.^2+vx.^2).^3+vtx^3);
    % fv=1./(sqrt(vz.^2+vx.^2).^3+vtz^3);
    % fv=(sqrt(vz.^2+vx.^2)<=vdz)./(sqrt(vz.^2+vx.^2).^3+vtz^3);
end
sumfv=2*pi*sum(sum(fv.*vx))*dvx*dvz; % should be normalized to 1

fvdat.vx=vx; fvdat.vz=vz; fvdat.fv=fv; fvdat.dvx=dvx; fvdat.dvz=dvz;
fvdat.vtx=vtx; fvdat.vtz=vtz; fvdat.vdx=vdx; fvdat.vdz=vdz;
save(['fvdata',num2str(sid),'.mat'],'fvdat');


%%
close all;
subplot(221);
surf(vz,vx,fv); colorbar;
xlabel('v_z');ylabel('v_x');
subplot(222);
contour(vz/vA,vx/vA,fv,150); colorbar;
xlabel('v_z');ylabel('v_x');

subplot(223);
surf(vz,vx,log(fv)); colorbar;
xlabel('v_z');ylabel('v_x');
subplot(224);
contour(vz/vA,vx/vA,log10(fv),-30:1:-15); colorbar;
xlabel('v_z');ylabel('v_x');
xlim([-4,4]);
ylim([0,4]);
