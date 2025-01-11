% Hua-sheng XIE, huashengxie@gmail.com, 2024-12-29 01:14
% Solve the kinetic dispersion relation use fsolve()
% for arbitrary distribution version with Hermite-Hermite expansion
% To update: 1.fsolve work not well for damped mode; 2. integral write out
% from fDrHH() to speed up

close all; clear; clc;
global S c2 wcs wps2 rhocs kx kz vtzs vtps vdsz ds As Nss aslm msmax lsmax;

% % Default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg

N=2;
setNs=[];
B0=0.1; % magnetic field (Tesla)
munit=mp; % set the default normalization for mass

% % read input parameters
pardat = importdata('bo.in', ' ', 1);
[S, col]=size(pardat.data);
if(col~=8)
    disp('Wrong input data !!!');
end

for s=1:S
    qs0(s)=pardat.data(s,1); % charge, q/e
    ms0(s)=pardat.data(s,2); % mass, m/mp
    ns0(s)=pardat.data(s,3); % desity, m^-3
    Tzs0(s)=pardat.data(s,4); % parallel temperature, eV
    Tps0(s)=pardat.data(s,5); % perp temperature, eV
    vdsz0(s)=pardat.data(s,6); % para drift velocity, vds_z/c
    vdsr0(s)=pardat.data(s,7); % perp ring beam drift velocity, vds_r/c
    ifv0(s)=pardat.data(s,8); % fv_type, 1=analy drift bi-Maxwellian ring beam, 0=grid
end

Qtotal=sum(qs0.*ns0);
Jztotal=sum(qs0.*ns0.*vdsz0);

if((Qtotal~=0) || (Jztotal~=0) )
    disp('Warning: Total charge or current not zero !!!');
    %input(['Warning: Total charge or current not zero!!',...
    %    'Press any key to continue.']);
end

qs=qs0*qe; % * electron charge, e -> C (coulomb)
% ms=ms0*mp; % * proton mass, m_p -> kg, 18-10-13 09:57
ms=ms0*munit; % * unit mass, -> kg, 18-12-23 22:59

% 24-12-19 19:42
lsmax=zeros(S,1); msmax=zeros(S,1);
for s=1:S
    if(ifv0(s)==1)
        lsmax(s)=0; msmax(s)=0;
        aslm{s}=1.0;
        % coef=1/(sqrt(pi^3)*vtzs(s)*vtps(s)^2*As(s));
        % aslm=coef;
    else
        % including lsmax, msmax, aslm(lsmax+1,msmax+1), fv_s(vpara,vperp)
        % and Tzs, Tps, vdsz, vdsr should be the same as in 'bo.in'
        % load(['../input/fvdat',num2str(s),'.mat']);

        % 24-12-21 15:21
        load(['../241221/genfv/fvceff',num2str(s),'.mat']);
        aslm{s}=fvc.alm;
        vdsz(s)=fvc.vdz;
        vdsr(s)=fvc.vdr;
        vtzs(s)=fvc.vtz;
        vtps(s)=fvc.vtp;
        vdsz0(s)=vdsz(s)/sqrt(c2);
        vdsr0(s)=vdsr(s)/sqrt(c2);
        Tzs0(s)=ms(s)*vtzs(s)^2/(2*qe);
        Tps0(s)=ms(s)*vtps(s)^2/(2*qe);
        lsmax(s)=size(aslm{s},1)-1; % 24-12-21 21:54
        % lsmax(s)=min(size(aslm{s},1)-1,18); % 24-12-29 15:47
        msmax(s)=size(aslm{s},2)-1;
    end
end

Tzs=Tzs0*qe/kB; % T//, eV -> K (eV -> J * J -> K)
Tps=Tps0*qe/kB; % Tpr, eV -> K
vdsz=vdsz0*sqrt(c2); % vds_z, speed of light c -> m/s
vdsr=vdsr0*sqrt(c2); % vds_r, speed of light c -> m/s

Ns=(2*N+1)+0.*(1:S); % # of sum_N harmonic for magnetized species
if(~isempty(setNs)) % 18-12-15 23:51
    for js=1:size(setNs,1)
        if(setNs(js,1)<=S)
            %       Ns(js)=2*setNs(js,2)+1; % change the N for some species
            Ns(setNs(js,1))=2*setNs(js,2)+1; % 20-11-14 09:14
        end
    end
end
Nss=floor((Ns-1)/2); % 18-11-21 12:18

vtzs=sqrt(2*kB*Tzs./ms); % para thermal velocity, note the sqrt(2)
vtps=sqrt(2*kB*Tps./ms); % perp thermal velocity, note the sqrt(2)
lambdaDs=sqrt(epsilon0*kB*Tzs./(ns0.*qs.^2)); % Debye length, Tzs
kDs=1./lambdaDs;
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wcs=B0*qs./ms; % cyclotron frequency
% rhocs=sqrt(kB*Tps./ms)./abs(wcs); % cyclotron radius
rhocs=sqrt(kB*Tps./ms)./wcs; % cyclotron radius, 2018-06-13 21:47

wps2=wps.^2;
lmdT=Tzs./Tps;

ds=vdsr./vtps;
% normalized coefficient for F(v_perp)
As=exp(-ds.^2)+sqrt(pi)*ds.*erfc(-ds);

betasz=2*mu0*kB.*ns0.*Tzs./B0^2; % beta_para
betasp=2*mu0*kB.*ns0.*Tps./B0^2; % beta_perp
vA=B0/sqrt(mu0*sum(ms.*ns0)); % Alfven speed

wn=abs(wcs(1)); % normalized omega by wn, modify here to use other wn
kn=wn/vA;

%%
jplt=1;
if(jplt==1) % set kperp range
    kk=0.1:0.025:0.4;
elseif(jplt==2)
    kk=10.^(-1:0.1:4);
end

theta=45*pi/180;
% theta=0*pi/180;

nk=length(kk);
wwx=zeros(nk,5);
runtime=cputime;
for jk=1:nk

    kx=kk(jk)*kn*sin(theta); % set kperp
    kz=kk(jk)*kn*cos(theta); % set kperp

    jk
    if(mod(jk,1)==0)
        if(jk==1)
            x01=(0.0+0.01i)*wn;
        else
            % x01=wx1;
        end
        %     warning off;
        options = optimset('Display','on','TolFun',1e-14,'Tolx',1e-14);
        %     options = optimset('Display','iter','TolFun',1e-100);
        wx1=fsolve(@fDrHH,x01,options);
        tmp=wx1/wn
        wwx(jk,1)=wx1/wn;
    else
        wwx(jk,1)=NaN+1i*NaN;
    end

end
runtime=cputime-runtime;
%%
close all;
h=figure('unit','normalized','Position',[0.01 0.2 0.7 0.5],...
    'DefaultAxesFontSize',15);

if(jplt==1)

    plot(kk,real(wwx(:,1)),'rx','linewidth',1); hold on;
    % ylim([-0.1,0.5]);
    xlim([0,max(kk)]);
    subplot(122);
    plot(kk,imag(wwx(:,1)),'rx','linewidth',1); hold on;
    xlim([0,max(kk)]);
    ylim([0,0.12]);

elseif(jplt==2)

    semilogx(kk,real(wwx(:,1)),'rx','linewidth',1); hold on;
    xlim([min(kk),max(kk)]);
    % ylim([0,10]);

    %     subplot(121);semilogx(kxx,kxx*kn*c/wn,'k:','linewidth',2); hold on;

end
subplot(121);
% title(['S=',num2str(S),', N=',num2str(N),', B_0=',num2str(B0,3),...
%     ', \theta=',num2str(theta*180/pi),'\circ, ',10,...
%     ', m_s=',num2str(ms.'/me,4),', q_s=',num2str(qs.'/qe,3)]);
xlabel(['k/k_n, runtime=',num2str(runtime),'s']);
ylabel('\omega_r/\omega_{n}');
subplot(122);
% title(['n_{s0}=',num2str(ns0.',3),', T_{sz}=',num2str(Tsz.'/(qe/kB),3),10,...
%     ', T_{sp}=',num2str(Tsp.'/(qe/kB),3),', nk=',num2str(nk)]);
xlabel(['k/k_n, k_n=',num2str(kn,3),', \rho_{cs}=',num2str(rhocs,3)]);
ylabel(['\omega_i/\omega_{n}, \omega_n=',num2str(wn,3)]);

% save figure
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);

print(gcf,'-dpng',['em3d_S=',num2str(S),',N=',num2str(N),',B0=',num2str(B0,3),...
    ',ms=',num2str(ms/me,4),',qs=',num2str(qs/qe,3),...
    ',ns0=',num2str(ns0,3),',theta=',num2str(theta*180/pi,3),',kn=',num2str(kn,3),...
    ',wn=',num2str(wn,3),',nk=',num2str(nk),'.png']);
