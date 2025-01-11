% Hua-sheng XIE, huashengxie@gmail.com, ENN, 2025-01-02 21:53
% Hot magnetized plasma dispersion relation, for arbitrary distribution
% with GPDF expansion for both para & perp directions
% To update: 1.fsolve work not well for damped mode; 2. integral write out
% from fDrGG() to speed up

function f=fDrGG(w) % DR

global S c2 wcs wps2 rhocs kx kz vtzs vtps vdsz ds As Nss aslm msmax lsmax ifv0;

% vts=sqrt(2*kB*Ts./ms); % rhocs=vts/(sqrt(2)*wcs);

if(kx==0) % 18-12-16 14:59, to remove the singularity when kx=0. To update.
    kx=1e-30;
end

if(kz==0)
    kz=0.0001*kx;
end
% k=sqrt(kz^2+kx^2);
as=kx*rhocs*sqrt(2); % as=kx*vtps./wcs, note a sqrt(2), 18-11-29 16:34
nx=sqrt(c2)*kx/w; nz=sqrt(c2)*kz/w;
% bs=0.5*(kx*rhocs).^2;

Kxx=1; Kxy=0; Kxz=0;
Kyx=0; Kyy=1; Kyz=0;
Kzx=0; Kzy=0; Kzz=1;

for s=1:S % species

    Nsss=-Nss(s):Nss(s); % magnetized species, for different N

    for n=Nsss % sum_N

        Xsn_xx=0; Xsn_xy=0; Xsn_xz=0; Xsn_yy=0; Xsn_yz=0; Xsn_zz=0;
        if(ifv0(s)==1) % Hermite-Hermite lowest order, drift bi-Maxwellian ring beam
            Ans=zeros(msmax(s)+4,2); Bns=0.*Ans; Cns=0.*Ans;
            for m=0:(msmax(s)+2)
                % 24-12-21 07:52, calculate An, Bn, Cn using numerical integral
                idm=m+2; % to also store m-1 with m=0
                Ans(idm,1)=2.0/As(s)*funAn(n,as(s),ds(s),m,0); % p=0
                Ans(idm,2)=2.0/As(s)*funAn(n,as(s),ds(s),m,1); % p=1
                Bns(idm,1)=2.0/As(s)*funBn(n,as(s),ds(s),m,1); % p=1
                Bns(idm,2)=2.0/As(s)*funBn(n,as(s),ds(s),m,2); % p=2
                Cns(idm,1)=2.0/As(s)*funCn(n,as(s),ds(s),m,2); % p=2
                Cns(idm,2)=2.0/As(s)*funCn(n,as(s),ds(s),m,3); % p=3
            end

            zetasn=(w-kz*vdsz(s)-n*wcs(s))/(kz*vtzs(s));
            Zlns=zeros(lsmax(s)+5,1);
            Il=zeros(lsmax(s)+5,1);
            for l=0:(lsmax(s)+3)
                Zlns(l+2)=funZl(zetasn,l);
                Il(l+2)=funIn(l);
            end

            for l=0:lsmax(s)
                for m=0:msmax(s)
                    Xsn_xx=Xsn_xx+aslm{s}(l+1,m+1)*((n*wcs(s)/(kz*vtzs(s))*Zlns(l+2)- ...
                        Il(l+2))*(2*Ans(m+1+2,1)-m*Ans(m-1+2,1))+ ...
                        Ans(m+2,2)*(vtps(s)/vtzs(s))^2*(2*Zlns(l+1+2)-l*Zlns(l-1+2)));
                    Xsn_xy=Xsn_xy+aslm{s}(l+1,m+1)*((n*wcs(s)/(kz*vtzs(s))*Zlns(l+2)- ...
                        Il(l+2))*(2*Bns(m+1+2,1)-m*Bns(m-1+2,1))+ ...
                        Bns(m+2,2)*(vtps(s)/vtzs(s))^2*(2*Zlns(l+1+2)-l*Zlns(l-1+2)));
                    Xsn_yy=Xsn_yy+aslm{s}(l+1,m+1)*((n*wcs(s)/(kz*vtzs(s))*Zlns(l+2)- ...
                        Il(l+2))*(2*Cns(m+1+2,1)-m*Cns(m-1+2,1))+ ...
                        Cns(m+2,2)*(vtps(s)/vtzs(s))^2*(2*Zlns(l+1+2)-l*Zlns(l-1+2)));
                    Xsn_xz=Xsn_xz+aslm{s}(l+1,m+1)*(n*wcs(s)/(kz*vtps(s) ...
                        )*(vdsz(s)/vtzs(s)*Zlns(l+2)+ ...
                        Zlns(l+1+2))*(2*Ans(m+1+2,1)-m*Ans(m-1+2,1))+ ...
                        Ans(m+2,2)*(vtps(s)/vtzs(s))*((2*Zlns(l+2+2)-l*Zlns(l+2))+ ...
                        vdsz(s)/vtzs(s)*(2*Zlns(l+1+2)-l*Zlns(l-1+2))));
                    Xsn_yz=Xsn_yz+aslm{s}(l+1,m+1)*((n*wcs(s)/(kz*vtps(s) ...
                        )*(vdsz(s)/vtzs(s)*Zlns(l+2)+ ...
                        Zlns(l+1+2)))*(2*Bns(m+1+2,1)-m*Bns(m-1+2,1))+ ...
                        Bns(m+2,2)*(vtps(s)/vtzs(s))*((2*Zlns(l+2+2)-l*Zlns(l+2))+ ...
                        vdsz(s)/vtzs(s)*(2*Zlns(l+1+2)-l*Zlns(l-1+2))));
                    Xsn_zz=Xsn_zz+aslm{s}(l+1,m+1)*((n*wcs(s)/(kz*vtzs(s) ...
                        )*((vdsz(s)/vtzs(s))^2*Zlns(l+2)+ ...
                        2*(vdsz(s)/vtzs(s))*Zlns(l+1+2)+Zlns(l+2+2)))*( ...
                        2*Ans(m+1+2,1)-m*Ans(m-1+2,1))+ ...
                        Ans(m+2,2)*(vtps(s)/vtzs(s))^2*((2*vdsz(s)/vtzs(s)*( ...
                        Zlns(l+2+2)-Il(l+1+2))+ ...
                        2*(Zlns(l+3+2)-Il(l+2+2))-l*vdsz(s)/vtzs(s)*( ...
                        Zlns(l+2)-Il(l-1+2))-l*(Zlns(l+1+2)-Il(l+2)))+ ...
                        vdsz(s)/vtzs(s)*(vdsz(s)/vtzs(s)*(2*Zlns(l+1+2)- ...
                        l*Zlns(l-1+2))+2*(Zlns(l+2+2)-l*Zlns(l+2)) )));
                end
            end

        else % 25-01-02 21:55 GPDF-GPDF

            Ans=zeros(2*msmax(s)+6,2); Bns=0.*Ans; Cns=0.*Ans;
            mmin=-msmax(s);
            for m=-(msmax(s)+1):(msmax(s)+1)
                idm=m+2-mmin; % to also store m-1 with m=-mmax
                Ans(idm,1)=funGnmp(as(s),n,m,0,1);%funAn(n,as(s),ds(s),m,0); % p=0
                Ans(idm,2)=funGnmp(as(s),n,m,1,1); % p=1
                Bns(idm,1)=funGnmp(as(s),n,m,1,2); % p=1
                Bns(idm,2)=funGnmp(as(s),n,m,2,2); % p=2
                Cns(idm,1)=funGnmp(as(s),n,m,2,3); % p=2
                Cns(idm,2)=funGnmp(as(s),n,m,3,3); % p=3
            end

            zetasn=(w-n*wcs(s))/(kz*vtzs(s));
            Zlpns=zeros(2*lsmax(s)+5,3);
            lmin=-lsmax(s);
            for l=-(lsmax(s)+1):(lsmax(s)+1)
                Zlpns(l-lmin+2,1)=funZlp(zetasn,l,0);
                Zlpns(l-lmin+2,2)=funZlp(zetasn,l,1);
                Zlpns(l-lmin+2,3)=funZlp(zetasn,l,2);
            end

            for l=-lsmax(s):lsmax(s)
                for m=-msmax(s):msmax(s)
                    Xsn_xx=Xsn_xx+aslm{s}(l-lmin+1,m-mmin+1)*1i*((n*wcs(s)/(kz*vtzs(s))*Zlpns(l-lmin+2,1)+ ...
                        (l==0))*((m-2)*Ans(m-1-mmin+2,1)+2*m*Ans(m-mmin+2,1)+(m+2)*Ans(m+1-mmin+2,1))+ ...
                        Ans(m-mmin+2,2)*(vtps(s)/vtzs(s))^2*((l-1)*Zlpns(l-1-lmin+2,1)+ ...
                        2*l*Zlpns(l-lmin+2,1)+(l+1)*Zlpns(l+1-lmin+2,1)));
                    Xsn_xy=Xsn_xy+aslm{s}(l-lmin+1,m-mmin+1)*1i*((n*wcs(s)/(kz*vtzs(s))*Zlpns(l-lmin+2,1) ...
                        )*((m-2)*Bns(m-1-mmin+2,1)+2*m*Bns(m-mmin+2,1)+(m+2)*Bns(m+1-mmin+2,1))+ ...
                        Bns(m-mmin+2,2)*(vtps(s)/vtzs(s))^2*((l-1)*Zlpns(l-1-lmin+2,1)+ ...
                        2*l*Zlpns(l-lmin+2,1)+(l+1)*Zlpns(l+1-lmin+2,1)));
                    Xsn_yy=Xsn_yy+aslm{s}(l-lmin+1,m-mmin+1)*1i*((n*wcs(s)/(kz*vtzs(s))*Zlpns(l-lmin+2,1)+ ...
                        (l==0))*((m-2)*Cns(m-1-mmin+2,1)+2*m*Cns(m-mmin+2,1)+(m+2)*Cns(m+1-mmin+2,1))+ ...
                        Cns(m-mmin+2,2)*(vtps(s)/vtzs(s))^2*((l-1)*Zlpns(l-1-lmin+2,1)+ ...
                        2*l*Zlpns(l-lmin+2,1)+(l+1)*Zlpns(l+1-lmin+2,1)));
                    Xsn_xz=Xsn_xz+aslm{s}(l-lmin+1,m-mmin+1)*1i*((n*wcs(s)/(kz*vtps(s))*Zlpns(l-lmin+2,2) ...
                        )*((m-2)*Ans(m-1-mmin+2,1)+2*m*Ans(m-mmin+2,1)+(m+2)*Ans(m+1-mmin+2,1))+ ...
                        Ans(m-mmin+2,2)*(vtps(s)/vtzs(s))*((l-1)*Zlpns(l-1-lmin+2,2)+ ...
                        2*l*Zlpns(l-lmin+2,2)+(l+1)*Zlpns(l+1-lmin+2,2)));
                    Xsn_yz=Xsn_yz+aslm{s}(l-lmin+1,m-mmin+1)*1i*((n*wcs(s)/(kz*vtps(s))*Zlpns(l-lmin+2,2) ...
                        )*((m-2)*Bns(m-1-mmin+2,1)+2*m*Bns(m-mmin+2,1)+(m+2)*Bns(m+1-mmin+2,1))+ ...
                        Bns(m-mmin+2,2)*(vtps(s)/vtzs(s))*((l-1)*Zlpns(l-1-lmin+2,2)+ ...
                        2*l*Zlpns(l-lmin+2,2)+(l+1)*Zlpns(l+1-lmin+2,2)));
                    Xsn_zz=Xsn_zz+aslm{s}(l-lmin+1,m-mmin+1)*1i*((n*wcs(s)/(kz*vtzs(s))*Zlpns(l-lmin+2,3) ...
                        )*((m-2)*Ans(m-1-mmin+2,1)+2*m*Ans(m-mmin+2,1)+(m+2)*Ans(m+1-mmin+2,1))+ ...
                        Ans(m-mmin+2,2)*(vtps(s)/vtzs(s))^2*((l-1)*Zlpns(l-1-lmin+2,3)+ ...
                        2*l*Zlpns(l-lmin+2,3)+(l+1)*Zlpns(l+1-lmin+2,3)+2i*(l==0)));
                end
            end
        end

        Xsn_xx=(n*wcs(s)/(kx*vtps(s)))^2*Xsn_xx;
        Xsn_xy=1i*(n*wcs(s)/(kx*vtps(s)))*Xsn_xy;
        % Xsn_xz=(vtzs(s)/vtps(s))*(n*wcs(s)/(kx*vtps(s)))*Xsn_xz;
        Xsn_xz=(n*wcs(s)/(kx*vtps(s)))*Xsn_xz;
        % Xsn_yy=Xsn_yy;
        Xsn_yz=-1i*Xsn_yz; % to check '+/-'
        Xsn_zz=(vtzs(s)/vtps(s))^2*Xsn_zz;

        Xsn_yx=-Xsn_xy;
        Xsn_zx=Xsn_xz;
        Xsn_zy=-Xsn_yz;

        term1=wps2(s)/w^2;
        Kxx=Kxx+term1*Xsn_xx;
        Kxy=Kxy+term1*Xsn_xy;
        Kxz=Kxz+term1*Xsn_xz;
        Kyx=Kyx+term1*Xsn_yx;
        Kyy=Kyy+term1*Xsn_yy;
        Kyz=Kyz+term1*Xsn_yz;
        Kzx=Kzx+term1*Xsn_zx;
        Kzy=Kzy+term1*Xsn_zy;
        Kzz=Kzz+term1*Xsn_zz;
    end
end

Dxx=Kxx-nz^2; Dyx=Kyx; Dxy=Kxy;
Dxz=Kxz+nx*nz; Dyy=Kyy-(nx^2+nz^2); Dzx=Kzx+nx*nz;
Dzy=Kzy; Dyz=Kyz; Dzz=Kzz-nx^2;

DD=Dxx*Dyy*Dzz+Dyx*Dzy*Dxz+Dzx*Dyz*Dxy-...
    Dxz*Dyy*Dzx-Dyz*Dzy*Dxx-Dzz*Dyx*Dxy;
f=DD;
% f=DD^(1/3);
% tmp=max([abs(Dxx),abs(Dyy),abs(Dzz),abs(Dxy),abs(Dxz),abs(Dyz)]);
% tmp=mean([abs(Dxx*Dyy*Dzz),abs(Dyx*Dzy*Dxz),abs(Dzx*Dyz*Dxy)]);
% f=DD/(nx^2+nz^2); % 24-12-31 15:31
% f=DD/abs(Dzz);

end
