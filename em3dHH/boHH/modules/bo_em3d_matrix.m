% 18-12-16 08:58, Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, China.
% This file set the BO (PDRK) kernel matrix elements, of electromagnetic
% 3D case.
%
% This new version support: theta=pi/2, kz=0; theta=0, kx=0; and also kz<0.
%
% 18-12-18 21:21 We try two methods to calculate b11, b22, b33 terms, it
% seems that imethod=1 works better, which convergent easily for small N.
% To understand why.
%
% 18-12-28 21:16 A bug is fixed for EM3D-U p13 & p31 terms
% 19-01-29 12:33 Another bug is fixed for EM3D-U p13 terms
% 20-10-25 21:40 Another bug is fixed for EM3D-M p32 vdsy term
% 20-11-14 08:32 Update to can handle the calculation of density and
% velocity perturbation/polarization of each species.
% 21-10-08 09:44 fixed bug for imethod=1 (to check)
% 
% 24-12-19 14:29 remove vdx,vdy,collision,unmagnetized and loss cone
% 24-12-21 07:24 update for arbitrary distribution with Hermite a_{s,lm}
% 24-12-21 11:59 test not ok, sum11tmp should not only for j=1
% 24-12-21 13:36 can run for anayltic alm=1, ok. small difference
% 14:04 fix P32 bug, agree BO2019
% 23:00 for grid distribution, also agree. However, slight difference for
% z-direction expansion
% 24-12-22 fix some bugs, benchmark agree

% % -- main program begin to set the matrix elements --

if(kx==0) % 18-12-16 14:59, to remove the singularity when kx=0. To update.
    kx=1e-30;
end
% if(kz==0) % not needed
%   kz=1e-30;
% end

k=sqrt(kz^2+kx^2);
as=kx*rhocs*sqrt(2); % as=kx*vtps./wcs, note a sqrt(2), 18-11-29 16:34

M=sparse(NN,NN);
snj=0;

if(kz<0) % 18-12-22 00:54
    czjj=-czj;
else
    czjj=czj;
end

czjj_l=zeros(length(czjj),max(lsmax)+5);
for l=0:(max(lsmax)+3) % 24-12-21 08:29
    % to keep l<0 czjj^l=0, i.e., czjj_l(:,1)=czjj^(-1)=0
    czjj_l(:,l+2)=czjj.^l;
end

% initialize
% b11=0; b12=0; b13=0; b21=0; b22=0; b23=0; b31=0; b32=0; b33=0;
% 20-11-14 10:35 update b_ij to b_{ij,s} for calcualte dJ_s
b11s=zeros(1,S); b12s=0.*b11s; b13s=0.*b11s; b21s=0.*b11s; b22s=0.*b11s;
b23s=0.*b11s; b31s=0.*b11s; b32s=0.*b11s; b33s=0.*b11s;
% csnj=zeros(1,3*SNJ);
csnj=zeros(1,SNJ); % 18-12-17 22:42
b11snj=csnj.*0; b12snj=csnj.*0; b13snj=csnj.*0;
b21snj=csnj.*0; b22snj=csnj.*0; b23snj=csnj.*0;
b31snj=csnj.*0; b32snj=csnj.*0; b33snj=csnj.*0;

% 18-12-18 21:17 !!!%%%
imethod=1; % seems =0 not works well or need large N to convergent, % 19-11-30 12:56 especially for ion Bernstein wave
% imethod=0;

for s=1:S % species

    Nsss=-Nss(s):Nss(s); % magnetized species, for different N

    for n=Nsss % sum_N

        for j=1:J % poles of Z(zeta)
            snj=snj+1;

            % magnetized species
            %
            % 18-12-20 17:08 Note the czjj(j) term in csnj, which is to support
            % also kz<0 in Z function, the previous one only correct for kz>0.
            % csnj(snj)=czj(j)*kz*vtzs(s)+kz*vdsz(s)+kx*vdsx(s)-1i*nus(s)+n*wcs(s);
            csnj(snj)=czjj(j)*kz*vtzs(s)+kz*vdsz(s)+n*wcs(s);
            cnj=csnj(snj);

            if(j==1) % 18-12-17 00:08 only need calculate once
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

                if(1==1)
                    % 18-12-02 10:00, calculate An, Bn, Cn using numerical integral
                    Anabb=4.0/As(s)*funAn(n,as(s),ds(s),1,0);
                    Anab0=4.0/As(s)*funAn(n,as(s),ds(s),0,1);
                    Bnabb=4.0/As(s)*funBn(n,as(s),ds(s),1,1);
                    Bnab0=4.0/As(s)*funBn(n,as(s),ds(s),0,2);
                    Cnabb=4.0/As(s)*funCn(n,as(s),ds(s),1,2);
                    Cnab0=4.0/As(s)*funCn(n,as(s),ds(s),0,3);
                end
            end

            sum11tmp1=0; sum11tmp2=0; sum11tmp3=0;
            sum12tmp1=0; sum12tmp2=0;
            sum22tmp1=0; sum22tmp2=0; sum22tmp3=0;
            sum13tmp1=0; sum13tmp2=0;
            sum23tmp1=0; sum23tmp2=0;
            sum32tmp1=0; sum32tmp2=0;
            sum33tmp1=0; sum33tmp2=0; sum33tmp3=0;

            if(1==1)
                for l=0:lsmax(s)
                    for m=0:msmax(s)
                        sum11tmp1=sum11tmp1+aslm{s}(l+1,m+1)*czjj_l(j,l+2)*( ...
                            2*Ans(m+1+2,1)-m*Ans(m-1+2,1));
                        sum11tmp2=sum11tmp2+aslm{s}(l+1,m+1)*( ...
                            2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2))*Ans(m+2,2);

                        sum12tmp1=sum12tmp1+aslm{s}(l+1,m+1)*czjj_l(j,l+2)*( ...
                            2*Bns(m+1+2,1)-m*Bns(m-1+2,1));
                        sum12tmp2=sum12tmp2+aslm{s}(l+1,m+1)*( ...
                            2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2))*Bns(m+2,2);

                        sum22tmp1=sum22tmp1+aslm{s}(l+1,m+1)*czjj_l(j,l+2)*( ...
                            2*Cns(m+1+2,1)-m*Cns(m-1+2,1));
                        sum22tmp2=sum22tmp2+aslm{s}(l+1,m+1)*( ...
                            2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2))*Cns(m+2,2);

                        sum13tmp1=sum13tmp1+aslm{s}(l+1,m+1)*( ...
                            vdsz(s)/vtzs(s)*czjj_l(j,l+2)+czjj_l(j,l+1+2))*( ...
                            2*Ans(m+1+2,1)-m*Ans(m-1+2,1));
                        sum13tmp2=sum13tmp2+aslm{s}(l+1,m+1)*( ...
                            (2*czjj_l(j,l+2+2)-l*czjj_l(j,l+2))+vdsz(s)/vtzs(s)*( ...
                            2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2)))*Ans(m+2,2);

                        % 24-12-21 09:25 to check p31=p13
                        % sum31tmp1=sum31tmp1+aslm{s}(l+1,m+1)*( ...
                        %     vdsz(s)/vtzs(s)*czjj_l(j,l+2)+czjj_l(j,l+1+2))*( ...
                        %     2*Ans(m+1+2,1)-m*Ans(m-1+2,1));
                        % sum31tmp2=sum31tmp2+aslm{s}(l+1,m+1)*( ...
                        %     (2*czjj_l(j,l+2+2)-l*czjj_l(j,l+2))+vdsz(s)/vtzs(s)*( ...
                        %     2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2)))*Ans(m+2,2);

                        sum23tmp1=sum23tmp1+aslm{s}(l+1,m+1)*( ...
                            vdsz(s)/vtzs(s)*czjj_l(j,l+2)+czjj_l(j,l+1+2))*( ...
                            2*Bns(m+1+2,1)-m*Bns(m-1+2,1));
                        sum23tmp2=sum23tmp2+aslm{s}(l+1,m+1)*( ...
                            (2*czjj_l(j,l+2+2)-l*czjj_l(j,l+2))+vdsz(s)/vtzs(s)*( ...
                            2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2)))*Bns(m+2,2);

                        % sum32tmp1=sum32tmp1+aslm{s}(l+1,m+1)*( ...
                        %     vdsz(s)/vtzs(s)*czjj_l(j,l+2)+czjj_l(j,l+1+2))*( ...
                        %     2*Ans(m+1+2,2)-m*Ans(m-1+2,2));
                        % sum32tmp2=sum32tmp2+aslm{s}(l+1,m+1)*( ...
                        %     (2*czjj_l(j,l+2+2)-l*czjj_l(j,l+2))+vdsz(s)/vtzs(s)*( ...
                        %     2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2)))*Ans(m+2,3);
                        sum32tmp1=sum32tmp1+aslm{s}(l+1,m+1)*( ...
                            vdsz(s)/vtzs(s)*czjj_l(j,l+2)+czjj_l(j,l+1+2))*( ...
                            2*Bns(m+1+2,1)-m*Bns(m-1+2,1));
                        sum32tmp2=sum32tmp2+aslm{s}(l+1,m+1)*( ...
                            (2*czjj_l(j,l+2+2)-l*czjj_l(j,l+2))+vdsz(s)/vtzs(s)*( ...
                            2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2)))*Bns(m+2,2);

                        sum33tmp1=sum33tmp1+aslm{s}(l+1,m+1)*( ...
                            (vdsz(s)/vtzs(s))^2*czjj_l(j,l+2)+ ...
                            2*(vdsz(s)/vtzs(s))*czjj_l(j,l+1+2)+czjj_l(j,l+2+2))*( ...
                            2*Ans(m+1+2,1)-m*Ans(m-1+2,1));
                        sum33tmp2=sum33tmp2+aslm{s}(l+1,m+1)*( ...
                            (vdsz(s)/vtzs(s))^2*(2*czjj_l(j,l+1+2)-l*czjj_l(j,l-1+2))+ ...
                            2*(vdsz(s)/vtzs(s))*(2*czjj_l(j,l+2+2)-l*czjj_l(j,l+2))+ ...
                            (2*czjj_l(j,l+3+2)-l*czjj_l(j,l+1+2))...
                            )*Ans(m+2,2);

                        if(j==1)
                            sum11tmp3=sum11tmp3+aslm{s}(l+1,m+1)*funIn(l)*( ...
                                2*Ans(m+1+2,1)-m*Ans(m-1+2,1));
                            sum22tmp3=sum22tmp3+aslm{s}(l+1,m+1)*funIn(l)*( ...
                                2*Cns(m+1+2,1)-m*Cns(m-1+2,1));
                            sum33tmp3=sum33tmp3+aslm{s}(l+1,m+1)*( ...
                                (vdsz(s)/vtzs(s))*(2*funIn(l+1)-l*funIn(l-1))+ ...
                                (2*funIn(l+2)-l*funIn(l)))*Ans(m+2,2);
                        end

                    end
                end
            else
                sum11tmp1=Anabb;
                sum11tmp2=Anab0*czjj(j);
                sum12tmp1=Bnabb;
                sum12tmp2=Bnab0*czjj(j);

                sum22tmp1=Cnabb;
                sum22tmp2=Cnab0*czjj(j);

                sum13tmp1=Anabb*(czjj(j)+vdsz(s)/vtzs(s));
                sum13tmp2=Anab0*(czjj(j)+vdsz(s)/vtzs(s))*czjj(j);

                sum23tmp1=Bnabb*(czjj(j)+vdsz(s)/vtzs(s));
                sum23tmp2=Bnab0*(czjj(j)+vdsz(s)/vtzs(s))*czjj(j);

                sum32tmp1=Bnabb*(czjj(j)+vdsz(s)/vtzs(s));
                sum32tmp2=Bnab0*(czjj(j)+vdsz(s)/vtzs(s))*czjj(j);

                sum33tmp1=Anabb*(czjj(j)+vdsz(s)/vtzs(s))^2;
                sum33tmp2=Anab0*(czjj(j)+vdsz(s)/vtzs(s))^2*czjj(j);

                sum11tmp3=Anabb;
                sum22tmp3=Cnabb;
                sum33tmp3=Anab0/2;
            end

            % 18-12-18 17:05
            if(imethod==1 && j==1)
                % p11snj_b=(n*wcs(s)/kx)*n*wcs(s)/kx*Anabb/vtps(s)^2;
                % p22snj_b=Cnabb;
                % p33snj_b=Anab0/2; % 18-12-18 23:55

                % 24-12-21 11:22
                p11snj_b=(n*wcs(s)/(kx*vtps(s)))^2*sum11tmp3;
                p22snj_b=sum22tmp3;
                p33snj_b=sum33tmp3;
                b11s(s)=b11s(s)-wps2(s)*p11snj_b; % 20-11-14 10:44 update
                b22s(s)=b22s(s)-wps2(s)*p22snj_b;
                b33s(s)=b33s(s)-wps2(s)*p33snj_b;
            end

            tmp=wps2(s)*bzj(j)/cnj;

            % p11snj=(n*wcs(s)/kx)*(n*wcs(s)/kx*(n*wcs(s)...
            %   )*Anabb/vtps(s)^2 + Anab0*(n*wcs(s)/kx...
            %   )*kz*czjj(j)/vtzs(s));
            p11snj=(n*wcs(s)/(kx*vtps(s)))^2*(n*wcs(s)*sum11tmp1+ ...
                (kz*vtzs(s))*(vtps(s)/vtzs(s))^2*sum11tmp2);

            % p12snj=(n*wcs(s)/kx)*(1i*((n*wcs(s))*Bnabb/...
            %   vtps(s)+Bnab0*vtps(s)/vtzs(s)*kz*czjj(j)));
            p12snj=1i*(n*wcs(s)/(kx*vtps(s)))*(n*wcs(s)*sum12tmp1+ ...
                (kz*vtzs(s))*(vtps(s)/vtzs(s))^2*sum12tmp2);

            % p21snj=n*wcs(s)/kx*(n*wcs(s))*(-1i* ...
            %   vtps(s)*Bnabb)/vtps(s)^2 + ...
            %   (n*wcs(s)/kx)/vtzs(s)*(-1i*vtps(s)*Bnab0)*kz*czjj(j);
            p21snj=-p12snj;

            % p22snj=n*wcs(s)*Cnabb + vtps(s)/vtzs(s)*(...
            %   vtps(s)*Cnab0)*kz*czjj(j);
            p22snj=n*wcs(s)*sum22tmp1+ ...
                (kz*vtzs(s))*(vtps(s)/vtzs(s))^2*sum22tmp2;

            % p13snj=(n*wcs(s)/kx)*(n*wcs(s)*(czjj(j)*vtzs(s)+...
            %   vdsz(s))*Anabb/vtps(s)^2 + Anab0*(kz*czjj(j)^2+(...
            %   kz*vdsz(s))/vtzs(s)*czjj(j)));
            p13snj=(n*wcs(s)/(kx*vtps(s)))*((vtzs(s)/vtps(s))*n*wcs(s)*sum13tmp1+ ...
                (kz*vtps(s))*sum13tmp2);

            % p31snj=n*wcs(s)/kx*(n*wcs(s))*(czjj(j)*vtzs(s)+...
            %   vdsz(s))*Anabb/vtps(s)^2 +Anab0*(n*wcs(s)/kx)*(czjj(j)+ ...
            %   vdsz(s)/vtzs(s))*kz*czjj(j);
            p31snj=p13snj;

            % p23snj=(-1i*vtps(s)*Bnabb)*(czjj(j)*vtzs(s)+vdsz(s) ...
            %   )*n*wcs(s)/vtps(s)^2 + (-1i*vtps(s)*Bnab0 ...
            %   )*(kz*czjj(j)+(kz*vdsz(s))/vtzs(s))*czjj(j);
            p23snj=-1i*((vtzs(s)/vtps(s))*n*wcs(s)*sum23tmp1+ ...
                (kz*vtps(s))*sum23tmp2);

            % fixed a bug, 20-10-25 21:15
            % p32snj= 1i*((n*wcs(s))* ...
            %   Bnabb*vtzs(s)/vtps(s)+Bnab0*kz*vtps(s)*czjj(j))*(czjj(j)+...
            %   vdsz(s)/vtzs(s));
            p32snj=1i*((vtzs(s)/vtps(s))*n*wcs(s)*sum32tmp1+ ...
                (kz*vtps(s))*sum32tmp2);

            % p33snj=( n*wcs(s)*(vtzs(s)^2*czjj(j)+vtzs(s)*vdsz(s))*Anabb/ ...
            %   vtps(s)^2 + Anab0*czjj(j)*(czjj(j)*kz*vtzs(s)+(kz*vdsz(s)...
            %   )) )*(czjj(j)+ vdsz(s)/vtzs(s));
            p33snj=(vtzs(s)/vtps(s))*((vtzs(s)/vtps(s))*n*wcs(s)*sum33tmp1+ ...
                (kz*vtps(s))*sum33tmp2);

            b11snj(snj)=b11snj(snj)+tmp*p11snj;
            %     b11=b11-tmp*p11snj; % the extra wps2(s) term will be add in final step
            b11s(s)=b11s(s)-tmp*p11snj; %  20-11-14  10:45 update

            b12snj(snj)=b12snj(snj)+tmp*p12snj;
            %     b12=b12-tmp*p12snj;
            b12s(s)=b12s(s)-tmp*p12snj; %  20-11-14  10:45 update

            b21snj(snj)=b21snj(snj)+tmp*p21snj;
            b21s(s)=b21s(s)-tmp*p21snj; %

            b22snj(snj)=b22snj(snj)+tmp*p22snj;
            b22s(s)=b22s(s)-tmp*p22snj; % the extra wps2(s) term will be add in final step

            b13snj(snj)=b13snj(snj)+tmp*p13snj;
            b13s(s)=b13s(s)-tmp*p13snj; %

            b31snj(snj)=b31snj(snj)+tmp*p31snj;
            b31s(s)=b31s(s)-tmp*p31snj; %

            b23snj(snj)=b23snj(snj)+tmp*p23snj;
            b23s(s)=b23s(s)-tmp*p23snj; %

            b32snj(snj)=b32snj(snj)+tmp*p32snj;
            b32s(s)=b32s(s)-tmp*p32snj; %

            b33snj(snj)=b33snj(snj)+tmp*p33snj;
            b33s(s)=b33s(s)-tmp*p33snj; % the extra wps2(s) term will be add in final step

        end
    end

    % extra term for diag elements, 20-11-14 10:49 update
    if(imethod==0) % seems need large N to convergent
        b11s(s)=b11s(s)-wps2(s);
        b22s(s)=b22s(s)-wps2(s);
        b33s(s)=b33s(s)-wps2(s);
        %  else % for unmagnetized species
        %      b11s(s)=b11s(s)-wps2(intersect(s,jds)); % to check
        %      b22s(s)=b22s(s)-wps2(intersect(s,jds));
        %      b33s(s)=b33s(s)-wps2(intersect(s,jds));
        %  elseif(exist('intersect','var')) % for unmagnetized species, 21-09-18 16:07 fixed bug for imethod=1
    elseif(imethod==1 && exist('intersect','var')) % for unmagnetized species, 21-10-08 09:46 fixed bug for imethod=1
        b11s(s)=b11s(s)-wps2(intersect(s,jds)); % to check
        b22s(s)=b22s(s)-wps2(intersect(s,jds));
        b33s(s)=b33s(s)-wps2(intersect(s,jds));
    end

end

for snj=1:SNJ % set the eigen matrix
    jjx=snj+0*SNJ1;
    jjy=snj+1*SNJ1;
    jjz=snj+2*SNJ1;
    % v_snjx
    M=M+sparse(jjx,jjx,csnj(snj),NN,NN)+...
        sparse(jjx,SNJ3+1,b11snj(snj),NN,NN)+...
        sparse(jjx,SNJ3+2,b12snj(snj),NN,NN)+...
        sparse(jjx,SNJ3+3,b13snj(snj),NN,NN);

    % v_snjy
    M=M+sparse(jjy,jjy,csnj(snj),NN,NN)+...
        sparse(jjy,SNJ3+1,b21snj(snj),NN,NN)+...
        sparse(jjy,SNJ3+2,b22snj(snj),NN,NN)+...
        sparse(jjy,SNJ3+3,b23snj(snj),NN,NN);

    % v_snjz
    M=M+sparse(jjz,jjz,csnj(snj),NN,NN)+...
        sparse(jjz,SNJ3+1,b31snj(snj),NN,NN)+...
        sparse(jjz,SNJ3+2,b32snj(snj),NN,NN)+...
        sparse(jjz,SNJ3+3,b33snj(snj),NN,NN);

end

% E(J), J_{x,y,z}=j_{x,y,z}+sum(v_snj{x,y,z})
tp=-1;
jj=(0*SNJ1+1):(1*SNJ1); ii=0.*jj+SNJ3+1; M=M+sparse(ii,jj,tp,NN,NN);
jj=(1*SNJ1+1):(2*SNJ1); ii=0.*jj+SNJ3+2; M=M+sparse(ii,jj,tp,NN,NN);
jj=(2*SNJ1+1):(3*SNJ1); ii=0.*jj+SNJ3+3; M=M+sparse(ii,jj,tp,NN,NN);

% % jx(E), jy(E), jz(E)
% M=M+sparse(1*SNJ1,SNJ3+1,b11,NN,NN)+...
%   sparse(1*SNJ1,SNJ3+2,b12,NN,NN)+...
%   sparse(1*SNJ1,SNJ3+3,b13,NN,NN)+...
%   sparse(2*SNJ1,SNJ3+1,b21,NN,NN)+...
%   sparse(2*SNJ1,SNJ3+2,b22,NN,NN)+...
%   sparse(2*SNJ1,SNJ3+3,b23,NN,NN)+...
%   sparse(3*SNJ1,SNJ3+1,b31,NN,NN)+...
%   sparse(3*SNJ1,SNJ3+2,b32,NN,NN)+...
%   sparse(3*SNJ1,SNJ3+3,b33,NN,NN);

% jsx(E), jsy(E), jsz(E), 20-11-14 11:02 separate each species
for s=1:S
    M=M+sparse(1*SNJ1-S+s,SNJ3+1,b11s(s),NN,NN)+...
        sparse(1*SNJ1-S+s,SNJ3+2,b12s(s),NN,NN)+...
        sparse(1*SNJ1-S+s,SNJ3+3,b13s(s),NN,NN)+...
        sparse(2*SNJ1-S+s,SNJ3+1,b21s(s),NN,NN)+...
        sparse(2*SNJ1-S+s,SNJ3+2,b22s(s),NN,NN)+...
        sparse(2*SNJ1-S+s,SNJ3+3,b23s(s),NN,NN)+...
        sparse(3*SNJ1-S+s,SNJ3+1,b31s(s),NN,NN)+...
        sparse(3*SNJ1-S+s,SNJ3+2,b32s(s),NN,NN)+...
        sparse(3*SNJ1-S+s,SNJ3+3,b33s(s),NN,NN);
end

% E(B)
M=M+sparse(SNJ3+1,SNJ3+5,c2*kz,NN,NN)+...
    sparse(SNJ3+2,SNJ3+4,-c2*kz,NN,NN)+...
    sparse(SNJ3+2,SNJ3+6,c2*kx,NN,NN)+...
    sparse(SNJ3+3,SNJ3+5,-c2*kx,NN,NN);

% B(E)
M=M+sparse(SNJ3+4,SNJ3+2,-kz,NN,NN)+...
    sparse(SNJ3+5,SNJ3+1,kz,NN,NN)+...
    sparse(SNJ3+5,SNJ3+3,-kx,NN,NN)+...
    sparse(SNJ3+6,SNJ3+2,kx,NN,NN);

% for Darwin model, use lambda*MB*X=M*X, 18-12-17 22:55
if(iem==2) % to check
    MB=speye(NN);
    % modify the omega*E equation
    MB=MB+sparse(SNJ3+1,SNJ3+1,kx^2/k^2,NN,NN)+...
        sparse(SNJ3+3,SNJ3+1,kx*kz/k^2,NN,NN)+...
        sparse(SNJ3+1,SNJ3+3,kx*kx/k^2,NN,NN)+...
        sparse(SNJ3+3,SNJ3+3,kz^2/k^2,NN,NN)-...
        sparse(SNJ3+1,SNJ3+1,1,NN,NN)-...
        sparse(SNJ3+2,SNJ3+2,1,NN,NN)-...
        sparse(SNJ3+3,SNJ3+3,1,NN,NN);

end

% % -- main program end of set the matrix elements --
