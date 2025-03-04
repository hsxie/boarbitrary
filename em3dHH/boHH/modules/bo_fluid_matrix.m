% 20-12-22 22:34, Hua-sheng XIE, huashengxie@gmail.com, ENN-FTRC, China.
% This file set the BO (PDRK & PDRF) kernel matrix elements, of
% electromagnetic 3D fluid case.
% This is rewritten/update version of [Xie2014] PDRF:
% [Xie2014] H. S. Xie, PDRF: A general dispersion relation solver for 
% magnetized multi-fluid plasma, Comput. Phys. Comm. 185 (2014) 670-675.
% Note: There are several differences from 2014 version.
% 24-12-19 14:08 remove vdsx&vdsy

vdsx=0.*vdsz; vdsy=0.*vdsz; % 24-12-19 14:08

% % -- main program begin to set the matrix elements --

k=sqrt(kz^2+kx^2);
kvds=kx*vdsx+kz*vdsz;
if (B0~=0)
    PspB0=Psp/B0;
    PszB0=Psz/B0;
else
    PspB0=0.0.*Psp;
    PszB0=0.0.*Psz;
end

% eigen matrix
M=sparse(NN,NN);
for s=1:S
  
  ind=(s-1)*4;
  
  % dn ~ n & v
  M=M+sparse(ind+1,ind+1,kvds(s),NN,NN)+...  
     sparse(ind+1,ind+2,kx*ns0(s),NN,NN)+... 
     sparse(ind+1,ind+4,kz*ns0(s),NN,NN);

  % dv ~ n & v
  M=M+sparse(ind+2,ind+1,kx*csp(s)^2/ns0(s),NN,NN)+...  
     sparse(ind+4,ind+1,kz*csz(s)^2/ns0(s),NN,NN)+...
     sparse(ind+2,ind+2,kvds(s),NN,NN)+...
     sparse(ind+3,ind+3,kvds(s),NN,NN)+...
     sparse(ind+4,ind+4,kvds(s),NN,NN)+...
     sparse(ind+3,ind+2,-1i*wcs(s),NN,NN)+...
     sparse(ind+2,ind+3,1i*wcs(s),NN,NN);

  % dv ~ E
  M=M+sparse(ind+2,SJ+1,1i*qs(s)/ms(s),NN,NN)+...
     sparse(ind+3,SJ+2,1i*qs(s)/ms(s),NN,NN)+...
     sparse(ind+4,SJ+3,1i*qs(s)/ms(s),NN,NN);

  % dv ~ B, 20-12-25 12:27 update
%   M=M+sparse(ind+2,SJ+4,kz*(PszB0(s)-PspB0(s))/rhoms(s),NN,NN)+...
%      sparse(ind+2,SJ+5,-1i*qs(s)/ms(s)*vdsz(s),NN,NN)+...
%      sparse(ind+2,SJ+6,1i*qs(s)/ms(s)*vdsy(s)+...
%       (ifluidmodel==2)*kx*PspB0(s)*(gammaps(s)-1)/rhoms(s),NN,NN)+...
%      sparse(ind+3,SJ+4,1i*qs(s)/ms(s)*vdsz(s),NN,NN)+...
%      sparse(ind+3,SJ+5,kz*(PszB0(s)-PspB0(s))/rhoms(s),NN,NN)+...
%      sparse(ind+3,SJ+6,-1i*qs(s)/ms(s)*vdsx(s),NN,NN)+...
%      sparse(ind+4,SJ+4,-1i*qs(s)/ms(s)*vdsy(s)+...
%       kx*(PszB0(s)-PspB0(s))/rhoms(s)*2,NN,NN)+... % *2
%      sparse(ind+4,SJ+5,1i*qs(s)/ms(s)*vdsx(s),NN,NN)+...
%      sparse(ind+5,SJ+6,-(ifluidmodel==2)*kz*PszB0(s)*(gammazs(s)-1)/rhoms(s),NN,NN);
  M=M+sparse(ind+2,SJ+4,kz*(PszB0(s)-PspB0(s))/rhoms(s),NN,NN)+...
     sparse(ind+2,SJ+5,-1i*qs(s)/ms(s)*vdsz(s),NN,NN)+...
     sparse(ind+2,SJ+6,1i*qs(s)/ms(s)*vdsy(s)+...
      (ifluidmodel==2)*kx*PspB0(s)*(gammaps(s)-1)/rhoms(s),NN,NN)+...
     sparse(ind+3,SJ+4,1i*qs(s)/ms(s)*vdsz(s),NN,NN)+...
     sparse(ind+3,SJ+5,kz*(PszB0(s)-PspB0(s))/rhoms(s),NN,NN)+...
     sparse(ind+3,SJ+6,-1i*qs(s)/ms(s)*vdsx(s),NN,NN)+...
     sparse(ind+4,SJ+4,-1i*qs(s)/ms(s)*vdsy(s)+...
      kx*(PszB0(s)-PspB0(s))/rhoms(s)*1,NN,NN)+... % *2
     sparse(ind+4,SJ+5,1i*qs(s)/ms(s)*vdsx(s),NN,NN)+...
     sparse(ind+5,SJ+6,-(ifluidmodel==2)*kz*PszB0(s)*(gammazs(s)-1)/rhoms(s),NN,NN);
 
 
%   if(ifluidmodel==2) % fluid closure model: =2 double polytrope
%     csp=sqrt(Psp./rhoms); % Psp=csp2.*rhos0;
%   else %  =1, Xie2014 adiabatic
%     csp=sqrt(gammaps.*Psp./rhoms); % Psp=csp2.*rhos0;
%   end

  % dE ~ n
  M=M+sparse(SJ+1,ind+1,-1i*qs(s)*vdsx(s)/epsilon0,NN,NN)+...
     sparse(SJ+2,ind+1,-1i*qs(s)*vdsy(s)/epsilon0,NN,NN)+...
     sparse(SJ+3,ind+1,-1i*qs(s)*vdsz(s)/epsilon0,NN,NN);
  % dE ~ v
  M=M+sparse(SJ+1,ind+2,-1i*qs(s)*ns0(s)/epsilon0,NN,NN)+...
     sparse(SJ+2,ind+3,-1i*qs(s)*ns0(s)/epsilon0,NN,NN)+...
     sparse(SJ+3,ind+4,-1i*qs(s)*ns0(s)/epsilon0,NN,NN);
 
end

% E(B)
M=M+sparse(SJ+1,SJ+5,c2*kz,NN,NN)+...
  sparse(SJ+2,SJ+4,-c2*kz,NN,NN)+...
  sparse(SJ+2,SJ+6,c2*kx,NN,NN)+...
  sparse(SJ+3,SJ+5,-c2*kx,NN,NN);

% B(E)
M=M+sparse(SJ+4,SJ+2,-kz,NN,NN)+...
  sparse(SJ+5,SJ+1,kz,NN,NN)+...
  sparse(SJ+5,SJ+3,-kx,NN,NN)+...
  sparse(SJ+6,SJ+2,kx,NN,NN);

% % -- main program end of set the matrix elements --
