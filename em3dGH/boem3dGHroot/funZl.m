function Z=funZl(z,l,N)
% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-05-26 23:37
% Calculate GPDF, Generalized Plasma Dispersion Function, see [Xie2013]
%
% Modify the parameters in the code for your own usage
%
% Ref: [1] H. S. Xie, Generalized Plasma Dispersion Function:
% One-Solve-All Treatment, Visualizations and Application to
% Landau Damping, PoP, 2013. (Or http://arxiv.org/abs/1305.6476)
%
% Examples:
%  1. z=-1:0.01:1; [Zp,Z]=zetaph(z); plot(z,real(Zp),z,imag(Zp));
%  2. z=-1:0.01:1; F = '1./(1+v.^2)/pi'; [Zp,Z]=zetaph(z,0,F);
%     plot(z,real(Z),z,imag(Z));
%  3. [x,y]=meshgrid(-1:0.1:1,-1:0.1:1); z=x+1i*y; [Zp,Z]=zetaph(z,1);
%     subplot(121);surf(x,y,real(Zp)); subplot(122);surf(x,y,imag(Zp));
% 24-12-28 07:58 modify for v^l*exp(-v^2)/sqrt(pi)
if nargin<3, N = 128; end

% Default for usual Plasma Dispersion Function
if (nargin<2)
    l=0;
end
F = ['v.^',num2str(l),'.*exp(-v.^2)/sqrt(pi)'];
Z=calZ(z,F,N);
end

function Z=calZ(z,F,N,del)
if nargin<4, del = 1; end
Z=hilb(z,F,N,del);
ind1=find(isnan(Z));
z1=z(ind1)+1e-10;         % avoid NaN of higher order singular point
Z(ind1)=hilb(z1,F,N,del);% e.g., z=ia for Lorentzian F = '1./(a^2+v.^2)'
end

function Z=hilb(z,F,N,del,L)
% note: 1. in fact, a_n need calcualte only once for a fixed F, so you can
%   rewrite the code to speed up.
%       2. f(z) in analytic continuation can also be replaced by
%   sum{a_n*rho_n}
%       3. usually, del=1, but for flat-top and triangular distributions,
%   it seems we should set del=0 for correct analytic continuation

if nargin<5, L = sqrt(N/sqrt(2)); end  % optimal choice of L
%     L=10;
if nargin<4, del = 1; end

% 1. Define initial parameters
Z = zeros(size(z)); % initialize output

% 2. Calculate
idx=find(imag(z) == 0);
idx1=find(imag(z) ~= 0);
% 2.1 real line
for id=idx
    n  = [-N:N-1]';                   % Compute the collocation points
    v  = L*tan(pi*(n+1/2)/(2*N));
    FF = eval(F);             % sample the function
    FF(isnan(FF))=0;
    Fz = eval(['@(v)',F]);    % define F(z), for analytic continuation

    % Weideman95 method to calculate the real line Hilbert transform
    a  = fft(fftshift(FF.*(L-1i*v)));     % These three lines compute
    a  = exp(-1i*n*pi/(2*N)).*fftshift(a);   % expansion coefficients
    a  = flipud(1i*(sign(n+1/2).*a))/(2*N);
    t  = (L+1i*z(id))./(L-1i*z(id)); % The evaluation of the transform
    h  = polyval(a,t)./(t.^N.*(L-1i*z(id)));   % reduces to polynomial
    % evaluation in the variable t
    Z(id) = h + 1i.*Fz(z(id));
end
% 2.2. upper and lower half plane
for id=idx1
    M = 2*N;  M2 = 2*M;
    k = [-M+1:1:M-1]';% M2 = no. of sampling points
    theta = k*pi/M; v = L*tan(theta/2);     % define theta & v
    FF = eval(F);             % sample the function
    FF(isnan(FF))=0;
    Fz = eval(['@(v)',F]);    % define F(z), for analytic continuation

    % Weideman94 method to calculate upper half plane Hilbert transform
    W = (L^2+v.^2);                         % default weight function
    FF = FF.*W; FF = [0; FF];             % function to be transformed
    a = (fft(fftshift(FF)))/M2;           % coefficients of transform
    a0 = a(1); a = flipud(a(2:N+1));            % reorder coefficients
    z1 = (imag(z(id))>0).*z(id)+(imag(z(id))<0).*conj(z(id));
    t = (L+1i*z1)./(L-1i*z1); p = polyval(a,t); % polynomial evaluation
    h = 1i*(2*p./(L-1i*z1).^2+(a0/L)./(L-1i*z1));  % Evaluate h(f(v))
    % convert upper half-plane results to lower half-plane if necesary
    %         Z(id) = h.*(imag(z(id))>0)+conj(h-2i.*Fz(z1)).*(imag(z(id))<0);
    Z(id) = h.*(imag(z(id))>0)+(conj(h)+...
        del*2i.*Fz(z(id))).*(imag(z(id))<0);
end
Z=Z.*pi;
end
