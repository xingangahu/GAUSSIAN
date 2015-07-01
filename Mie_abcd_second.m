function resault = Mie_abcd_second(nmax,m,x)
%% 计算an,bn,cn,dn
% m 复折射率
% x 尺度参数
nmax = fix(nmax);
n = (1:nmax);nu = (n+0.5);z = m.*x;m2 = m.*m;
sqx = sqrt(0.5*pi./x);sqz = sqrt(0.5*pi./z);
bx = besselj(nu,x).*sqx;
bz = besselj(nu,z).*sqz;
yx = bessely(nu,x).*sqx;
hx =bx+1i*yx;

h2x =bx-1i*yx;

b1x =[sin(x)/x,bx(1:nmax-1)];
b1z =[sin(z)/z,bz(1:nmax-1)];
y1x =[-cos(x)/x,yx(1:nmax-1)];
h1x = b1x+1i*y1x;

h12x = b1x-1i*y1x;

ax = x.*b1x-n.*bx;
az = z.*b1z-n.*bz;
ahx = x.*h1x-n.*hx;

a2hx = x.*h12x-n.*h2x;%递推公式不是很理解

an =(x.*bx.*az-m.*ax.*z.*bz)./(x.*h2x.*az-m.*a2hx.*z.*bz);
bn = (m.*x.*bx.*az-ax.*z.*bz)./(z.*h2x.*az-a2hx.*z.*bz);

cn = m.*(x.*h2x.*ax-a2hx.*x.*bx)./(x.*h2x.*az-m.*a2hx.*z.*bz);
dn = m2.*(x.*h2x.*ax-a2hx.*x.*bx)./(m.*x.*h2x.*az-a2hx.*z.*bz);

% resault =[an;bn;cn;dn];
resault =[real(an)-1i*imag(an);real(bn)-1i*imag(bn);real(cn)-1i*imag(cn);real(dn)-1i*imag(dn)];















