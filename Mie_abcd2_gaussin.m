function resault = Mie_abcd2_gaussin(nmax,m,x)
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
b1x =[sin(x)/x,bx(1:nmax-1)];
b1z =[sin(z)/z,bz(1:nmax-1)];
y1x =[-cos(x)/x,yx(1:nmax-1)];
h1x = b1x+1i*y1x;
ax = x.*b1x-n.*bx;
az = z.*b1z-n.*bz;
ahx = x.*h1x-n.*hx;
an =(m2.*bz.*ax-bx.*az)./(m2.*bz.*ahx-hx.*az);
bn =(bz.*ax-bx.*az)./(bz.*ahx-hx.*az);
cn =(bx.*ahx-hx.*ax)./(bz.*ahx-hx.*az);
dn =m.*(bx.*ahx-hx.*ax)./(m2.*bz.*ahx-hx.*az);
resault =[an;bn;cn;dn];