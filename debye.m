function resault = debye(n,m,x)

syms alpha beta p

TNX = sqrt(pi*x/2)*(besselj(n+0.5,x)+1i*bessely(n+0.5,x));
TNMX = sqrt(pi*m*x/2)*(besselj(n+0.5,m*x)+1i*bessely(n+0.5,m*x));
TNPX = (pi^(1/2)*(besselj(n + 1/2, x) + bessely(n + 1/2, x)*1i))/(4*(x/2)^(1/2)) +...
    pi^(1/2)*(x/2)^(1/2)*(besselj(n - 1/2, x) + bessely(n - 1/2, x)*1i - ...
    (besselj(n + 1/2, x)*(n + 1/2))/x - (bessely(n + 1/2, x)*(n + 1/2)*1i)/x);
TNPMX =pi^(1/2)*((m*x)/2)^(1/2)*(m*besselj(n - 1/2, m*x) + m*bessely(n - 1/2, m*x)*1i ...
    - (besselj(n + 1/2, m*x)*(n + 1/2))/x - (bessely(n + 1/2, m*x)*(n + 1/2)*1i)/x) +...
    (pi^(1/2)*m*(besselj(n + 1/2, m*x) + bessely(n + 1/2, m*x)*1i))/(4*((m*x)/2)^(1/2));

Epsnx= sqrt(x*pi/2)*(besselj(n+1/2,x)-1i*bessely(n+1/2,x));
Epsnmx= sqrt(m*x*pi/2)*(besselj(n+1/2,m*x)-1i*bessely(n+1/2,m*x));
Epsnxp =pi^(1/2)*(x/2)^(1/2)*(besselj(n - 1/2, x) - bessely(n - 1/2, x)*...
    1i - (besselj(n + 1/2, x)*(n + 1/2))/x + (bessely(n + 1/2, x)*(n + 1/2)*1i)/x) + ...
    (pi^(1/2)*(besselj(n + 1/2, x) - bessely(n + 1/2, x)*1i))/(4*(x/2)^(1/2));
Epsnmxp = pi^(1/2)*((m*x)/2)^(1/2)*(m*besselj(n - 1/2, m*x) - ...
    m*bessely(n - 1/2, m*x)*1i - (besselj(n + 1/2, m*x)*(n + 1/2))/x + ...
    (bessely(n + 1/2, m*x)*(n + 1/2)*1i)/x) + (pi^(1/2)*m*(besselj(n + 1/2, m*x)...
    - bessely(n + 1/2, m*x)*1i))/(4*((m*x)/2)^(1/2));

real(Epsnxp*TNMX)
-alpha*real(Epsnxp*TNMX)
beta*Epsnx*TNPMX
D = -alpha*Epsnxp*TNMX+beta*Epsnx*TNPMX;
RN212 =(alpha*TNPX*TNMX-beta*TNX*TNPMX)/D;
RN121 = (alpha*Epsnxp*Epsnmx-beta*Epsnx*Epsnmxp)/D;


%TM²¨
TN21 = 2*1i/D;
TN12 = m*2*1i/D;
TN12NUM = subs(TN12,[alpha beta],[m 1]);
TN21NUM = subs(TN21,[alpha beta],[m 1]);
RN212NUM = subs(RN212,[alpha beta],[m 1]);
RN121NUM = subs(RN121,[alpha beta],[m 1]);
an = double(0.5*(1-RN212NUM-symsum(TN21NUM*RN121NUM^(p-1)*TN12NUM,p,1,100)));
cn = double(m*symsum(TN21NUM*RN121NUM^(p-1),p,1,100));

%TE²¨
TN21 =m*2*1i/D;
TN12 = 2*1i/D;
TN12NUM = subs(TN12,[alpha beta],[1 m]);
TN21NUM = subs(TN21,[alpha beta],[1 m]);
RN212NUM = subs(RN212,[alpha beta],[1 m]);
RN121NUM = subs(RN121,[alpha beta],[1 m]);
bn = double(0.5*(1-RN212NUM-symsum(TN21NUM*RN121NUM^(p-1)*TN12NUM,p,1,100)));
dn = double((m^2)*symsum(TN21NUM*RN121NUM^(p-1),p,1,100));

resault = [an,bn,cn,dn];









