function resault = INTERNALFILED(E0,lemat,r,theta,fea,x0,y0,z0,l,w0,complex,a)
%% 计算球体内部场的电场强度
% E0 振幅
% lemat 波长
% r 球体的半径
% (x0,y0,z0)  束腰中心坐标
% w0 束腰半径
% l  扩散长度
% theta,fea 角度
% complex 复折射率
% 这里磁导率取1

k = 2*pi/lemat;
R = k*r;
Nmax = R +4.05*R^(1/3)+2;



EraTemp02 = 0; %第二级Ers临时值
EathetaTemp02 = 0; %第二级theta临时值
Eafea02 = 0; %第二级fea临时值

% ABCD = Mie_abcd2_gaussin(Nmax,complex,k*a);

for n=1:Nmax
    Cnpw = 1/(1i*k)*(-1i)^n*(2*n+1)/(n*(n+1));
    Feanr = sqrt(pi*R/2)*besselj(n+1/2,R);
    Feanrp =pi^(1/2)*(besselj(n - 1/2, R) - (besselj(n + 1/2, R)*(n + 1/2))/R)*(R/2)^(1/2) + (pi^(1/2)*besselj(n + 1/2, R))/(4*(R/2)^(1/2));
    
%     CN =ABCD(3,n);
%     DN =ABCD(4,n);
    
     reasultan = debye(n,complex,k*a);
    
    %     CN =real(ABCD(3,n))-1i*imag(ABCD(3,n));
    %     DN =real(ABCD(3,n))-1i*imag(ABCD(3,n));
    
    EaTemp01 = 0;  %第一级Ers临时值
    EathetaTemp01 = 0; %第一级theta临时值
    Eafea01 = 0; %第一级fea临时值
    
    for m =-n:n
        
        pit_resault = PITAUNM(theta,n,abs(m));
        
        pinm = pit_resault(1);
        taunm = pit_resault(2);
        
        %         Cnm =CN*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,1);
        %         Dnm = DN*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,0);
        
       
        
        Cnm =reasultan(3)*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,1);
        Dnm = reasultan(4)*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,0);
        
        EaTemp01 = EaTemp01+Cnpw*n*(n+1)*...
            Feanr*Cnm*pinm*sin(theta)*exp(1i*m*fea);
        
        EathetaTemp01 = EathetaTemp01+...
            Cnpw*(Cnm*Feanrp*taunm+...
            m*(1/complex)*Dnm*Feanr*pinm)*exp(1i*m*fea);
        
        Eafea01 = Eafea01+Cnpw*(m*Cnm*Feanrp*pinm+...
            (1/complex)*Dnm*Feanr*taunm)*exp(1i*m*fea);
        
    end
    EraTemp02 =EraTemp02 +EaTemp01;
    EathetaTemp02 = EathetaTemp02 +EathetaTemp01;
    Eafea02 =Eafea02 + Eafea01;
end
Era =(E0/(R^2))*EraTemp02;
Eatheta = (E0/R)*EathetaTemp02;
Eafea = 1i*(E0/R)*Eafea02;

resault = [Era,Eatheta,Eafea];
