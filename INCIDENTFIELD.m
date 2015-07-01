function resault =INCIDENTFIELD(E0,lemat,r,theta,fea,x0,y0,z0,l,w0)
%% 计算入射场的电场强度
% E0 振幅
% lemat 波长
% r 球体的半径
% (x0,y0,z0)  束腰中心坐标
% w0 束腰半径
% l  扩散长度
% theta,fea 角度

k = 2*pi/lemat;
R = k*r;
Nmax = R +4.05*R^(1/3)+2;

EriTemp02 = 0; %第二级Eri临时值
ErthetaTemp02 = 0; %第二级theta临时值
Erifea02 = 0; %第二级fea临时值
for n=1:Nmax
    Cnpw = 1/(1i*k)*(-1i)^n*(2*n+1)/(n*(n+1));
    Feanr = sqrt(pi*R/2)*besselj(n+1/2,R);
    Feanrp =pi^(1/2)*(besselj(n - 1/2, R) - (besselj(n + 1/2, R)*(n + 1/2))/R)*(R/2)^(1/2) + (pi^(1/2)*besselj(n + 1/2, R))/(4*(R/2)^(1/2));
    EriTemp01 = 0;  %第一级Eri临时值
    ErthetaTemp01 = 0; %第一级theta临时值
    Erifea01 = 0; %第一级fea临时值
    for m =-n:n
        
        pit_resault = PITAUNM(theta,n,abs(m));  
        
        pinm = pit_resault(1);
        taunm = pit_resault(2);
        
        EriTemp01 = EriTemp01+Cnpw*n*(n+1)*...
            Feanr*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,1)...
            *pinm*sin(theta)*exp(1i*m*fea);
        
        ErthetaTemp01 = ErthetaTemp01+...
            Cnpw*(Feanrp*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,1)*taunm+...
            m*Feanr*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,0)*pinm)*exp(1i*m*fea);
        
        Erifea01 = Erifea01+Cnpw*(m*Feanrp*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,1)*pinm+...
            Feanr*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,0)*taunm)*exp(1i*m*fea);
        
    end
    EriTemp02 =EriTemp02 +EriTemp01;
    ErthetaTemp02 = ErthetaTemp02 +ErthetaTemp01;
    Erifea02 =Erifea02 + Erifea01;
end
Eri =(E0/(R^2))*EriTemp02;
Eitheta = (E0/R)*ErthetaTemp02;
Eifea = (1i*E0/R)*Erifea02;

resault = [Eri,Eitheta,Eifea];