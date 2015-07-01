function resault =INCIDENTFIELD(E0,lemat,r,theta,fea,x0,y0,z0,l,w0)
%% �������䳡�ĵ糡ǿ��
% E0 ���
% lemat ����
% r ����İ뾶
% (x0,y0,z0)  ������������
% w0 �����뾶
% l  ��ɢ����
% theta,fea �Ƕ�

k = 2*pi/lemat;
R = k*r;
Nmax = R +4.05*R^(1/3)+2;

EriTemp02 = 0; %�ڶ���Eri��ʱֵ
ErthetaTemp02 = 0; %�ڶ���theta��ʱֵ
Erifea02 = 0; %�ڶ���fea��ʱֵ
for n=1:Nmax
    Cnpw = 1/(1i*k)*(-1i)^n*(2*n+1)/(n*(n+1));
    Feanr = sqrt(pi*R/2)*besselj(n+1/2,R);
    Feanrp =pi^(1/2)*(besselj(n - 1/2, R) - (besselj(n + 1/2, R)*(n + 1/2))/R)*(R/2)^(1/2) + (pi^(1/2)*besselj(n + 1/2, R))/(4*(R/2)^(1/2));
    EriTemp01 = 0;  %��һ��Eri��ʱֵ
    ErthetaTemp01 = 0; %��һ��theta��ʱֵ
    Erifea01 = 0; %��һ��fea��ʱֵ
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