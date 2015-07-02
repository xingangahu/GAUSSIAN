function resault = SCATTEREDFIELD(E0,lemat,r,theta,fea,x0,y0,z0,l,w0,complex,a)
%% ����ɢ�䳡�ĵ糡ǿ��
% E0 ���
% lemat ����
% a ����İ뾶 r ����뾶
% (x0,y0,z0)  ������������
% w0 �����뾶
% l  ��ɢ����
% theta,fea �Ƕ�
% complex ��������

k = 2*pi/lemat;
R = k*r;
Nmax = R +4.05*R^(1/3)+2;

ErsTemp02 = 0; %�ڶ���Ers��ʱֵ
EsthetaTemp02 = 0; %�ڶ���theta��ʱֵ
Esfea02 = 0; %�ڶ���fea��ʱֵ

% ABCD = Mie_abcd2_gaussin(Nmax,complex,k*a);

for n=1:Nmax
    Cnpw = 1/(1i*k)*(-1i)^n*(2*n+1)/(n*(n+1));
    Epsnr = sqrt(R*pi/2)*(besselj(n+1/2,R)-1i*bessely(n+1/2,R));
    
    Epsnrp = sqrt(pi/2)*(1/2)*R^(-1/2)*(besselj(n+1/2,R)...
        -1i*bessely(n+1/2,R))+sqrt(pi*R/2)*0.5*...
        ((besselj(n-1/2,R)-besselj(n+1.5,R))-1i*(bessely(n-1/2,R)-bessely(n+1.5,R)));
    
%     Epsnrp =pi^(1/2)*(R/2)^(1/2)*(besselj(n - 1/2, R) - bessely(n - 1/2, R)*...
%         1i - (besselj(n + 1/2, R)*(n + 1/2))/R + (bessely(n + 1/2, R)*(n + 1/2)*1i)/R) + ...
%         (pi^(1/2)*(besselj(n + 1/2, R) - bessely(n + 1/2, R)*1i))/(4*(R/2)^(1/2));
    
%     AN =ABCD(1,n);
%     BN =ABCD(2,n);
    
    %     AN =-ABCD(1,n);
    %     BN =-ABCD(2,n);
    
    reasultan = debye(n,complex,k*a);
    
    EsiTemp01 = 0;  %��һ��Ers��ʱֵ
    EsthetaTemp01 = 0; %��һ��theta��ʱֵ
    Esfea01 = 0; %��һ��fea��ʱֵ
    
    for m =-n:n
        
        pit_resault = PITAUNM(theta,n,abs(m));
        
        pinm = pit_resault(1);
        taunm = pit_resault(2);        
        
        
        Anm =reasultan(1)*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,1);
        Bnm = reasultan(2)*BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,0);
        
        EsiTemp01 = EsiTemp01+Cnpw*n*(n+1)*...
            Epsnr*Anm*pinm*sin(theta)*exp(1i*m*fea);
        
        EsthetaTemp01 = EsthetaTemp01+...
            Cnpw*(Anm*Epsnrp*taunm+...
            m*Bnm*Epsnr*pinm)*exp(1i*m*fea);
        
        Esfea01 = Esfea01+Cnpw*(m*Anm*Epsnrp*pinm+...
            Bnm*Epsnr*taunm)*exp(1i*m*fea);
        
    end
    ErsTemp02 =ErsTemp02 +EsiTemp01;
    EsthetaTemp02 = EsthetaTemp02 +EsthetaTemp01;
    Esfea02 =Esfea02 + Esfea01;
end
Ers = (-E0/(R^2))*ErsTemp02;
Estheta =(-E0/R)*EsthetaTemp02;
Esfea = -1i*(E0/R)*Esfea02;

resault = [Ers,Estheta,Esfea];