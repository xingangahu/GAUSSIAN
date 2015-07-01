function resault=BSCOEFFICIENTS(x0,y0,z0,l,k,w0,m,n,flag)
%% �����˹������״ϵ��
% (x0,y0,z0) :������������
% l: ��ɢ����
% w0:�����뾶
% m,n:�������
% flag :TE,TM����� 0����TE�� 1����TM��
Q =1/(1i-2*z0/l);
Znm = (-2*1i/((2*n+1)))^(abs(m)-1);
S = 1/(k*w0);
RouN=(n+1/2)*S;
Eps0 =x0/w0;
Eta0 = y0/w0;
Rou0 = sqrt(Eps0^2+Eta0^2);
Gamma = sqrt(RouN^2+Rou0^2);
if(Eta0==0)
    Fea0=0;
else
    Fea0 =Eps0/Eta0;
end

switch flag
    case 0 %TE��        
        resault =Q*Znm*exp(-1i*Q*Gamma^2+1i*k*z0)...
            *(1/2)*(exp(1i*(m-1)*Fea0)*besselj(m-1,2*Q*RouN*Rou0)...
            -exp(1i*(m+1)*Fea0)*besselj(m+1,2*Q*RouN*Rou0));
    case 1 %TM��
        resault =1i*Q*Znm*exp(-1i*Q*Gamma^2+1i*k*z0)...
            *(1/2)*(exp(1i*(m-1)*Fea0)*besselj(m-1,2*Q*RouN*Rou0)...
            +exp(1i*(m+1)*Fea0)*besselj(m+1,2*Q*RouN*Rou0));
end
