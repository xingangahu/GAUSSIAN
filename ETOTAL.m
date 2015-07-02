%% �����˹����������糡ǿ��
% E0 �糡��� lemat ����
% r theta fea ��������İ뾶
% a ����İ뾶
% x0,y0,z0 �������İ뾶
% l ƫ�񳤶�
% w0 �����뾶
% complex ��������
% ���ȵ�λΪum
% ETOTALINFILE �ܳ�

clc;clear;
E0 = 1;
lemat = 0.6328;
fea = 0;
x0 = 0 ;y0 = 0;z0 = 0;
w0 = 100;
l = (2*pi/lemat)*(w0^2);
complex =1.33;
a = 0.5;

r = 10^(-8):0.05:1;
theta = 10^(-8):pi/180:2*pi;

ETOTALINFILE = zeros(length(r),length(theta));

mflag = 1;
for rtemp = r
    nflag = 1;
    disp(['���ڼ���뾶Ϊ' num2str(rtemp) '����ֵ']);
    
    for thetatemp = theta
        if(rtemp>=a)            
            incident_resault =INCIDENTFIELD(E0,lemat,rtemp,thetatemp,fea,x0,y0,z0,l,w0);
            scater_resault = SCATTEREDFIELD(E0,lemat,rtemp,thetatemp,fea,x0,y0,z0,l,w0,complex,a);
        else
            intern_resault = INTERNALFILED(E0,lemat,rtemp,thetatemp,fea,x0,y0,z0,l,w0,complex,a);
        end
        
        if(rtemp<a)
            ETOTALINFILE(mflag,nflag) =abs(intern_resault)*(abs(intern_resault))';
        else
            ETOTALINFILE(mflag,nflag) =abs(incident_resault+scater_resault)*(abs(incident_resault+scater_resault))';
        end
        nflag =nflag+1;
        
        
    end
    mflag =mflag+1;
end

% xlswrite('d:\ETOTALINFILE.xls', ETOTALINFILE, 'sheet1');

% pcolor(r'*cos(theta),r'*sin(theta),ETOTALINFILE);
% shading interp
pcolor([r'*cos(theta),r'*cos(theta(1))],[r'*sin(theta),r'*sin(theta(1))],[ETOTALINFILE,ETOTALINFILE(:,1)]);
title('��˹��������������ɢ������')
shading interp
colorbar;








