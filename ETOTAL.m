%% 计算高斯波束下球体电场强度
% E0 电场振幅 lemat 波长
% r theta fea 计算区域的半径
% a 球体的半径
% x0,y0,z0 束腰中心半径
% l 偏振长度
% w0 束腰半径
% complex 复折射率
% 长度单位为um
% ETOTALINFILE 总场

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
    disp(['正在计算半径为' num2str(rtemp) '处的值']);
    
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
title('高斯波束下球形粒子散射特性')
shading interp
colorbar;








