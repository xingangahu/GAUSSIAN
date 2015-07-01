function resault = PITAUNMLEGENDRE(theta,Nstop,Mstop)
%% 直接利用勒让德函数计算角函数

cost = cos(theta);
sint = sin(theta);

if(cost==0)
    cost = 1e-5;
end
if(sint==0)
    sint = 1e-5;
end

if(Mstop>Nstop)
    disp('错误的参数输入Mstop>Nstop')
    resault =[0,0];
    return
end
tempLEGENDRE01 = legendre(Nstop,cost,'norm');

PINM =tempLEGENDRE01(1,1)/sint;

if(Nstop==Mstop)
    PINNM=0;
else
    tempLEGENDRE02 = legendre(Nstop-1,cost,'norm');
    PINNM =tempLEGENDRE02(1,1)/sint;
end

TAUNM = Nstop*cost*PINM-(Nstop+Mstop)*PINNM;

resault = [PINM,TAUNM];
end