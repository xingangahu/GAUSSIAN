function resault = PITAUNM(theta,nStop,mStop)
%% 计算角函数
% theta :角度
% nStop,mStop :n,m的最大值
%% 变量声明
PInm =zeros(nStop+1,mStop+1);
TAUnm = zeros(nStop+1,mStop+1);
%%
Cost = cos(theta);
Sint = sin(theta);
if(mStop<0||nStop<0||mStop>nStop)
    disp('bad argument in pitaunm')
end
if(abs(Sint)<1e-99)
    Sint = 1e-99;
end

PInm(1,1) = 1/Sint;
TAUnm(1,1) = 0;
Factm = 1;
PInm(2,1) = Cost/Sint;
TAUnm(2,1) =Cost*PInm(2,1)-PInm(1,1);

if(nStop==1&&mStop==0)
    resault = [0,0];
    return
end

ntab = 3; %matlab中的元素坐标
for n=2:nStop
    PInm(ntab,1) = ((2*n-1)*Cost*PInm(ntab-1,1)-(n-1)*PInm(ntab-2,1))/n;
    TAUnm(ntab,1) = n*(Cost*PInm(ntab,1)-PInm(ntab-1,1));
    ntab =ntab+1;
end

mtab = 2; %matlab中的元素坐标
for m = 1:mStop
    PInm(mtab,mtab) =-PInm(mtab-1,mtab-1)*Factm*Sint;
    %     PInm(mtab,mtab) =-PInm(mtab-1,mtab-1)*(2*m+1)*abs(Sint);
    TAUnm(mtab,mtab) = m*Cost*PInm(mtab,mtab);
    %      TAUnm(mtab,mtab) = (m*Cost-2*m/(Cost*(2*m+1)))*PInm(mtab,mtab);
    PInm(mtab+1,mtab) =(2*m+1)*Cost*PInm(mtab,mtab);
    TAUnm(mtab+1,mtab) = (m+1)*Cost*PInm(mtab+1,mtab)-(2*m+1)*PInm(mtab,mtab);
    %     Factm = Factm+2;
    
    ntab = m+3; %matlab中的元素坐标
    for n =m+2:nStop
        PInm(ntab,mtab) = ((2*n-1)*Cost*PInm(ntab-1,mtab)-(n+m-1)*PInm(ntab-2,mtab))/(n-m);
        TAUnm(ntab,mtab) = n*Cost*PInm(ntab,mtab)-(n+m)*PInm(ntab-1,mtab);
        ntab =ntab+1;
    end
    
    mtab =mtab+1;
end

% PInm
% TAUnm
% nStop
% mStop
resault = [PInm(nStop+1,mStop+1),TAUnm(nStop+1,mStop+1)];

%save Pi.txt -ascii PInm
%save Tau.txt -ascii TAUnm


