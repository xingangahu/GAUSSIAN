%% 验证 an，bn，cn，dn
% 参数说明
% m 复折射率
% x 尺度参数
% r 半径
% lemat 波长 单位为um
clc;clear

m = 1.33;
lemat = 0.6;
r = 5;
x = 2*pi*r/lemat;
nmax =100;

resaultan = zeros(1,100);
for n = 1:100
    resault =debye(n,m,x);
resaultan(n) =resault(4);
end

plot(1:100,real(resaultan),'b*-')
hold on
plot(1:100,imag(resaultan),'r*-')

% resault = Mie_abcd_second(50,m,x);
% % resault = Mie_abcd2_gaussin(100,m,x);
% 
% resault(:,50)
% [ran,col] =size(resault);
% 
% plot(1:col,real(resault(3,:))','b*-')
% hold on
% plot(1:col,imag(resault(3,:))','r*-')

title('an的实部与虚部与n的关系 采用debye级数获取')
ylabel('an 的实部与虚部')
xlabel('n') 
grid on