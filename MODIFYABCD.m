%% ��֤ an��bn��cn��dn
% ����˵��
% m ��������
% x �߶Ȳ���
% r �뾶
% lemat ���� ��λΪum
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

title('an��ʵ�����鲿��n�Ĺ�ϵ ����debye������ȡ')
ylabel('an ��ʵ�����鲿')
xlabel('n') 
grid on