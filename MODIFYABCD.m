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

% resault = Mie_abcd_second(nmax,m,x);
resault = Mie_abcd2_gaussin(nmax,m,x);

[ran,col] =size(resault);

plot(1:col,real(resault(3,:))','b*-')
hold on
plot(1:col,imag(resault(3,:))','r*-')

title('an��ʵ�����鲿��n�Ĺ�ϵ')
ylabel('an ��ʵ�����鲿')
xlabel('n') 
grid on