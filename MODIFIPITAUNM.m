clc;clear

m = 1.33;
lemat = 0.6;
r = 5;
x = 2*pi*r/lemat;
theta = 10^(-8):pi/180:pi/2;
resaultpi = zeros(6,length(theta));
resaulttau = zeros(6,length(theta));

mflag = 1;
for ntemp = 1:6
    mflag = 1;
    for tempTheta = theta
%         pitaunm= PITAUNMLEGENDRE(tempTheta,ntemp,1);
         pitaunm=PITAUNM(tempTheta,ntemp,1);
        resaultpi(ntemp,mflag) = pitaunm(1);
        resaulttau(ntemp,mflag) = pitaunm(2);
        mflag = mflag+1;
    end
%     plot(theta,resaultpi(ntemp,:),'b*-')
    plot(theta,resaulttau(ntemp,:),'rp-')
    hold on
end


% plot(theta,resaultpi(1,:),'b*-')
% hold on
% plot(theta,resaultpi(2,:),'b*-')
% hold on
% plot(theta,resaultpi(3,:),'b*-')
% hold on
% plot(theta,resaultpi(4,:),'b*-')
% hold on
% plot(theta,resaultpi(5,:),'b*-')
% hold on


title('pi与theta的关系')
ylabel('pi')
xlabel('theta')
grid on