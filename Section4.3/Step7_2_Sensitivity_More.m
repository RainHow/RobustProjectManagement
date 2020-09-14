clear;
clc

FontSize = 50;
BigFontSize = 56;

for instance_type = 1:2
    eval(['instance' num2str(instance_type) ' = load(''DATA/DC' num2str(instance_type) '/InstanceAnalysis.mat'');']);
    eval(['result' num2str(instance_type) '= load(''DATA/DC' num2str(instance_type) '/CapacityAnalysis.mat'');']);
end

M1 = length(result1.Reservation_x);
M2 = length(result2.Reservation_x);

N = [instance1.N(1:M1); instance2.N(1:M2)];
N_pr = [instance1.N_pr(1:M1); instance2.N_pr(1:M2)];
Gamma = [instance1.Gamma(1:M1); instance2.Gamma(1:M2)];
LongestChain = [instance1.LongestChain(1:M1); instance2.LongestChain(1:M2)];
Reservation_x = [result1.Reservation_x(1:M1); result2.Reservation_x(1:M2)];
Reservation_y = [result1.Reservation_y(1:M1); result2.Reservation_y(1:M2)];

clear instance* result*

Amount_x = zeros(M1+M2,1);
Amount_y = zeros(M1+M2,1);

for i = 1 : (M1+M2)
    Amount_x(i) = [1 2 3 4 5] * double(Reservation_x{i});
    Amount_y(i) = sum(Reservation_y{i});
end

Gamma_XAmount = sortrows([Gamma Amount_x]);
N_XAmount = sortrows([N Amount_x]);
N_pr_XAmount = sortrows([N_pr Amount_x]);

Gamma_YAmount = sortrows([Gamma Amount_y]);
N_YAmount = sortrows([N Amount_y]);
N_pr_YAmount = sortrows([N_pr Amount_y]);

subplot(1,2,1)
plot(Gamma_XAmount(:,1),Gamma_XAmount(:,2))
title('General capacity reservation');
xlabel('\Gamma')
ylabel('Amount of general capacity')
set(gca,'FontSize',FontSize)
subplot(1,2,2)
plot(Gamma_YAmount(:,1),Gamma_YAmount(:,2))
title('Customized capacity reservation');
xlabel('\Gamma')
ylabel('Amount of specific capacity')
set(gca,'FontSize',FontSize)

eval([ 'saveas(gcf,''DATA/section5_3.fig'')' ])


% 
% 
% subplot(2,3,1)
% plot(Gamma_XAmount(:,1),Gamma_XAmount(:,2))
% title('Reservation-x on Gamma');
% xlabel('Gamma')
% ylabel('amount of x')
% set(gca,'FontSize',FontSize)
% 
% 
% subplot(2,3,2)
% plot(N_XAmount(:,1),N_XAmount(:,2))
% title('Reservation-x on Task');
% xlabel('# of tasks')
% ylabel('amount of x')
% set(gca,'FontSize',FontSize)
% 
% 
% subplot(2,3,3)
% plot(N_pr_XAmount(:,1),N_pr_XAmount(:,2))
% title('Reservation-x on Precedence');
% xlabel('# of precedence relationship')
% ylabel('amount of x')
% set(gca,'FontSize',FontSize)
% 
% 
% subplot(2,3,4)
% plot(Gamma_YAmount(:,1),Gamma_YAmount(:,2))
% title('Reservation-y on Gamma');
% xlabel('Gamma')
% ylabel('amount of y')
% set(gca,'FontSize',FontSize)
% 
% 
% subplot(2,3,5)
% plot(N_YAmount(:,1),N_YAmount(:,2))
% title('Reservation-y on Task');
% xlabel('# of tasks')
% ylabel('amount of y')
% set(gca,'FontSize',FontSize)
% 
% 
% subplot(2,3,6)
% plot(N_pr_YAmount(:,1),N_pr_YAmount(:,2))
% title('Reservation-y on Precedence');
% xlabel('# of precedence relationship')
% ylabel('amount of y')
% set(gca,'FontSize',FontSize)