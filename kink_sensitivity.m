close all
clear all
file=["2400.mat","3100.mat","3400.mat","3700.mat","4500.mat"];
figure
hold on
for i = 1:length(file)
load(file(i))
t_vec_0=[0,t_vec];
plot(t_vec_0,P_block)
end
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([2000,4600])
title('Well Block Pressure')
legend("2400","3100","3400","3700","4500")
hold off