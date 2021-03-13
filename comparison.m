load('IMPES mod.mat')

figure
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
hold on
plot(t_vec_0,P_block)
set(gca,'fontsize',18)
% plot(BPR_1_t,BPR_1_p)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([2800,4000])
title('Well Block Pressure')

load('Implicit.mat')
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
plot(t_vec_0,P_block)
set(gca,'fontsize',18)
% plot(BPR_1_t,BPR_1_p)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([2800,4000])
title('Well Block Pressure')

load('implicit_benchmark.mat')
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
plot(t_vec_0,P_block)
set(gca,'fontsize',18)
plot(BPR_1_t,BPR_1_p)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([2800,4000])
title('Well Block Pressure')

legend('modified IMPES','implicit','implicit benchmark','Eclipse')
