[P_res_2D_view,S_res_2D_view]=to2D(info,P_matrix_n1_k);
figure
imagesc(P_res_2D_view)
% surf(P_res_2D_view)
set(gca,'fontsize',18)
title('Block Pressure Map (psi)')
colorbar
% surf(P_res_2D_view)


figure
imagesc(S_res_2D_view)
set(gca,'fontsize',18)
title('Gas Saturation Map')
colorbar

figure
hold on
plot(t_vec,CFL_max_vec)
% plot(CFL_max_vec)
set(gca,'fontsize',18)
% xlabel('time (days)')
xlabel('time step')
ylabel('CFL number')
% title('CFL number VS time')
title('CFL number VS time step')
hold off

figure
[row_CFL,col_CFL]=size(CFL_matrix_out);
imagesc(CFL_matrix_out(2:row_CFL-1,2:col_CFL-1))
set(gca,'fontsize',18)
title('CFL Map')
colorbar

figure
hold on
plot(t_vec,iter_newton_vec)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('Newton iteration number')
title('Newton iteration number VS time')
hold off

figure
hold on
plot(t_vec,delta_t_vec)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('time step interval')
title('time step interval VS time')
hold off


figure
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
hold on
plot(t_vec_0,P_block)
set(gca,'fontsize',18)
plot(BPR_1_t,BPR_1_p)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([2800,4000])
title('Well Block Pressure')
hold off

figure
hold on
% t_vec=[0:1:max_iter]*info.delta_t;
plot(t_vec_0,P_FPR)
plot(FPR_1_t,FPR_1_p)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([3300,4500])
legend('MATLAB output','Eclipse Output')
title('Field Pressure')
hold off

figure
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
plot(t_vec,Qg_out*5.615/10^3)
plot(FGPR_1_t(2:length(FGPR_1_t)),FGPR_1_p(2:length(FGPR_1_p)))
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('Qg (MSCF/days)')
% xlim([0,1600])
% ylim[(0.5*10^5,2.5*10^5)]
title('Gas flow rate: Qg')
hold off

figure
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
plot(t_vec,Qo_out)
plot(FOPR_1_t(2:length(FOPR_1_t)),FOPR_1_p(2:length(FOPR_1_p)))
set(gca,'fontsize',18)
xlim([0,1600])
ylim([0,5000])
xlabel('time (days)')
ylabel('Qo (STB/days)')
title('Oil flow rate: Qo')
hold off
% ****************************************************************************