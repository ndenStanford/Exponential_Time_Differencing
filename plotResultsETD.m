%% plots

% [P_res_2D_view,S_res_2D_view]=to2D(info,P_matrix_n1_k);
figure (1)
% imagesc(P_res_2D_view)
% imagesc(Pressure_matrix_n1_2D)
imagesc(Pressure_matrix_n_2D)
% imagesc(P_matrix_n1_k_2D_P)
% surf(P_res_2D_view)
set(gca,'fontsize',18)
% title('Block Pressure Map (psi)')
title('Pressure (psi)')
colorbar
% surf(P_res_2D_view)

figure (2)
% imagesc(S_res_2D_view)
% imagesc(Saturation_matrix_n1_2D)
imagesc(Saturation_matrix_n_2D)
% imagesc(P_matrix_n1_k_2D_S)
set(gca,'fontsize',18)
% title('Gas Saturation Map')
title('Saturation')
colorbar

figure (3)
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

figure (4)
[row_CFL,col_CFL]=size(CFL_matrix_out);
imagesc(CFL_matrix_out(2:row_CFL-1,2:col_CFL-1))
set(gca,'fontsize',18)
title('CFL Map')
colorbar

% figure
% hold on
% plot(t_vec,iter_newton_vec)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% ylabel('Newton iteration number')
% title('Newton iteration number VS time')
% hold off

figure (5)
hold on
plot(t_vec,delta_t_vec)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('time step interval')
title('time step interval VS time')
hold off

% figure
% hold on
% plot(t_vec,P_BHP)
% plot(BHP_1_t,BHP_1_p)
% xlabel('time (days)')
% ylabel('pressure (psi)')
% title('bottom hole pressure')
% hold off

% figure (6)
% % t_vec=[0:1:max_iter]*info.delta_t;
% t_vec_0=[0,t_vec];
% hold on
% plot(t_vec_0,P_block,'LineWidth',2)
% set(gca,'fontsize',18)
% plot(BPR_1_t,BPR_1_p,'LineWidth',2)
% xlabel('time (days)')
% ylabel('pressure (psi)')
% % xlim([0,1600])
% xlim([0,600])
% ylim([2800,4000])
% legend('MATLAB output','Eclipse Output')
% legend boxoff
% % title('Well Block Pressure')
% hold off
% 
% figure (7)
% hold on
% % t_vec=[0:1:max_iter]*info.delta_t;
% plot(t_vec_0,P_FPR,'LineWidth',2)
% plot(FPR_1_t,FPR_1_p,'LineWidth',2)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% ylabel('pressure (psi)')
% % xlim([0,1600])
% xlim([0,600])
% ylim([3300,4500])
% legend('MATLAB output','Eclipse Output')
% legend boxoff
% % title('Field Pressure')
% hold off

figure (6)
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
hold on
plot(t_vec_0,P_block,'-','LineWidth',2)
set(gca,'fontsize',18)
plot(BPR_1_t,BPR_1_p,'-.','LineWidth',2)

plot(t_vec_0,P_FPR,':','LineWidth',2)
plot(FPR_1_t,FPR_1_p,'--','LineWidth',2)
xlabel('time (days)')
ylabel('pressure (psi)')
% xlim([0,1600])
xlim([0,600])
ylim([2800,4500])
legend('Well Block: ETDRK4','Well Block: Eclipse','Field : ETDRK4','Field : Eclipse')
legend boxoff
% title('Field Pressure')
hold off


% figure (8)
% % t_vec=[1:1:max_iter]*info.delta_t;
% hold on
% plot(t_vec,Qg_out*5.615/10^3,'LineWidth',2)
% plot(FGPR_1_t(2:length(FGPR_1_t)),FGPR_1_p(2:length(FGPR_1_p)),'LineWidth',2)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% ylabel('Qg (MSCF/days)')
% xlim([0,600])
% % xlim([0,1600])
% % ylim[(0.5*10^5,2.5*10^5)]
% legend('MATLAB output','Eclipse Output')
% legend boxoff
% % title('Gas flow rate: Qg')
% hold off
% 
% figure (9)
% % t_vec=[1:1:max_iter]*info.delta_t;
% hold on
% plot(t_vec,Qo_out,'LineWidth',2)
% plot(FOPR_1_t(2:length(FOPR_1_t)),FOPR_1_p(2:length(FOPR_1_p)),'LineWidth',2)
% set(gca,'fontsize',18)
% % xlim([0,1600])
% xlim([0,600])
% ylim([2000,8000])
% xlabel('time (days)')
% ylabel('Qo (STB/days)')
% legend('MATLAB output','Eclipse Output')
% legend boxoff
% % title('Oil flow rate: Qo')
% hold off

figure (8)
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
yyaxis left
plot(t_vec,Qg_out*5.615/10^3,'-','LineWidth',2)
plot(FGPR_1_t(2:length(FGPR_1_t)),FGPR_1_p(2:length(FGPR_1_p)),'-.','LineWidth',2)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('Qg (MSCF/days)')
xlim([0,600])

yyaxis right
plot(t_vec,Qo_out,':','LineWidth',2)
plot(FOPR_1_t(2:length(FOPR_1_t)),FOPR_1_p(2:length(FOPR_1_p)),'--','LineWidth',2)
set(gca,'fontsize',18)
% xlim([0,1600])
xlim([0,600])
ylim([2000,8000])
xlabel('time (days)')
ylabel('Qo (STB/days)')
legend('Qg: ETDRK4','Qg: Eclipse','Qo: ETDRK4','Qo: Eclipse')
legend boxoff
% title('Oil flow rate: Qo')
hold off
% ****************************************************************************

