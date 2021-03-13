% figure
% hold on
% load('case2_2000.mat')
% plot(t_vec,CFL_max_vec)
% % plot(CFL_max_vec)
% load('case2_orig.mat')
% plot(t_vec,CFL_max_vec)
% % plot(CFL_max_vec)
% load('case2_5000.mat')
% plot(t_vec,CFL_max_vec)
% % plot(CFL_max_vec)
% % load('case2_7000.mat')
% % plot(CFL_max_vec)
% 
% set(gca,'fontsize',18)
% xlabel('time (days)')
% xlim([0,1600])
% % xlabel('time step')
% ylabel('maximum CFL number')
% title('maximum CFL number VS time (days)')
% legend('2000 STB/day','3000 STB/day','5000 STB/day')
% % title('CFL number VS time step')
% hold off

% figure
% hold on
% load('case1_13.mat')
% plot(t_vec,CFL_max_vec)
% % plot(CFL_max_vec)
% load('case1_orig.mat')
% plot(t_vec,CFL_max_vec)
% % plot(CFL_max_vec)
% load('case1_45.mat')
% plot(t_vec,CFL_max_vec)
% % plot(CFL_max_vec)
% % load('case2_7000.mat')
% % plot(CFL_max_vec)
% 
% set(gca,'fontsize',18)
% xlabel('time (days)')
% xlim([0,1600])
% % xlabel('time step')
% ylabel('maximum CFL number')
% title('maximum CFL number VS time (days)')
% legend('13x13x1','23x23x1','45x45x1')
% % title('CFL number VS time step')
% hold off


% figure
% hold on
% load('case1_13.mat')
% plot(t_vec,iter_newton_vec)
% load('case1_orig.mat')
% plot(t_vec,iter_newton_vec)
% load('case1_45.mat')
% for i=1:length(iter_newton_vec)
%     if iter_newton_vec(i)>15
%         iter_newton_vec(i)=15;
%     end
% end
% plot(t_vec,iter_newton_vec)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% xlim([0,1600])
% ylabel('Newton iteration number')
% title('Newton iteration number VS time')
% legend('13x13x1','23x23x1','45x45x1')
% hold off

% figure
% hold on
% load('case2_2000.mat')
% plot(t_vec,iter_newton_vec)
% load('case2_orig.mat')
% plot(t_vec,iter_newton_vec)
% load('case2_5000.mat')
% % for i=1:length(iter_newton_vec)
% %     if iter_newton_vec(i)>15
% %         iter_newton_vec(i)=15;
% %     end
% % end
% plot(t_vec,iter_newton_vec)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% xlim([0,1600])
% ylabel('Newton iteration number')
% title('Newton iteration number VS time')
% legend('2000 STB/day','3000 STB/day','5000 STB/day')
% hold off

% % grid refinement
% figure
% hold on
% % t_vec=[0:1:max_iter]*info.delta_t;
% load('case2_t2.mat')
% plot(t_vec_0,P_FPR)
% load('case2_orig.mat')
% plot(t_vec_0,P_FPR)
% load('case1_t10.mat')
% plot(t_vec_0,P_FPR)
% plot(FPR_1_t,FPR_1_p)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% ylabel('pressure (psi)')
% % xlim([0,1600])
% % ylim([3300,4500])
% legend('MATLAB output','Eclipse Output')
% % legend('MATLAB 13x13x1','MATLAB 23x23x1','MATLAB 45x45x1','Eclipse')
% legend('max timestep = 2 days','max timestep = 5 days','max timestep = 10 days','Eclipse')
% title('Field Pressure')
% hold off
% 
% figure
% % t_vec=[1:1:max_iter]*info.delta_t;
% hold on
% load('case2_t2.mat')
% plot(t_vec,Qg_out*5.615/10^3)
% load('case2_orig.mat')
% plot(t_vec,Qg_out*5.615/10^3)
% load('case1_t10.mat')
% plot(t_vec,Qg_out*5.615/10^3)
% plot(FGPR_1_t(2:length(FGPR_1_t)),FGPR_1_p(2:length(FGPR_1_p)))
% set(gca,'fontsize',18)
% xlabel('time (days)')
% ylabel('Qg (MSCF/days)')
% % legend('MATLAB 13x13x1','MATLAB 23x23x1','MATLAB 45x45x1','Eclipse')
% legend('max timestep = 2 days','max timestep = 5 days','max timestep = 10 days','Eclipse')
% % xlim([0,1600])
% % ylim[(0.5*10^5,2.5*10^5)]
% title('Gas flow rate: Qg')
% hold off
% 
% figure
% % t_vec=[1:1:max_iter]*info.delta_t;
% hold on
% load('case2_t2.mat')
% plot(t_vec,Qo_out)
% load('case2_orig.mat')
% plot(t_vec,Qo_out)
% load('case1_t10.mat')
% plot(t_vec,Qo_out)
% plot(FOPR_1_t(2:length(FOPR_1_t)),FOPR_1_p(2:length(FOPR_1_p)))
% set(gca,'fontsize',18)
% xlim([0,1600])
% ylim([2000,3500])
% xlabel('time (days)')
% ylabel('Qo (STB/days)')
% % legend('MATLAB 13x13x1','MATLAB 23x23x1','MATLAB 45x45x1','Eclipse')
% legend('max timestep = 2 days','max timestep = 5 days','max timestep = 10 days','Eclipse')
% title('Oil flow rate: Qo')
% hold off
% 
% figure
% % t_vec=[0:1:max_iter]*info.delta_t;
% % t_vec=[0,t_vec];
% hold on
% load('case2_t2.mat')
% plot(t_vec,P_BHP)
% load('case2_orig.mat')
% plot(t_vec,P_BHP)
% load('case1_t10.mat')
% plot(t_vec,P_BHP)
% plot(BHP_1_t,BHP_1_p)
% set(gca,'fontsize',18)
% xlabel('time (days)')
% ylabel('pressure (psi)')
% title('bottom hole pressure')
% % legend('MATLAB 13x13x1','MATLAB 23x23x1','MATLAB 45x45x1','Eclipse')
% legend('max timestep = 2 days','max timestep = 5 days','max timestep = 10 days','Eclipse')
% hold off

% [P_res_2D_view,S_res_2D_view]=to2D(info,P_matrix_n1_k);
% figure
% load('case1_perm')
% imagesc(P_res_2D_view)
% surf(P_res_2D_view)
% set(gca,'fontsize',18)
% title('Block Pressure Map (psi)')
% colorbar

figure
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
hold on
load('run2.mat')
plot(t_vec_0,P_block)
set(gca,'fontsize',18)
% plot(BPR_1_t,BPR_1_p)
xlabel('time (days)')
ylabel('pressure (psi)')
legend('channel flow','original case 2')
% xlim([0,1600])
% ylim([2800,4000])
title('Well Block Pressure')
hold off

figure
hold on
% t_vec=[0:1:max_iter]*info.delta_t;
load('run2.mat')
plot(t_vec_0,P_FPR)
% load('case2_orig.mat')
% plot(t_vec_0,P_FPR)
% plot(FPR_1_t,FPR_1_p)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('pressure (psi)')
% xlim([0,1600])
% ylim([3300,4500])
legend('channel flow','original case 2')
title('Field Pressure')
hold off

figure
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
load('run2.mat')
plot(t_vec,Qg_out*5.615/10^3)
% load('case2_orig.mat')
% plot(t_vec,Qg_out*5.615/10^3)
% plot(FGPR_1_t(2:length(FGPR_1_t)),FGPR_1_p(2:length(FGPR_1_p)))
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('Qg (MSCF/days)')
legend('channel flow','original case 2')
% xlim([0,1600])
% ylim[(0.5*10^5,2.5*10^5)]
title('Gas flow rate: Qg')
hold off

figure
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
load('run2.mat')
plot(t_vec,Qo_out)
% load('case2_orig.mat')
% plot(t_vec,Qo_out)
% plot(FOPR_1_t(2:length(FOPR_1_t)),FOPR_1_p(2:length(FOPR_1_p)))
set(gca,'fontsize',18)
% xlim([0,1600])
% ylim([2000,3500])
xlabel('time (days)')
ylabel('Qo (STB/days)')
legend('channel flow','original case 2')
title('Oil flow rate: Qo')
hold off
