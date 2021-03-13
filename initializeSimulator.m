function [delta_t_n,P_matrix_n,Pressure_matrix_n,Saturation_matrix_n,P_BHP,P_block,P_FPR,Qo_out,Qg_out,Qo_out_2,Qg_out_2,t_total,t_vec,CFL_max_vec,iter_newton_vec,delta_t_vec,p_well,iter_newton_max,R_vec]=initializeSimulator(info)
% initialize first time step
delta_t_n=info.delta_t;
% make the first initial guess P matrix
P_matrix_n =P_matrix(info.P_res_orig,info.S_res_orig);
Pressure_matrix_n=P_matrix_n(1:2:end-1);
Saturation_matrix_n=P_matrix_n(2:2:end);

% march pressure implicitly
% newton iteration set-up
P_BHP=[];
P_block=[info.Pinit];
P_FPR=[info.Pinit];
Qo_out=[];
Qg_out=[];
Qo_out_2=[];
Qg_out_2=[];

t_total=0;
t_vec=[];
CFL_max_vec=[];
iter_newton_vec=[];
delta_t_vec=[];
% fake condition to get into the loop
p_well=info.Pinit;
iter_newton_max=30;
% save residual
R_vec=[];
end