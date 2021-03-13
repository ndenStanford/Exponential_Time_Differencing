%% initialize the reservoir
format long
clear all
close all
clc

[BHP_1_t,BHP_1_p,BPR_1_t,BPR_1_p,FPR_1_t,FPR_1_p,FOPR_1_t,FOPR_1_p,FGPR_1_t,FGPR_1_p]=importEclipse();
info=initializeParameters();
scheme="ETDRK2";%["ETDRK2","ETD1","ETD2","IFRK2","RK2"];
for delta_t_loop=1%[0.01,0.05,0.1,0.5,1,1.5,2,10]
    info.delta_t=delta_t_loop;
    [delta_t_n,P_matrix_n,Pressure_matrix_n,Saturation_matrix_n,P_BHP,P_block,P_FPR,Qo_out,Qg_out,Qo_out_2,Qg_out_2,t_total,t_vec,CFL_max_vec,iter_newton_vec,delta_t_vec,p_well,iter_newton_max,R_vec]=initializeSimulator(info);
    tic
    mainSimulatorRK2
    time_toc=toc;
    save("plot_"+scheme+"_"+num2str(delta_t_loop)+".mat")
    plotResultsETD
end 


