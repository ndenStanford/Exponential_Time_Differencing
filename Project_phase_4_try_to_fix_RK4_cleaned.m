%% initialize the reservoir
format long
clear all
close all
clc

[BHP_1_t,BHP_1_p,BPR_1_t,BPR_1_p,FPR_1_t,FPR_1_p,FOPR_1_t,FOPR_1_p,FGPR_1_t,FGPR_1_p]=importEclipse();
info=initializeParameters();
% ETDRK4 Lie group goes up until 3.0
for scheme=["ETDRK2"]%["ETDRK2","ETD1","ETD2","IFRK2","RK2","ETDRK4_TF"]
    for delta_t_loop=[10]%[1]%[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.5]
        %explicit scheme
        % constant time step
        info.delta_t=delta_t_loop;
        %implicit scheme
        % variable time step
        %         info.delta_t_max=delta_t_loop;
        [delta_t_n,P_matrix_n,Pressure_matrix_n,Saturation_matrix_n,P_BHP,P_block,P_FPR,Qo_out,Qg_out,Qo_out_2,Qg_out_2,t_total,t_vec,CFL_max_vec,iter_newton_vec,delta_t_vec,p_well,iter_newton_max,R_vec]=initializeSimulator(info);
        tic
        select_L="high";
        mainSimulatorRK4
        %mainSimulatorImplicit
        %mainSimulatorIMPES
        time_toc=toc;
        if select_L=="mid" || select_L=="low"
            midfix="_"+select_L+"_";
        else
            midfix="_";
        end
%         save("plot_"+scheme+midfix+num2str(delta_t_loop)+".mat")
%         save("plot_Implicit"+midfix+num2str(delta_t_loop)+".mat")
%         save("plot_IMPES"+midfix+num2str(delta_t_loop)+".mat")
        plotResultsETD
%         plotResultsImplicit
%         disp(num_slash);
    end
end


