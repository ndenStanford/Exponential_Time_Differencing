close all
clc
clear all
% get the ground truth solution
load("plot_ETDRK2_0.001.mat")
figure (1)
t_vec_0=[0,t_vec];
% plot(t_vec_0,P_block,'ground truth')
n_interpl=10000;
t_vec_0_interpl_truth =linspace(0,t_vec_0(end),n_interpl);
P_block_interpl_truth = interp1(t_vec_0,P_block,t_vec_0_interpl_truth);
trapz_truth=trapz(P_block_interpl_truth);
trapz_truth_before_kink=trapz(P_block_interpl_truth(1:281));

for select_L = ["high","mid","low"]
    if select_L=="mid" || select_L=="low"
        midfix="_"+select_L+"_";
    else
        midfix="_";
    end
    for scheme=["ETDRK2"]%["ETDRK4_TF","ETDRK4_Lie","RK4","ETDRK2","ETD2","RK2","ETD1"]
        k=1;
        time_toc_vec=[];
        time_step_vec=[];
        P_block_interpl={};
        
        if scheme=="ETDRK4_Lie"
            t_vec_loop=[0.05,0.1,0.5,1.0,2.0,3.0];
        end
        if scheme=="ETDRK4_TF"
            t_vec_loop=[0.05,0.1,0.5,1.0];
        end
        if scheme=="RK4"
            t_vec_loop=[0.01,0.05,0.1,0.2,0.25,0.3];
        end
        if scheme=="ETDRK2"
            t_vec_loop=[0.001,0.01,0.05,0.1,0.2,0.4,0.5,1.2,1.5];
        end
        if scheme=="ETD2"
            t_vec_loop=[0.01,0.05,0.1,0.5];
        end
        if scheme=="RK2"
            t_vec_loop=[0.01,0.05,0.1];
        end
        if scheme=="ETD1"
            t_vec_loop=[0.001,0.01,0.05,0.1,0.5,1,1.5];
        end
        if scheme=="IFRK2"
            t_vec_loop=0.01;%[0.01,0.05];
        end
        if scheme=="Implicit"
            t_vec_loop=[0.01,0.1,1,1.5,2,5,50];
        end
        
        for t=t_vec_loop
            file_target="plot_"+scheme+midfix+num2str(t)+".mat";
            if exist(file_target, 'file') == 2
                %         load("plot_"+scheme+"_mid_"+num2str(t)+".mat")
                
                load(file_target)
                %         load("plot_"+scheme+"_"+num2str(t)+".mat")
                if scheme~="Implicit"
                    time_step_vec=[time_step_vec, info.delta_t];
                end
                if scheme=="Implicit"
                    time_step_vec=[time_step_vec, info.delta_t_max];
                end
                time_toc_vec=[time_toc_vec, time_toc];
                figure (1)
                % t_vec=[0:1:max_iter]*info.delta_t;
                t_vec_0=[0,t_vec];
                hold on
                plot(t_vec_0,P_block,'DisplayName',scheme+" "+num2str(t)+" "+select_L)
                set(gca,'fontsize',18)
                % plot(BPR_1_t,BPR_1_p)
                xlabel('time (days)')
                ylabel('pressure (psi)')
                xlim([0,1600])
                ylim([2800,4000])
                title('Well Block Pressure')
                legend show
                % hold off
                
                t_vec_0_interpl{k} =linspace(0,t_vec_0(end),n_interpl);
                P_block_interpl{k} = interp1(t_vec_0,P_block,t_vec_0_interpl{k});
                k=k+1;
            end
        end
        
        trapz_diff=[];
        trapz_diff_before_kink=[];
        for i=1:length(P_block_interpl)
            diff=abs(P_block_interpl{i}-P_block_interpl_truth);
            diff_normalized=diff./P_block_interpl_truth;
            figure (5)
            hold on
            plot(linspace(0,1500,length(diff)),diff_normalized,'DisplayName',scheme+" "+t_vec_loop(i)+" "+select_L)
            trapz_diff=[trapz_diff,trapz(diff)];
            trapz_diff_before_kink=[trapz_diff_before_kink,trapz(diff(1:281))];
            xlabel('time')
            ylabel('diffrence from finiest timestep')
            legend show
        end
        %
        %     figure (6)
        %     hold on
        %     plot(time_toc_vec(2:end),trapz_diff)
        %     time_toc_vec
        %     legend("ETD","Implicit")
        %     xlabel('total simulation time')
        %     ylabel('total diffrence from finiest timestep')
        %
        
        figure (7)
        hold on
        plot(time_step_vec(1:end),trapz_diff./trapz_truth,'DisplayName',scheme)
        %     loglog(time_step_vec(1:end),trapz_diff_before_kink./trapz_truth_before_kink,'DisplayName',scheme)
        xlabel('time step')
        ylabel('total diffrence from ground truth')
        % ylabel('total diffrence from ground truth before kink')
        legend show
        hold off
        
        figure (8)
        hold on
        plot(time_step_vec(1:end),trapz_diff_before_kink./trapz_truth_before_kink,'DisplayName',scheme)
        xlabel('time step')
        ylabel('total diffrence from ground truth before kink')
        legend show
        hold off
    end
end