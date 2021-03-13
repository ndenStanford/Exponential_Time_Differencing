close all
% clc
% clear all
set(gca,'FontSize', 18)
% get the ground truth solution
load("plot_ETDRK2_0.001.mat")
figure (1)
t_vec_0=[0,t_vec];
t_vec_0=t_vec_0(1:750002);
P_block=P_block(1:750002);
plot(t_vec_0,P_block,'DisplayName','ground truth','LineWidth',2)
n_interpl=10000;
t_vec_0_interpl_truth =linspace(0,t_vec_0(end),n_interpl);
% P_block_interpl_truth = interp1(t_vec_0,P_block,t_vec_0_interpl_truth);
P_block_interpl_truth = interp1(t_vec_0,P_block,t_vec_0_interpl_truth);
P_block_truth = P_block;
trapz_truth=trapz(P_block_interpl_truth);
% trapz_truth_600=trapz(P_block_interpl_truth(1:6001));
trapz_truth_600=trapz(P_block_interpl_truth(1:4001));
% trapz_truth_before_kink=trapz(P_block_interpl_truth(1:281));
trapz_truth_before_kink=trapz(P_block_interpl_truth(1:299));
color=['y','m','c','r','g','b','k'];
symbol=["o","s","d","p","o","s","d","p"];
line=["--",":","-.","-","--",":","-.","-","--",":","-."];
color_idx=1;

for select_L = ["high"]%,"mid","low"]
    if select_L=="mid" || select_L=="low"
        midfix="_"+select_L+"_";
    else
        midfix="_";
    end
    for scheme=["ETDRK4_Lie"]%["RK2","ETD2","ETDRK2","RK4","ETDRK4_Lie","ETDRK4_TF","Implicit"]
        k=1;
        time_toc_vec=[];
        time_step_vec=[];
        sum_time_vec=[];
        P_block_interpl={};
        
        if scheme=="ETDRK4_Lie"
            t_vec_loop=[1,2,3];%[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5,2,2.5,3,3.5];
        end
        if scheme=="ETDRK4_TF"
            t_vec_loop=[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.5,0.75,1,1.25];
        end
        if scheme=="RK4"
            t_vec_loop=[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3];
        end
        if scheme=="ETDRK2"
            t_vec_loop=[0.1,0.5,1];%[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.5,0.75,1,1.25,1.5,2];
        end
        if scheme=="ETD2"
            t_vec_loop=[0.5,0.75];%[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.3,0.4,0.5,0.75,1];
        end
        if scheme=="RK2"
            t_vec_loop=[0.15];%[0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175];
        end
        if scheme=="ETD1"
            t_vec_loop=[0.001,0.01,0.05,0.1,0.5,1,1.5];
        end
        if scheme=="IFRK2"
            t_vec_loop=0.01;%[0.01,0.05];
        end
        if scheme=="Implicit"
            t_vec_loop=[0.1,1,10,50];
        end
        if scheme=="IMPES"
            t_vec_loop=[0.1,0.5,1,2,5,10,50];
        end
        
        for t=t_vec_loop
            file_target="plot_"+scheme+midfix+num2str(t)+".mat";
            if exist(file_target, 'file') == 2
                load(file_target)
                % change select_L for plot legend
                scheme_plot=scheme;
                if scheme=="ETDRK4_Lie"
                    scheme_plot="ETDRK4-Lie";
                end
                if scheme=="ETDRK4_TF"
                    scheme_plot="ETDRK4-TF";
                end
                if scheme~="Implicit"
                    time_step_vec=[time_step_vec, info.delta_t];
                end
                if scheme=="Implicit"
                    time_step_vec=[time_step_vec, info.delta_t_max];
                end
                time_toc_vec=[time_toc_vec, time_toc];
                sum_time_vec=[sum_time_vec, sum_time];
                
%                 figure (1)
%                 % t_vec=[0:1:max_iter]*info.delta_t;
                t_vec_0=[0,t_vec];
%                 hold on
%                 plot(t_vec_0,P_block,'DisplayName',scheme_plot+" "+num2str(t)+" "+select_L,'LineWidth',2)
%                 scatter(t_vec_0,P_block,'DisplayName',scheme_plot+" "+num2str(t)+" "+select_L,'LineWidth',2)
%                 set(gca,'fontsize',18)
%                 %                 plot(BPR_1_t,BPR_1_p)
%                 xlabel('time (days)')
%                 ylabel('pressure (psi)')
%                 xlim([0,600])
%                 
%                 ylim([2800,4000])
%                 %title('Well Block Pressure')
%                 legend boxoff
%                 legend show
%                 % hold off
                
                %                 t_vec_0_interpl{k} = linspace(0,t_vec_0(end),n_interpl);
                t_vec_0_interpl{k} = linspace(0,750,n_interpl);
                P_block_interpl{k} = interp1(t_vec_0,P_block,t_vec_0_interpl{k});
                %                 figure
                %                 hold on
                %                 plot(t_vec_0_interpl{k},P_block_interpl{k})
                %                 plot(t_vec_0,P_block)
                %                 hold off
                k=k+1;
                
                
            end
        end
        
        trapz_diff=[];
        trapz_diff_before_kink=[];
        trapz_diff_600=[];
        
        for i=1:length(P_block_interpl)
            figure (1)
            set(gca,'FontSize', 18)
            hold on
            %             x=t_vec_0;
            %             y=P_block;
            %             xi = linspace(0,max(x),150);
            %             xi = linspace(0,750,length(P_block_interpl{i}));
            %             yi = interp1(x,y,xi,'linear');
            %             plot(xi,yi,'pr')
            
            %scatter(linspace(0,750,length(P_block_interpl{i})),P_block_interpl{i},'DisplayName',scheme_plot+" "+t_vec_loop(i)+" "+select_L);
            plot(linspace(0,750,length(P_block_interpl{i})),P_block_interpl{i},line(color_idx),'DisplayName',scheme_plot+" "+t_vec_loop(i)+" "+select_L,'LineWidth',3);
            xlim([0,600])
            ylim([2800,4000])
            legend boxoff
            legend show
            xlabel('time (days)','FontSize', 18)
            ylabel('pressure (psi)','FontSize', 18)
            
            diff=abs(P_block_interpl{i}-P_block_interpl_truth);
            %diff_normalized=diff./P_block_interpl_truth;
            diff_normalized=diff;
%             figure (5)
%             hold on
%             xxx=linspace(0,750,length(diff));
%             scatter(linspace(0,750,length(diff)),diff_normalized,'DisplayName',scheme_plot+" "+t_vec_loop(i)+" "+select_L)
%             %plot(linspace(0,750,length(diff)),diff_normalized,'DisplayName',scheme_plot+" "+t_vec_loop(i)+" "+select_L)
            trapz_diff=[trapz_diff,trapz(diff)];
%             % till day 600
%             %             trapz_diff_600=[trapz_diff_before_kink,trapz(diff(1:6001))];
            trapz_diff_600=[trapz_diff_600,trapz(diff(1:4001))];
%             %             trapz_diff_before_kink=[trapz_diff_before_kink,trapz(diff(1:281))];
            trapz_diff_before_kink=[trapz_diff_before_kink,trapz(diff(1:299))];
%             xlabel('time')
%             ylabel('diffrence from finiest timestep')
%             legend show
            color_idx=color_idx+1;
        end
        
        figure (7)
        set(gca,'FontSize', 12)
        hold on
        plot(time_step_vec,trapz_diff_before_kink./trapz_truth_before_kink,line(color_idx)+symbol(color_idx),'MarkerSize',10,'DisplayName',scheme_plot+" "+select_L,'LineWidth',3)
%         scatter(time_step_vec,trapz_diff_before_kink./trapz_truth_before_kink,color(color_idx),'d','filled','HandleVisibility','off')
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('time step','FontSize', 18)
        ylabel('normalized error','FontSize', 18)
%         title('Error plot of before kink','FontSize', 18)
        legend('Location','southeast')
        legend boxoff
        legend show
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 18)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 18)
        hold off
        
        figure (8)
        set(gca,'FontSize', 12)
        hold on
        plot(time_step_vec,trapz_diff_600./trapz_truth_600,line(color_idx)+symbol(color_idx),'MarkerSize',10,'DisplayName',scheme_plot+" "+select_L,'LineWidth',3)
%         scatter(time_step_vec,trapz_diff_600./trapz_truth_600,color(color_idx),'filled','d','HandleVisibility','off')
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('time step','FontSize', 18)
        ylabel('normalized error','FontSize', 18)
%         title('Error plot of entire simulation','FontSize', 18)
        legend('Location','southeast')
        legend boxoff
        legend show
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 18)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 18)
        hold off
        
        figure (9)
        set(gca,'FontSize', 12)
        hold on
        plot(sum_time_vec,trapz_diff_600./trapz_truth_600,line(color_idx)+symbol(color_idx),'MarkerSize',10,'DisplayName',scheme_plot+" "+select_L,'LineWidth',3)
%         scatter(sum_time_vec,trapz_diff_600./trapz_truth_600,'filled','d','HandleVisibility','off')
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('total solver computation time','FontSize', 18)
        ylabel('normalized error','FontSize', 18)
        xlim([5*10^-1,10^3])
%         title('Error plot of entire simulation','FontSize', 18)
        legend('Location','southeast')
        legend boxoff
        legend show
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 18)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 18)
        hold off
        
        figure (10)
        set(gca,'FontSize', 12)
        hold on
        plot(time_step_vec,sum_time_vec,line(color_idx)+symbol(color_idx),'MarkerSize',10,'DisplayName',scheme_plot+" "+select_L,'LineWidth',3)
%         scatter(time_step_vec,sum_time_vec,'filled','d','HandleVisibility','off')
        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('time step','FontSize', 18)
        ylabel('total solver computation time','FontSize', 18)
        legend('Location','southeast')
        legend boxoff
        legend show
        xAX = get(gca,'XAxis');
        set(xAX,'FontSize', 18)
        yAX = get(gca,'YAxis');
        set(yAX,'FontSize', 18)
        hold off
        
        %         figure (9)
        %         set(gca,'FontSize', 18)
        %         hold on
        %         set(gca, 'XScale', 'log', 'YScale', 'log');
        %         plot(time_step_vec(1:end),trapz_diff./trapz_truth,'DisplayName',scheme+" "+select_L)
        %         xlabel('time step')
        %         ylabel('normalized error')
        %         %title('Error plot of simulation before kink')
        %         legend('Location','southeast')
        %         legend boxoff
        %         legend show
        %         hold off
        color_idx=color_idx+1;
    end
end



