% change this to 'p' for the first case and 'q' for the second case
control='p';
sum_time=0;

% number of A\b performed
num_slash=0;
while (t_total<info.t_runtime) %% && (p_well>info.p_switch)
    P_matrix_n1_k=P_matrix_n;
    % newton iteration
    % tol1=10^-5;
    tol1=10^-3;
    tol2=10^-2;
    tol3=10^-3;
    
    % fake condition to get into the loop
    loop=true;
    cond1=true;
    cond2=true;
    cond3=true;
    iter_newton=0;
    
    while  loop==true || (cond1 ==true) || (cond2 ==true) || (cond3 ==true)
        loop=false;
        cond1=false;
        cond2=false;
        cond3=false;
        % cut time step
        if iter_newton>iter_newton_max
            %     info.delta_t=info.delta_t/2;
            delta_t_n=delta_t_n/2;
            iter_newton=0;
            P_matrix_n1_k=P_matrix_n;
        end
        
        % make a connection list
        [Px_upwinded,Py_upwinded,Sx_upwinded,Sy_upwinded,Hx_oil,Hy_oil,Hx_gas,Hy_gas]=connection_list(info,P_matrix_n1_k);
        
        % make D matrix
        D_matrix_out=D_matrix(info,P_matrix_n,P_matrix_n1_k,delta_t_n);
        % D_matrix_out=D_matrix(info,P_matrix_n,P_matrix_n1_k);
        
        % make T matrix
        T_matrix_out=T_matrix(info,P_matrix_n1_k);
        
        % make Q matrix
        % fix k to be the geometric average
        Q_matrix_out=Q_matrix(info,P_matrix_n1_k,control);
        
        % make R matrix
        tic
        R_matrix=(T_matrix_out*P_matrix_n1_k-D_matrix_out*(P_matrix_n1_k-P_matrix_n)+Q_matrix_out);
        sum_time=sum_time+toc;
        
        % make J matrix
        J_matrix_out=J_matrix(info,P_matrix_n1_k,control,delta_t_n);
        % J_matrix_out=J_matrix(info,P_matrix_n1_k,control);
        tic
        delta_matrix=-J_matrix_out\R_matrix;
        num_slash=num_slash+1;
        P_matrix_n1_k1=P_matrix_n1_k+delta_matrix;
        sum_time=sum_time+toc;
        
        % normR=norm(R_matrix);
        % [P_res_2D_view,S_res_2D_view]=to2D(info,P_matrix_n1_k);
        % figure
        % imagesc(P_res_2D_view)
        % title('P')
        % colorbar
        % figure
        % imagesc(S_res_2D_view)
        % title('S')
        % colorbar
        
        [delta_matrix_2D_P,delta_matrix_2D_S]=to2D(info,delta_matrix);
        [P_matrix_n_P,P_matrix_n_S]=to2D(info,P_matrix_n);
        [P_matrix_n1_k_2D_P,P_matrix_n1_k_2D_S]=to2D(info,P_matrix_n1_k);
        [P_matrix_n1_k1_2D_P,P_matrix_n1_k1_2D_S]=to2D(info,P_matrix_n1_k1);
        
        % normR=max(abs((R_matrix)));
        
        % calculate R normalized matrix
        R_norm=zeros(2*info.Nx*info.Ny,1);
        for ii=1:info.Ny
            for jj=1:info.Nx
                l=(ii-1)*info.Nx+jj;
                p=P_matrix_n1_k1_2D_P(ii,jj);
                %         R_norm(2*l-1,1)=abs(5.615*Bg(info,p)*info.delta_t/(info.phi*info.V)*R_matrix(2*l-1,1));
                %         R_norm(2*l,1)=abs(5.615*Bo(info,p)*info.delta_t/(info.phi*info.V)*R_matrix(2*l,1));
                R_norm(2*l-1,1)=abs(5.615*Bg(info,p)*delta_t_n/(info.phi*info.V)*R_matrix(2*l-1,1));
                R_norm(2*l,1)=abs(5.615*Bo(info,p)*delta_t_n/(info.phi*info.V)*R_matrix(2*l,1));
            end
        end
        
        % cond1=max(abs((R_matrix)))>tol1;
        cond1=max(abs((R_norm)))>tol1;
        % this might be wrong, wrong place
        cond2=max(max(abs(P_matrix_n1_k1_2D_S-P_matrix_n1_k_2D_S)))>tol2;
        cond3=max(max(abs((P_matrix_n1_k1_2D_P-P_matrix_n1_k_2D_P)/mean(P_matrix_n1_k_2D_P))))>tol3;
        % cond3=max(max(abs((P_matrix_n1_k1_2D_P-P_matrix_n1_k_2D_P))))>tol3;
        
        
        % CFL matrix
        T_matrix_CFL_out=T_matrix_CFL(info,P_matrix_n1_k);
        % CFL_matrix=(T_matrix_CFL_out)*P_matrix_n1_k*info.delta_t/(info.phi*info.V);
        CFL_matrix=(T_matrix_CFL_out)*P_matrix_n1_k*delta_t_n/(info.phi*info.V);
        [CFL_matrix_G,CFL_matrix_O]=to2D(info,CFL_matrix);
        CFL_matrix_out=CFL_matrix_G+CFL_matrix_O;
        CFL_max=max(max(CFL_matrix_out));
        
        % update the P matrix for the next newton iteration
        P_matrix_n1_k=P_matrix_n1_k1;
        iter_newton=iter_newton+1;
        
    end % end newton iteration
    
    delta_s_max=max(max(abs(P_matrix_n1_k_2D_S-P_matrix_n_S)));
    delta_p_max=max(max(abs((P_matrix_n1_k_2D_P-P_matrix_n_P))));
    
    CFL_max_vec=[CFL_max_vec,CFL_max];
    iter_newton_vec=[iter_newton_vec,iter_newton];
    delta_t_vec=[delta_t_vec,delta_t_n];
    
    % turn this section off for part 1
    % check if the BHP drops below 2000 for case 2
    [Q_matrix_2D_Qg,Q_matrix_2D_Qo]=to2D(info,Q_matrix_out);
    p_wellblock=P_matrix_n1_k_2D_P(info.xw,info.yw);
    s_wellblock=P_matrix_n1_k_2D_S(info.xw,info.yw);
    qwo_wellblock=info.qwo_control;
    % qwo_wellblock=-Q_matrix_2D_Qo(info.xw,info.yw);
    To_check=info.WI*(kro(s_wellblock)/(info.mu_oil*Bo(info,p_wellblock)));
    p_well=p_wellblock-(qwo_wellblock)/To_check;
    %
    if control =='q'
        if (p_well>info.p_switch)
            control='q';
            P_BHP=[P_BHP,p_well];
        else
            control='p';
            P_BHP=[P_BHP,info.p_switch];
        end
    end
    
    % add data to plot
    [Q_matrix_2D_Qg,Q_matrix_2D_Qo]=to2D(info,Q_matrix_out);
    P_block=[P_block,P_matrix_n1_k_2D_P(info.xw,info.yw)];
    P_FPR=[P_FPR,mean(mean(P_matrix_n1_k_2D_P))];
    Qg_out=[Qg_out,-Q_matrix_2D_Qg(info.xw,info.yw)];
    Qo_out=[Qo_out,-Q_matrix_2D_Qo(info.xw,info.yw)];
    
    % save info to debug
    save("debug_plots\Implicit_"+num2str(t_total)+".mat");
%     disp(num2str(t_total));
    
    % auto-time stepping
    eta_s=0.05;
    eta_p=50;
    omega=0.5;
    % pressure
    delta_t_n1_p=delta_t_n*(1+omega)*eta_p/(delta_p_max+omega*eta_p);
    % saturation
    delta_t_n1_s=delta_t_n*(1+omega)*eta_s/(delta_s_max+omega*eta_s); % why these two are the same
    delta_t_n1_all=[info.delta_t_max,delta_t_n1_p,delta_t_n1_s];
    delta_t_n1=min(delta_t_n1_all);
    t_total=t_total+delta_t_n;
    
    t_vec=[t_vec, t_total];
    
    % update pressure and time for the next iteration
    % constant time step
%     delta_t_n=delta_t_n1;
    P_matrix_n=P_matrix_n1_k1;
end