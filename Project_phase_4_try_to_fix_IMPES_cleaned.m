%% initialize the reservoir
tic
format long
clear all
clc

[BHP_1_t,BHP_1_p,BPR_1_t,BPR_1_p,FPR_1_t,FPR_1_p,FOPR_1_t,FOPR_1_p,FGPR_1_t,FGPR_1_p]=importEclipse();
info=initializeParameters();
[delta_t_n,P_matrix_n,Pressure_matrix_n,Saturation_matrix_n,P_BHP,P_block,P_FPR,Qo_out,Qg_out,Qo_out_2,Qg_out_2,t_total,t_vec,CFL_max_vec,iter_newton_vec,delta_t_vec,p_well,iter_newton_max,R_vec]=initializeSimulator(info);

% change this to 'p' for the first case and 'q' for the second case
control='p';

while (t_total<info.t_runtime) %% && (p_well>info.p_switch)
    P_matrix_n1_k=P_matrix_n;
    Pressure_matrix_n1_k=P_matrix_n1_k(1:2:end-1);
    Saturation_matrix_n1_k=P_matrix_n1_k(2:2:end);
    % newton iteration
    tol1=10^-3;
    tol2=10^-2;
    tol3=10^-3;
    
    loop=true;
    cond1=false;
    cond3=false;
    iter_newton=0;
    
    while  loop==true || (cond1 ==false)  || (cond3 ==false)
        loop=false;
        % cut time step
        if iter_newton>iter_newton_max
            delta_t_n=delta_t_n/2;
            iter_newton=0;
            Pressure_matrix_n1_k=Pressure_matrix_n;
        end
        
        % make a connection list
        [Px_upwinded,Py_upwinded,Sx_upwinded,Sy_upwinded,Hx_oil,Hy_oil,Hx_gas,Hy_gas]=connection_list(info,P_matrix_n1_k);
        % make D matrix
        [D_matrix_out,d11_matrix,d12_matrix,d21_matrix,d22_matrix]=D_matrix_IMPES(info,P_matrix_n,Pressure_matrix_n1_k,delta_t_n);
        % make T matrix
        T_matrix_out=T_matrix_IMPES(info,P_matrix_n1_k,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
        
        % make Q matrix
        % fix k to be the geometric average
        Q_matrix_out=Q_matrix_IMPES(info,P_matrix_n1_k,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
        % make R matrix
        % R_matrix=(T_matrix_out*P_matrix_n1_k-D_matrix_out*(P_matrix_n1_k-P_matrix_n)+Q_matrix_out);
        R_matrix=(T_matrix_out*Pressure_matrix_n1_k-D_matrix_out*(Pressure_matrix_n1_k-Pressure_matrix_n)+Q_matrix_out);
        % make J matrix
        J_matrix_out=J_matrix_IMPES(info,P_matrix_n1_k,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix,delta_t_n);
        delta_matrix=-J_matrix_out\R_matrix;
        % P_matrix_n1_k1=P_matrix_n1_k+delta_matrix;
        Pressure_matrix_n1_k1=Pressure_matrix_n1_k+delta_matrix;
        
        % calculate R normalized matrix % now what to do?
        R_norm=zeros(2*info.Nx*info.Ny,1);
        for ii=1:info.Ny
            for jj=1:info.Nx
                l=(ii-1)*info.Nx+jj;
                R_norm(l,1)=abs(5.615*delta_t_n/(info.phi*info.V)*R_matrix(1,1));
            end
        end
        
        cond1=max(abs((R_norm)))<=tol1;
        cond3=max(max(abs((Pressure_matrix_n1_k1-Pressure_matrix_n1_k)/mean(Pressure_matrix_n1_k))))<tol3;
        
        % CFL matrix
        T_matrix_CFL_out=T_matrix_CFL(info,P_matrix_n1_k);
        CFL_matrix=(T_matrix_CFL_out)*P_matrix_n1_k*info.delta_t/(info.phi*info.V);
        [CFL_matrix_G,CFL_matrix_O]=to2D(info,CFL_matrix);
        CFL_matrix_out=CFL_matrix_G+CFL_matrix_O;
        CFL_max=max(max(CFL_matrix_out));
        
        % update the P matrix for the next newton iteration
        Pressure_matrix_n1_k=Pressure_matrix_n1_k1;
        iter_newton=iter_newton+1;
        
    end % end newton iteration
    
    % march saturation explicitly
    [D_matrix_out_Sat,d11_matrix,d12_matrix,d21_matrix,d22_matrix]=D_matrix_IMPES_Sat(info,P_matrix_n,Pressure_matrix_n1_k,delta_t_n);
    T_matrix_out_Sat=T_matrix_IMPES_Sat(info,P_matrix_n1_k,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
    Q_matrix_out_Sat=Q_matrix_IMPES_Sat(info,P_matrix_n1_k,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
    Saturation_matrix_n1_k1=Saturation_matrix_n1_k+(T_matrix_out_Sat*Pressure_matrix_n1_k-D_matrix_out_Sat*(Pressure_matrix_n1_k-Pressure_matrix_n)+Q_matrix_out_Sat);
    Pressure_matrix_n_2D=transpose(reshape(Pressure_matrix_n,[info.Nx,info.Ny]));
    Saturation_matrix_n_2D=transpose(reshape(Saturation_matrix_n,[info.Nx,info.Ny]));
    Pressure_matrix_n1_k1_2D=transpose(reshape(Pressure_matrix_n1_k1,[info.Nx,info.Ny]));
    Saturation_matrix_n1_k1_2D=transpose(reshape(Saturation_matrix_n1_k1,[info.Nx,info.Ny]));
    
    CFL_max_vec=[CFL_max_vec,CFL_max];
    iter_newton_vec=[iter_newton_vec,iter_newton];
    delta_t_vec=[delta_t_vec,delta_t_n];
    
    p_wellblock=Pressure_matrix_n1_k1_2D(info.xw,info.yw);
    s_wellblock=Saturation_matrix_n1_k1_2D(info.xw,info.yw);
    qwo_wellblock=info.qwo_control;
    To_check=info.WI*(kro(s_wellblock)/(info.mu_oil*Bo(info,p_wellblock)));
    p_well=p_wellblock-(qwo_wellblock)/To_check;
    
    if control =='q'
        % check if the BHP drops below 2000 for case 2
        if (p_well>info.p_switch)
            control='q';
            P_BHP=[P_BHP,p_well];
        else
            control='p';
            P_BHP=[P_BHP,info.p_switch];
        end
    end
   
    % add data to plot
    Q_matrix_all =Q_matrix_all_IMPES(info,p_wellblock,s_wellblock,control);
    [Q_matrix_2D_Qg,Q_matrix_2D_Qo]=to2D(info,Q_matrix_all);
    % P_block=[P_block,P_matrix_n1_k_2D_P(info.xw,info.yw)];
    % P_FPR=[P_FPR,mean(mean(P_matrix_n1_k_2D_P))];
    P_block=[P_block,Pressure_matrix_n1_k1_2D(info.xw,info.yw)];
    P_FPR=[P_FPR,mean(mean(Pressure_matrix_n1_k1_2D))];
    Qg_out=[Qg_out,-Q_matrix_2D_Qg(info.xw,info.yw)];
    Qo_out=[Qo_out,-Q_matrix_2D_Qo(info.xw,info.yw)];  
    
    % auto-time stepping
    % delta_t_n=info.delta_t;
    delta_t_max=1;
    eta_s=0.05;
    % eta_p=50;
    eta_p=200;
    omega=0.5;
    % pressure
    % delta_p_max=max(max(abs((P_matrix_n1_k_2D_P-P_matrix_n_P))));
    delta_p_max=max(max(abs((Pressure_matrix_n1_k1_2D-Pressure_matrix_n_2D))));
    delta_t_n1_p=delta_t_n*(1+omega)*eta_p/(delta_p_max+omega*eta_p);
    % saturation
    % delta_s_max=max(max(abs(P_matrix_n1_k_2D_S-P_matrix_n_S)));
    % delta_t_n1_s=delta_t_n*(1+omega)*eta_s/(delta_s_max+omega*eta_s);
    % delta_t_n1_all=[delta_t_max,delta_t_n1_p,delta_t_n1_s];
    % delta_t_n1=min(delta_t_n1_all);
    delta_t_n1_all=[delta_t_max,delta_t_n1_p];
    delta_t_n1=min(delta_t_n1_all);
    % info.delta_t=delta_t_n1; % fix this for implicit case too
    t_total=t_total+delta_t_n;
    t_vec=[t_vec, t_total];
    delta_t_n=delta_t_n1;
    R_vec=[R_vec,R_matrix];
    
    % update pressure for the next iteration
    % P_matrix_n=P_matrix_n1_k1;
    P_matrix_n(1:2:end-1)=Pressure_matrix_n1_k1;
    P_matrix_n(2:2:end)=Saturation_matrix_n1_k1;
    Pressure_matrix_n=P_matrix_n(1:2:end-1);
    Saturation_matrix_n=P_matrix_n(2:2:end);
end
toc

plotResults