% change this to 'p' for the first case and 'q' for the second case
% for i=1:max_iter
control='p';
% fake condition to get into the loop
delta_p_max=0;
% load L0
load('step0.mat')
% after setting L0 to a constant then can pull many terms out of
% time stepping loop
delta_t_n=info.delta_t;
pinv_L=pinv(L0);
exp_Lh=expm(L0*delta_t_n);
[row,col]=size(L0);
I=eye(row,col);
M1_1=pinv_L*(exp_Lh-I);
M1_2=exp_Lh*(I-pinv_L*L0);
M2_1=pinv_L*pinv_L;
M2_2=(1/2)*exp_Lh*(I-pinv_L*L0);
while (t_total<info.t_runtime) %% && (p_well>info.p_switch)
    delta_t_n=info.delta_t;
    loop=true;
    %     while (delta_p_max > 10 && t_total > 40) || loop
    loop=false;
    % P_matrix_n1_k=P_matrix_n;
    Pressure_matrix_n=P_matrix_n(1:2:end-1);
    Saturation_matrix_n=P_matrix_n(2:2:end);
    
    if (scheme=="ETDRK2")
        % ETDRK2 : pseudoinverse
        % L=L_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n);
        M1=M1_1+delta_t_n*M1_2;
        M2=M2_1*(exp_Lh-(I+L0*delta_t_n))+delta_t_n^2*M2_2;
        an=exp_Lh*Pressure_matrix_n+M1*N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        Pressure_matrix_n1=an+delta_t_n^-1*M2*(N_matrix_IMPES(info,an,Saturation_matrix_n,delta_t_n,control,L0)-N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        
    elseif (scheme=="Euler")
        % Euler
        Pressure_matrix_n1=Pressure_matrix_n+delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        
    elseif (scheme=="RK2")
        % RK2
        k1=delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        k2=delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n+k1,Saturation_matrix_n,delta_t_n,control,L0));
        Pressure_matrix_n1=Pressure_matrix_n+(k1+k1)/2;
        
    elseif (scheme=="ETD1")
        % ETD1
        M1=M1_1+delta_t_n*M1_2;
        Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+M1*Fn;
        
    elseif (scheme=="ETD2")
        % first time step need to do ETD1 because there is no Fn-1
        if t_total==0
            % ETD1
            M1=M1_1+delta_t_n*M1_2;
            Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
            Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+M1*Fn;
        else           
            % ETD2
            M1=M1_1+delta_t_n*M1_2;
            M2=M2_1*(exp_Lh-(I+L0*delta_t_n))+delta_t_n^2*M2_2;
            Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
            Fn_minus_1=N_matrix_IMPES(info,Pressure_matrix_previous,Saturation_matrix_previous,delta_t_n,control,L0);
            Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+M1*Fn+M2*(Fn-Fn_minus_1)/delta_t_n;
        end
        
    elseif (scheme=="IFRK2")
        % IFRK2
        Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        un_prime=exp_Lh*(Pressure_matrix_n+delta_t_n*Fn);
        Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+(delta_t_n/2)*(exp_Lh*Fn+N_matrix_IMPES(info,un_prime,Saturation_matrix_n,delta_t_n,control,L0));
    end
    
    % march saturation explicitly
    [D_matrix_out_Sat,d11_matrix,d12_matrix,d21_matrix,d22_matrix]=D_matrix_IMPES_Sat(info,P_matrix_n,Pressure_matrix_n1,delta_t_n);
    T_matrix_out_Sat=T_matrix_IMPES_Sat(info,P_matrix_n,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
    Q_matrix_out_Sat=Q_matrix_IMPES_Sat(info,P_matrix_n,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix);
    Saturation_matrix_n1=Saturation_matrix_n+(T_matrix_out_Sat*Pressure_matrix_n1-D_matrix_out_Sat*(Pressure_matrix_n1-Pressure_matrix_n)+Q_matrix_out_Sat);
    Pressure_matrix_n_2D=transpose(reshape(Pressure_matrix_n,[info.Nx,info.Ny]));
    Saturation_matrix_n_2D=transpose(reshape(Saturation_matrix_n,[info.Nx,info.Ny]));
    Pressure_matrix_n1_2D=transpose(reshape(Pressure_matrix_n1,[info.Nx,info.Ny]));
    Saturation_matrix_n1_2D=transpose(reshape(Saturation_matrix_n1,[info.Nx,info.Ny]));
    
    p_wellblock=Pressure_matrix_n1_2D(info.xw,info.yw);
    s_wellblock=Saturation_matrix_n1_2D(info.xw,info.yw);
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
    
    %         % auto-time stepping
    %         % delta_t_n=info.delta_t;
    %         delta_t_max=1;
    %         eta_s=0.05;
    %         % eta_p=50;
    %         eta_p=100;% changed from 200
    %         omega=0.5;
    %         % pressure
    % delta_p_max=max(max(abs((P_matrix_n1_k_2D_P-P_matrix_n_P))));
    delta_p_max=max(max(abs((Pressure_matrix_n1_2D-Pressure_matrix_n_2D))));
    %         delta_t_n1_p=delta_t_n*(1+omega)*eta_p/(delta_p_max+omega*eta_p);
    %         % saturation
    %         % delta_s_max=max(max(abs(P_matrix_n1_k_2D_S-P_matrix_n_S)));
    %         % delta_t_n1_s=delta_t_n*(1+omega)*eta_s/(delta_s_max+omega*eta_s);
    %         % delta_t_n1_all=[delta_t_max,delta_t_n1_p,delta_t_n1_s];
    %         % delta_t_n1=min(delta_t_n1_all);
    %         delta_t_n1_all=[delta_t_max,delta_t_n1_p];
    %         delta_t_n1=min(delta_t_n1_all);
    %         delta_t_n=delta_t_n1;
    % info.delta_t=delta_t_n1; % fix this for implicit case too
    delta_t_plot=delta_t_n;
    %delta_t_n1=delta_t_n/2;
    delta_t_n=delta_t_n/2;
    %     end
    % CFL matrix
    T_matrix_CFL_out=T_matrix_CFL(info,P_matrix_n);
    CFL_matrix=(T_matrix_CFL_out)*P_matrix_n*info.delta_t/(info.phi*info.V);
    [CFL_matrix_G,CFL_matrix_O]=to2D(info,CFL_matrix);
    CFL_matrix_out=CFL_matrix_G+CFL_matrix_O;
    CFL_max=max(max(CFL_matrix_out));
    CFL_max_vec=[CFL_max_vec,CFL_max];
    
    t_total=t_total+delta_t_plot;
    delta_t_vec=[delta_t_vec,delta_t_plot];
    t_vec=[t_vec, t_total];
    %     delta_p_max_vec=[delta_p_max_vec,delta_p_max];
    
    %     % save info to debug
    %     save("debug_plots\"+num2str(t_total)+".mat");
    %     disp(num2str(t_total));
    
    % add data to plot
    Q_matrix_all =Q_matrix_all_IMPES(info,p_wellblock,s_wellblock,control);
    [Q_matrix_2D_Qg,Q_matrix_2D_Qo]=to2D(info,Q_matrix_all);
    P_block=[P_block,Pressure_matrix_n1_2D(info.xw,info.yw)];
    P_FPR=[P_FPR,mean(mean(Pressure_matrix_n1_2D))];
    Qg_out=[Qg_out,-Q_matrix_2D_Qg(info.xw,info.yw)];
    Qo_out=[Qo_out,-Q_matrix_2D_Qo(info.xw,info.yw)];
    
    % update pressure for the next iteration
    P_matrix_n(1:2:end-1)=Pressure_matrix_n1;
    P_matrix_n(2:2:end)=Saturation_matrix_n1;
    %     Pressure_matrix_n=P_matrix_n(1:2:end-1);
    %     Saturation_matrix_n=P_matrix_n(2:2:end);
    
    % store previous pressure for ETD2 scheme
    Pressure_matrix_previous=Pressure_matrix_n;
    Saturation_matrix_previous=Saturation_matrix_n;
end
