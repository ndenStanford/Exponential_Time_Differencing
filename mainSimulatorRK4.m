% change this to 'p' for the first case and 'q' for the second case
% for i=1:max_iter
control='p';
% fake condition to get into the loop
delta_p_max=0;
% load L0
sum_time=0;

if select_L=="mid" || select_L=="low"
    suffix="_"+select_L;
else
    suffix="";
end
load("step0"+suffix+".mat")
% after setting L0 to a constant then can pull many terms out of
% time stepping loop
tic
delta_t_n=info.delta_t;
pinv_L=pinv(L0);
exp_Lh=expm(L0*delta_t_n);
exp_half_Lh=expm(L0*delta_t_n/2);
[row,col]=size(L0);
I=eye(row,col);
M1_1=pinv_L*(exp_Lh-I);
M1_2=exp_Lh*(I-pinv_L*L0);
M2_1=pinv_L*pinv_L;
M2_2=(1/2)*exp_Lh*(I-pinv_L*L0);
if scheme=="ETD2" || scheme=="ETDRK2"
    sum_time=sum_time+toc;
%     % do this just to time
%     tic
%     M1=M1_1+delta_t_n*M1_2;
%     M2=M2_1*(exp_Lh-(I+L0*delta_t_n))+delta_t_n^2*M2_2;
%     toc
end       


tic
phi_1=(1/delta_t_n)*(delta_t_n*exp_Lh*(I-pinv_L*L0)+pinv_L*(exp_Lh-I));
phi_1_half=(2/delta_t_n)*(delta_t_n/2*exp_half_Lh*(I-pinv_L*L0)+pinv_L*(exp_half_Lh-I));
phi_2=(1/delta_t_n^2)*(delta_t_n^2*(1/2)*exp_Lh*(I-pinv_L*L0)+pinv_L^2*(exp_Lh-(I+L0*delta_t_n)));
phi_3=(1/(2*delta_t_n^3))*(delta_t_n^3*(1/3)*exp_Lh*(I-pinv_L*L0)+pinv_L^3*(2*exp_Lh-(2*I+L0*delta_t_n*(2*I+L0*delta_t_n))));
if scheme=="ETDRK4_Lie"
    sum_time=sum_time+toc;
end

tic
Z=0*L0;
M=32;
h=delta_t_n;
r = 15*exp(1i*pi*((1:M)-.5)/M);
A = h*L0;
p_1=Z;
p_1_half=Z;
p_2=Z;
p_3=Z;
p_4=Z;
f1 = Z; f2 = Z; f3 = Z; Q = Z;

A0=L0*h;% this won't work
p_1_1=pinv(L0)*(exp_Lh-I);
%p_1_2=inv(L0)*(exp_Lh-I);
p_1_half_1=pinv(L0)*(exp_half_Lh-I);
%p_1_half_2=inv(L0)*(exp_half_Lh-I);
p_2_1=pinv(L0)^2*(exp_Lh-(I+L0*delta_t_n));
%p_2_2=inv(L0)^2*(exp_Lh-(I+L0*delta_t_n));
p_3_1=pinv(L0)^3*(2*exp_Lh-(2*I+L0*delta_t_n*(2*I+L0*delta_t_n)));
p_3_2=inv(L0)^3*(2*exp_Lh-(2*I+L0*delta_t_n*(2*I+L0*delta_t_n)));
for j = 1:M
    z = r(j);
    zIA = inv(z*I-A);
    Q = Q + h*zIA*(exp(z/2)-1);
    f1 = f1 + h*zIA*(-4-z+exp(z)*(4-3*z+z^2))/z^2;
    f2 = f2 + h*zIA*(2+z+exp(z)*(z-2))/z^2;
    f3 = f3 + h*zIA*(-4-3*z-z^2+exp(z)*(4-z))/z^2;
    % z = r(j);
    % zIA = inv(z*I-A);
    % p_1 = p_1 + h*zIA*(exp(z)-1);
    % p_1_half = p_1_half + h*zIA*(exp(z/2)-1);
    % p_2=p_2+h*zIA*(exp(z)-(1+z))/z;
    % p_3=p_3+h*zIA*(2*exp(z)-(2+z*(2+z)))/(z^2);
end
% p_1 = real(p_1/M);
% p_1_half = real(p_1_half/M);
% p_2 = h*real(p_2/M);
% p_3 = h^2*real(p_3/M);
f1 = real(f1/M); f2 = real(f2/M); f3 = real(f3/M); Q = real(Q/M);
if scheme=="ETDRK4_TF"
    sum_time=sum_time+toc;
end

% f1_Lie=delta_t_n*((4*phi_3-3*phi_2+phi_1));
% f2_Lie=delta_t_n*(-4*phi_3+2*phi_2);
% f3_Lie=delta_t_n*(4*phi_3-phi_2);

% phi_1_TF=(1/delta_t_n)*(delta_t_n*exp_Lh*(I-pinv_L*L0)+p_1);
% phi_1_half_TF=(2/delta_t_n)*(delta_t_n/2*exp_half_Lh*(I-pinv_L*L0)+p_1_half);
% phi_2_TF=(1/delta_t_n^2)*(delta_t_n^2*(1/2)*exp_Lh*(I-pinv_L*L0)+pinv_L^2*(exp_Lh-(I+L0*delta_t_n)));
% phi_3_TF=(1/(2*delta_t_n^3))*(delta_t_n^3*(1/3)*exp_Lh*(I-pinv_L*L0)+pinv_L^3*(2*exp_Lh-(2*I+L0*delta_t_n*(2*I+L0*delta_t_n))));

while (t_total<info.t_runtime) %% && (p_well>info.p_switch)
    delta_t_n=info.delta_t;
    loop=true;
    %     while (delta_p_max > 10 && t_total > 40) || loop
    loop=false;
    % P_matrix_n1_k=P_matrix_n;
    Pressure_matrix_n=P_matrix_n(1:2:end-1);
    Saturation_matrix_n=P_matrix_n(2:2:end);
    
    %L=L_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n);
    
    if (scheme=="ETDRK2")
        % do not count because can pull off the loop
        % ETDRK2 : pseudoinverse
        M1=M1_1+delta_t_n*M1_2;
        M2=M2_1*(exp_Lh-(I+L0*delta_t_n))+delta_t_n^2*M2_2;
        %an=exp_Lh*Pressure_matrix_n+M1*N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        %Pressure_matrix_n1=an+delta_t_n^-1*M2*(N_matrix_IMPES(info,an,Saturation_matrix_n,delta_t_n,control,L0)-N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        N0=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        tic
        an=exp_Lh*Pressure_matrix_n+M1*N0;
        sum_time=sum_time+toc;
        Na=N_matrix_IMPES(info,an,Saturation_matrix_n,delta_t_n,control,L0);
        tic
        Pressure_matrix_n1=an+delta_t_n^-1*M2*(Na-N0);
        sum_time=sum_time+toc;
        
    elseif (scheme=="Euler")
        % Euler
        Pressure_matrix_n1=Pressure_matrix_n+delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        
    elseif (scheme=="RK2")
        % RK2
        N0=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        %         k1=delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        tic
        k1=delta_t_n*(L0*Pressure_matrix_n+N0);
        sum_time=sum_time+toc;
        %         k2=delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n+k1,Saturation_matrix_n,delta_t_n,control,L0));
        N1=N_matrix_IMPES(info,Pressure_matrix_n+k1,Saturation_matrix_n,delta_t_n,control,L0);
        %         k2=delta_t_n*(L0*(Pressure_matrix_n+k1)+N_matrix_IMPES(info,Pressure_matrix_n+k1,Saturation_matrix_n,delta_t_n,control,L0));
        tic
        k2=delta_t_n*(L0*(Pressure_matrix_n+k1)+N1);
        %   WRONG, NEED TO RERUN????
        Pressure_matrix_n1=Pressure_matrix_n+(k1+k1)/2;
        sum_time=sum_time+toc;
        
        % Lie group method
    elseif (scheme=="ETDRK4_Lie")
        tic
        un=Pressure_matrix_n;
        u1=un;
        sum_time=sum_time+toc;
        %         N1=L0*u1+N_matrix_IMPES(info,u1,Saturation_matrix_n,delta_t_n,control,L0);
        N1=N_matrix_IMPES(info,u1,Saturation_matrix_n,delta_t_n,control,L0);
        tic
        u2=exp_half_Lh*un+(delta_t_n/2)*phi_1_half*N1;
        sum_time=sum_time+toc;
        %         N2=L0*u2+N_matrix_IMPES(info,u2,Saturation_matrix_n,delta_t_n,control,L0);
        N2=N_matrix_IMPES(info,u2,Saturation_matrix_n,delta_t_n,control,L0);
        tic
        u3=exp_half_Lh*un+(delta_t_n/2)*phi_1_half*N2;
        sum_time=sum_time+toc;
        %         N3=L0*u3+N_matrix_IMPES(info,u3,Saturation_matrix_n,delta_t_n,control,L0);
        N3=N_matrix_IMPES(info,u3,Saturation_matrix_n,delta_t_n,control,L0);
        tic
        u4=exp_half_Lh*u2+delta_t_n*phi_1_half*(-N1/2+N3);
        sum_time=sum_time+toc;
        %         N4=L0*u4+N_matrix_IMPES(info,u4,Saturation_matrix_n,delta_t_n,control,L0);
        N4=N_matrix_IMPES(info,u4,Saturation_matrix_n,delta_t_n,control,L0);
        tic
        Pressure_matrix_n1=exp_Lh*un+delta_t_n*((4*phi_3-3*phi_2+phi_1)*N1+(-4*phi_3+2*phi_2)*(N2+N3)+(4*phi_3-phi_2)*N4);
        sum_time=sum_time+toc;
        
        %         Thrieftian method
    elseif (scheme=="ETDRK4_TF")
        tic
        un=Pressure_matrix_n;
        u1=un;
        sum_time=sum_time+toc;
        %N1 = L0*u1+N_matrix_IMPES(info,u1,Saturation_matrix_n,delta_t_n,control,L0);%Nu
        N1 = N_matrix_IMPES(info,u1,Saturation_matrix_n,delta_t_n,control,L0);%Nu
        tic
        u2 = exp_half_Lh*un + Q*N1;%a
        sum_time=sum_time+toc;
        %N2 = L0*u2+N_matrix_IMPES(info,u2,Saturation_matrix_n,delta_t_n,control,L0);%Na
        N2 = N_matrix_IMPES(info,u2,Saturation_matrix_n,delta_t_n,control,L0);%Na
        tic
        u3 = exp_half_Lh*un + Q*N2;%b
        sum_time=sum_time+toc;
        %N3 = L0*u3+N_matrix_IMPES(info,u3,Saturation_matrix_n,delta_t_n,control,L0);%Nb
        N3 = N_matrix_IMPES(info,u3,Saturation_matrix_n,delta_t_n,control,L0);%Nb
        tic
        u4 = exp_half_Lh*u2 + Q*(2*N3-N1);%c
        sum_time=sum_time+toc;
        %N4 = L0*u4+N_matrix_IMPES(info,u4,Saturation_matrix_n,delta_t_n,control,L0);%Nc
        N4 = N_matrix_IMPES(info,u4,Saturation_matrix_n,delta_t_n,control,L0);%Nc
        tic
        Pressure_matrix_n1 = exp_Lh*un + f1*N1 + 2*f2*(N2+N3) + f3*N4;
        sum_time=sum_time+toc;
        %                 disp('h')
        
        % fixed
    elseif (scheme=="RK4")
        % RK2
        N0=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        %         k1=delta_t_n*(L0*Pressure_matrix_n+N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0));
        tic
        k1=delta_t_n*(L0*Pressure_matrix_n+N0);
        sum_time=sum_time+toc;
        N1=N_matrix_IMPES(info,Pressure_matrix_n+k1/2,Saturation_matrix_n,delta_t_n,control,L0);
        %         k2=delta_t_n*(L0*(Pressure_matrix_n+k1/2)+N_matrix_IMPES(info,Pressure_matrix_n+k1/2,Saturation_matrix_n,delta_t_n,control,L0));
        tic
        k2=delta_t_n*(L0*(Pressure_matrix_n+k1/2)+N1);
        sum_time=sum_time+toc;
        N2=N_matrix_IMPES(info,Pressure_matrix_n+k2/2,Saturation_matrix_n,delta_t_n,control,L0);
        %         k3=delta_t_n*(L0*(Pressure_matrix_n+k2/2)+N_matrix_IMPES(info,Pressure_matrix_n+k2/2,Saturation_matrix_n,delta_t_n,control,L0));
        tic
        k3=delta_t_n*(L0*(Pressure_matrix_n+k2/2)+N2);
        sum_time=sum_time+toc;
        N3=N_matrix_IMPES(info,Pressure_matrix_n+k3,Saturation_matrix_n,delta_t_n,control,L0);
        %         k4=delta_t_n*(L0*(Pressure_matrix_n+k3)+N_matrix_IMPES(info,Pressure_matrix_n+k3,Saturation_matrix_n,delta_t_n,control,L0));
        tic
        k4=delta_t_n*(L0*(Pressure_matrix_n+k3)+N3);
        Pressure_matrix_n1=Pressure_matrix_n+(k1/6+k2/3+k3/3+k4/6);
        sum_time=sum_time+toc;
        
    elseif (scheme=="ETD1")
        % ETD1
        M1=M1_1+delta_t_n*M1_2;
        Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
        Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+M1*Fn;
        
    elseif (scheme=="ETD2")
        % first time step need to do ETD1 because there is no Fn-1
        if t_total==0
            % ETD1
            tic
            M1=M1_1+delta_t_n*M1_2;
            sum_time=sum_time+toc;
            Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
            tic
            Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+M1*Fn;
            sum_time=sum_time+toc;
        else
            % ETD2
            % do not count because can pull off the loop
            M1=M1_1+delta_t_n*M1_2;
            M2=M2_1*(exp_Lh-(I+L0*delta_t_n))+delta_t_n^2*M2_2;
            Fn=N_matrix_IMPES(info,Pressure_matrix_n,Saturation_matrix_n,delta_t_n,control,L0);
            Fn_minus_1=N_matrix_IMPES(info,Pressure_matrix_previous,Saturation_matrix_previous,delta_t_n,control,L0);
            tic
            Pressure_matrix_n1=exp_Lh*Pressure_matrix_n+M1*Fn+M2*(Fn-Fn_minus_1)/delta_t_n;
            sum_time=sum_time+toc;
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
