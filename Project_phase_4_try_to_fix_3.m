tic
format long
clear all
close all
clc
[BHP_1_t,BHP_1_p,BPR_1_t,BPR_1_p,FPR_1_t,FPR_1_p,FOPR_1_t,FOPR_1_p,FGPR_1_t,FGPR_1_p]=importEclipse();

% filename = 'BHP_case_1';
% BHP_1 = xlsread(filename);
% BHP_1_t=BHP_1(1:153,1:1);
% BHP_1_p=BHP_1(1:153,2:2);
% 
% filename = 'BPR_case_1';
% BPR_1 = xlsread(filename);
% BPR_1_t=BPR_1(1:153,1:1);
% BPR_1_p=BPR_1(1:153,2:2);
% 
% filename = 'FPR_case_1';
% FPR_1 = xlsread(filename);
% FPR_1_t=FPR_1(1:153,1:1);
% FPR_1_p=FPR_1(1:153,2:2);
% 
% filename = 'FOPR_case_1';
% FOPR_1 = xlsread(filename);
% FOPR_1_t=FOPR_1(1:153,1:1);
% FOPR_1_p=FOPR_1(1:153,2:2);
% 
% filename = 'FGPR_case_1';
% FGPR_1 = xlsread(filename);
% FGPR_1_t=FGPR_1(1:153,1:1);
% FGPR_1_p=FGPR_1(1:153,2:2);

% ****************************************************************************
% unit conversion factor
info.alpha=0.001127;

% Grid and simulation parameters
info.Nx=23;
info.Ny=23;
info.delta_t_max=5;
info.delta_t=0.01;
% info.delta_t=1;
info.t_runtime=10;

% sensitivity analysis
% three gird levels
% % grid
% info.Nx=45;
% info.Ny=45;
% % grid
% info.Nx=89;
% info.Ny=89;
% % grid
% info.Nx=13;
% info.Ny=13;

% well location
info.xw=(info.Nx+1)/2;
info.yw=(info.Ny+1)/2;
% info.xw=12;
% info.yw=12;

% % add additional well
% info.xwa=(info.Nx+1)/2+5;
% info.ywa=(info.Ny+1)/2+5;

% Reservoir description
info.Lx=6900;
info.Ly=6900;
info.Lz=100;
info.Dtop=1200;
info.kx=80;
info.kx_res_orig=info.kx*ones(info.Ny,info.Nx);
info.ky=120;
info.ky_res_orig=info.ky*ones(info.Ny,info.Nx);
info.phi=0.22;
info.cR=0;
info.Pinit=4500;
info.P_res_orig=info.Pinit*ones(info.Ny,info.Nx);
info.Sinit=0;
info.S_res_orig=info.Sinit*ones(info.Ny,info.Nx);
% info.P_res_n1=info.P_res_orig;
% info.S_res_n1=info.S_res_orig;


% Fluid properties
info.P_bub=3400;%*6894.76; % psi to Pa
info.P_atm=14.7;%*6894.76; % psi to Pa
info.co=(0.8e-5);%/6894.76; % psi-1 to Pa-1
info.rho_oil=49.1;
info.mu_oil=2.5;%*10^-3; % cp tp Pa*s
info.rho_gas=0.06055;

% % initialize the reservoir
info.delta_x=info.Lx/info.Nx;
info.delta_y=info.Ly/info.Ny;
info.V=info.Lz*info.delta_x*info.delta_y;

% well properties
info.rw=0.5/2;
info.s=0;
info.ro=0.28*(((info.ky/info.kx)^(1/2)*info.delta_x^2+(info.kx/info.ky)^(1/2)*info.delta_y^2)^(1/2))/((info.ky/info.kx)^(1/4)+(info.kx/info.ky)^(1/4));
info.WI=(info.alpha*2*pi*(info.kx*info.ky)^(1/2)*info.Lz/(log(info.ro/info.rw)+info.s));

% case 1 pressure control
info.p_bhp=2000;

% case 2 rate control
info.qwo_control=3000;
% info.qwo_control=10000;
info.p_switch=2000;

% make the first initial guess
% make P matrix
P_matrix_n =P_matrix(info.P_res_orig,info.S_res_orig);
P_matrix_n1_k =P_matrix_n;
delta_t_n=info.delta_t;

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
p_well=info.Pinit;
iter_newton_max=4;

% change this to 'p' for the first case and 'q' for the second case
control='p';

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
        R_matrix=(T_matrix_out*P_matrix_n1_k-D_matrix_out*(P_matrix_n1_k-P_matrix_n)+Q_matrix_out);
        
        % make J matrix
        J_matrix_out=J_matrix(info,P_matrix_n1_k,control,delta_t_n);
        % J_matrix_out=J_matrix(info,P_matrix_n1_k,control);
        
        delta_matrix=-J_matrix_out\R_matrix;
        P_matrix_n1_k1=P_matrix_n1_k+delta_matrix;
        
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
    
    % Qg_out_2=[Qg_out_2,-Q_matrix_2D_Qg(info.xwa,info.ywa)];
    % Qo_out_2=[Qo_out_2,-Q_matrix_2D_Qo(info.xwa,info.ywa)];
    
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
    delta_t_n=delta_t_n1;
    P_matrix_n=P_matrix_n1_k1;
end
toc
% now need to keep track of the total time



[P_res_2D_view,S_res_2D_view]=to2D(info,P_matrix_n1_k);
figure
imagesc(P_res_2D_view)
% surf(P_res_2D_view)
set(gca,'fontsize',18)
title('Block Pressure Map (psi)')
colorbar
% surf(P_res_2D_view)


figure
imagesc(S_res_2D_view)
set(gca,'fontsize',18)
title('Gas Saturation Map')
colorbar

figure
hold on
plot(t_vec,CFL_max_vec)
% plot(CFL_max_vec)
set(gca,'fontsize',18)
% xlabel('time (days)')
xlabel('time step')
ylabel('CFL number')
% title('CFL number VS time')
title('CFL number VS time step')
hold off

figure
[row_CFL,col_CFL]=size(CFL_matrix_out);
imagesc(CFL_matrix_out(2:row_CFL-1,2:col_CFL-1))
set(gca,'fontsize',18)
title('CFL Map')
colorbar

figure
hold on
plot(t_vec,iter_newton_vec)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('Newton iteration number')
title('Newton iteration number VS time')
hold off

figure
hold on
plot(t_vec,delta_t_vec)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('time step interval')
title('time step interval VS time')
hold off


figure
% t_vec=[0:1:max_iter]*info.delta_t;
t_vec_0=[0,t_vec];
hold on
plot(t_vec_0,P_block)
set(gca,'fontsize',18)
plot(BPR_1_t,BPR_1_p)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([2800,4000])
title('Well Block Pressure')
hold off

figure
hold on
% t_vec=[0:1:max_iter]*info.delta_t;
plot(t_vec_0,P_FPR)
plot(FPR_1_t,FPR_1_p)
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('pressure (psi)')
xlim([0,1600])
ylim([3300,4500])
legend('MATLAB output','Eclipse Output')
title('Field Pressure')
hold off

figure
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
plot(t_vec,Qg_out*5.615/10^3)
plot(FGPR_1_t(2:length(FGPR_1_t)),FGPR_1_p(2:length(FGPR_1_p)))
set(gca,'fontsize',18)
xlabel('time (days)')
ylabel('Qg (MSCF/days)')
% xlim([0,1600])
% ylim[(0.5*10^5,2.5*10^5)]
title('Gas flow rate: Qg')
hold off

figure
% t_vec=[1:1:max_iter]*info.delta_t;
hold on
plot(t_vec,Qo_out)
plot(FOPR_1_t(2:length(FOPR_1_t)),FOPR_1_p(2:length(FOPR_1_p)))
set(gca,'fontsize',18)
xlim([0,1600])
ylim([0,5000])
xlabel('time (days)')
ylabel('Qo (STB/days)')
title('Oil flow rate: Qo')
hold off
% ****************************************************************************



