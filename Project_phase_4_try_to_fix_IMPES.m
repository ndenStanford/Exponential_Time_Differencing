%% initialize the reservoir
tic
format long
clear all
clc

filename = 'BHP_case_1';
BHP_1 = xlsread(filename);
BHP_1_t=BHP_1(1:153,1:1);
BHP_1_p=BHP_1(1:153,2:2);

filename = 'BPR_case_1';
BPR_1 = xlsread(filename);
BPR_1_t=BPR_1(1:153,1:1);
BPR_1_p=BPR_1(1:153,2:2);

filename = 'FPR_case_1';
FPR_1 = xlsread(filename);
FPR_1_t=FPR_1(1:153,1:1);
FPR_1_p=FPR_1(1:153,2:2);

filename = 'FOPR_case_1';
FOPR_1 = xlsread(filename);
FOPR_1_t=FOPR_1(1:153,1:1);
FOPR_1_p=FOPR_1(1:153,2:2);

filename = 'FGPR_case_1';
FGPR_1 = xlsread(filename);
FGPR_1_t=FGPR_1(1:153,1:1);
FGPR_1_p=FGPR_1(1:153,2:2);

% ****************************************************************************
% unit conversion factor
info.alpha=0.001127;

% Grid and simulation parameters
info.Nx=23;
info.Ny=23;
% info.t=1;
info.delta_t=1e-3;
delta_t_n=info.delta_t;
info.t_runtime=15;%1500;

% well location
info.xw=(info.Nx+1)/2;
info.yw=(info.Ny+1)/2;

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
info.P_res_n1=info.P_res_orig;
info.S_res_n1=info.S_res_orig;

figure
imagesc(info.kx_res_orig)
set(gca,'fontsize',18)
title('X-axis permeability (md)')
colorbar

figure
imagesc(info.ky_res_orig)
set(gca,'fontsize',18)
title('Y-axis permeability (md)')
colorbar

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
Pressure_matrix_n=P_matrix_n(1:2:end-1);
Saturation_matrix_n=P_matrix_n(2:2:end);
% dP_matrix_n=zeros(2*info.Nx*info.Ny,1);
% for i=1:info.Nx*info.Ny
% dP_matrix_n(2*i-1)=-100;
% dP_matrix_n(2*i)=0.1;
% end

%% simulation

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


% change this to 'p' for the first case and 'q' for the second case
% for i=1:max_iter
control='p';

while (t_total<info.t_runtime) %% && (p_well>info.p_switch)
P_matrix_n1_k=P_matrix_n;
Pressure_matrix_n1_k=P_matrix_n1_k(1:2:end-1);
Saturation_matrix_n1_k=P_matrix_n1_k(2:2:end);
% newton iteration
% tol1=10^-5;
tol1=10^-3;
tol2=10^-2;
tol3=10^-3;

loop=true;
cond1=false;
% cond2=false;
cond3=false;
iter_newton=0;

% while  loop==true || (cond1 ==false) || (cond2 ==false) || (cond3 ==false) 
while  loop==true || (cond1 ==false)  || (cond3 ==false) 
loop=false;
% cut time step
if iter_newton>iter_newton_max
    delta_t_n=delta_t_n/2;
    iter_newton=0;
%     P_matrix_n1_k=P_matrix_n;
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
% m1=T_matrix_out*Pressure_matrix_n1_k;
% max_m1=max(abs(m1))
% m2=D_matrix_out*(Pressure_matrix_n1_k-Pressure_matrix_n);
% max_m2=max(abs(m2))
% make J matrix
J_matrix_out=J_matrix_IMPES(info,P_matrix_n1_k,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix,delta_t_n);
delta_matrix=-J_matrix_out\R_matrix;
% P_matrix_n1_k1=P_matrix_n1_k+delta_matrix;
Pressure_matrix_n1_k1=Pressure_matrix_n1_k+delta_matrix;

% [delta_matrix_2D_P,delta_matrix_2D_S]=to2D(info,delta_matrix);
% [P_matrix_n_P,P_matrix_n_S]=to2D(info,P_matrix_n);
% [P_matrix_n1_k_2D_P,P_matrix_n1_k_2D_S]=to2D(info,P_matrix_n1_k);
% [P_matrix_n1_k1_2D_P,P_matrix_n1_k1_2D_S]=to2D(info,P_matrix_n1_k1);

% calculate R normalized matrix % now what to do?
R_norm=zeros(2*info.Nx*info.Ny,1);
for ii=1:info.Ny
    for jj=1:info.Nx
        l=(ii-1)*info.Nx+jj;
%         p=P_matrix_n1_k1_2D_P(ii,jj);
%         R_norm(2*l-1,1)=abs(5.615*Bg(info,p)*delta_t_n/(info.phi*info.V)*R_matrix(2*l-1,1));
%         R_norm(2*l,1)=abs(5.615*Bo(info,p)*delta_t_n/(info.phi*info.V)*R_matrix(2*l,1));
        R_norm(l,1)=abs(5.615*delta_t_n/(info.phi*info.V)*R_matrix(1,1));
    end
end

% c1=max(abs((R_norm)))
cond1=max(abs((R_norm)))<=tol1;
% cond2=max(max(abs(P_matrix_n1_k1_2D_S-P_matrix_n1_k_2D_S)))<tol2;
% cond3=max(max(abs((P_matrix_n1_k1_2D_P-P_matrix_n1_k_2D_P)/mean(P_matrix_n1_k_2D_P))))<tol3;
% c3=max(max(abs((Pressure_matrix_n1_k1-Pressure_matrix_n1_k)/mean(Pressure_matrix_n1_k))))
cond3=max(max(abs((Pressure_matrix_n1_k1-Pressure_matrix_n1_k)/mean(Pressure_matrix_n1_k))))<tol3;

% CFL matrix
T_matrix_CFL_out=T_matrix_CFL(info,P_matrix_n1_k);
CFL_matrix=(T_matrix_CFL_out)*P_matrix_n1_k*info.delta_t/(info.phi*info.V);
[CFL_matrix_G,CFL_matrix_O]=to2D(info,CFL_matrix);
CFL_matrix_out=CFL_matrix_G+CFL_matrix_O;
CFL_max=max(max(CFL_matrix_out));

% update the P matrix for the next newton iteration
% P_matrix_n1_k=P_matrix_n1_k1;
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
% r1=T_matrix_out_Sat*Pressure_matrix_n1_k;
% r2=D_matrix_out_Sat*(Pressure_matrix_n1_k-Pressure_matrix_n);
% r3=Q_matrix_out_Sat;

CFL_max_vec=[CFL_max_vec,CFL_max];
iter_newton_vec=[iter_newton_vec,iter_newton];
delta_t_vec=[delta_t_vec,delta_t_n];


% [Q_matrix_2D_Qg,Q_matrix_2D_Qo]=to2D(info,Q_matrix_out);
% p_wellblock=P_matrix_n1_k_2D_P(info.xw,info.yw);
% s_wellblock=P_matrix_n1_k_2D_S(info.xw,info.yw);
p_wellblock=Pressure_matrix_n1_k1_2D(info.xw,info.yw);
s_wellblock=Saturation_matrix_n1_k1_2D(info.xw,info.yw);
qwo_wellblock=info.qwo_control;
% qwo_wellblock=-Q_matrix_2D_Qo(info.xw,info.yw);
To_check=info.WI*(kro(s_wellblock)/(info.mu_oil*Bo(info,p_wellblock)));
p_well=p_wellblock-(qwo_wellblock)/To_check;
% 

% % turn this section off for part 1 ****************************************
% % check if the BHP drops below 2000 for case 2
% if (p_well>info.p_switch)
%     control='q';
%     P_BHP=[P_BHP,p_well];
% else
%     control='p';
%     P_BHP=[P_BHP,info.p_switch];
% end
% % turn this section off for part 1 ****************************************
    

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
%% plots

% [P_res_2D_view,S_res_2D_view]=to2D(info,P_matrix_n1_k);
figure
% imagesc(P_res_2D_view)
imagesc(Pressure_matrix_n1_k1_2D)
% surf(P_res_2D_view)
set(gca,'fontsize',18)
title('Block Pressure Map (psi)')
colorbar
% surf(P_res_2D_view)

figure
% imagesc(S_res_2D_view)
imagesc(Saturation_matrix_n1_k1_2D)
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

% figure
% hold on
% plot(t_vec,P_BHP)
% plot(BHP_1_t,BHP_1_p)
% xlabel('time (days)')
% ylabel('pressure (psi)')
% title('bottom hole pressure')
% hold off


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
ylim([2000,8000])
xlabel('time (days)')
ylabel('Qo (STB/days)')
title('Oil flow rate: Qo')
hold off
% ****************************************************************************



