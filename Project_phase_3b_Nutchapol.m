clear all
clc
% % import result from Eclipse
% filename = 'Eclipse result.xlsx';
% E = xlsread(filename);
% [m,n]=size(E);
% Pressure_table=E(10:m,2:4);
% Eclipse_time=transpose(Pressure_table(:,1)); %days
% Eclipse_pressure=transpose(Pressure_table(:,3)); %psi

filename = 'reference_J_and_R.xlsx';
E = xlsread(filename);
R_excel=E(22:33,1:1);
J_excel=E(7:18,1:12);

% unit conversion factor
alpha=0.001127;

% Reservoir description
Lx=1500;%*0.3048; %ft to m
Ly=1500;%*0.3048; %ft to m
Lz=200;%*0.3048; %ft to m
% Dtop=5000*0.3048; %ft to m
% Dbot=5100*0.3048; %ft to m
% kx=200*9.869233*10^-16; % md to m^2
kx_res_orig=[50,50,50;150,150,150];%*9.869233*10^-16; % md to m^2
% ky=100*9.869233*10^-16; % md to m^2
ky_res_orig=[200,200,200;300,300,300];%*9.869233*10^-16; % md to m^2
phi=0.22;
cR=0;
% Pinit=6000*6894.76; % psi to PA

% Fluid properties
P_bub=3500;%*6894.76; % psi to Pa
P_atm=14.7;%*6894.76; % psi to Pa
co=(0.8e-5);%/6894.76; % psi-1 to Pa-1
% rho_oil=49.1*16.018463; % lb/ft^3 to kg/m^3 not sure if it's the same lb
mu_oil=2.5;%*10^-3; % cp tp Pa*s

% Formation Volume Factor

% Grid and simulation parameters
Nx=3;
Ny=2;
% Nz=1;
% q=2000*0.1589873/(24*60*60); % bbl/day to m^3/s
% t=400*(24*60*60); % day to s
delta_t=2;%*24*60*60;
P_res_orig=[2500,2525,2550;2450,2475,2500];%*6894.76; % psi to Pa
P_res_n1=[2505,2530,2555;2455,2480,2505];%*6894.76; % psi to Pa
S_res_orig=[0,0.1,0.2;0.3,0.4,0.5];
S_res_n1=[0.1,0.2,0.3;0.4,0.5,0.6];

% % ************************************************************************
% % Test case 2
% % % increasing refinement in x
% % Nx=70;
% % % increasing refinement in y
% % Ny=70;
% % increasing refinement in y
% delta_t=0.5*24*60*60;
% % ************************************************************************
%
%
%
% % initialize the reservoir
delta_x=Lx/Nx;
delta_y=Ly/Ny;
V=Lz*delta_x*delta_y;
% P_res_orig=Pinit*ones(Ny,Nx);

% make P matrix
P_matrix_n =P_matrix(P_res_orig,S_res_orig);
P_matrix_n1 =P_matrix(P_res_n1,S_res_n1);

% make a connection list
[Px_upwinded,Py_upwinded,Sx_upwinded,Sy_upwinded,Hx_oil,Hy_oil,Hx_gas,Hy_gas]=connection_list(P_res_n1,S_res_n1);

% make D matrix
D_matrix=D_matrix(P_res_orig,S_res_orig,P_res_n1,S_res_n1);

% make T matrix
T_matrix=T_matrix(P_res_n1,S_res_n1);

% make R matrix
R_matrix=(T_matrix*P_matrix_n1-D_matrix*(P_matrix_n1-P_matrix_n));%/(0.1589873/(24*60*60)); %  m^3/s to bbl/day

% make J matrix
% contribution from D matrix
J_D=zeros(2*Nx*Ny,2*Nx*Ny);

[Ny,Nx]=size(P_res_orig);
J_D=zeros(2*Nx*Ny,2*Nx*Ny);
for i=1:Ny
    for j=1:Nx
        Pn1=P_res_n1(i,j);
        Sn1=S_res_n1(i,j);
        %         bo_small_prime=(bo_small(Pn1)-bo_small(Pn))/delta_t;
        %         bg_small_prime=(bg_small(Pn1)-bg_small(Pn))/delta_t;
        %         Rs_prime=(Rs(Pn1)-Rs(Pn))/delta_t;
        d21=(1-Sn1)*(phi*bo_small_prime(Pn1));
        d22=-(phi*bo_small(Pn1));
        d11=(1-Sn1)*(phi*bo_small(Pn1)*Rs_prime(Pn1)+phi*bo_small_prime(Pn1)*Rs(Pn1))+Sn1*(phi*bg_small_prime(Pn1));
        d12=-phi*bo_small(Pn1)*Rs(Pn1)+phi*bg_small(Pn1);
        l=(i-1)*Nx+j;
        J_D(2*l-1,2*l-1)= d11*(V/(5.615*delta_t));
        J_D(2*l-1,2*l)= d12*(V/(5.615*delta_t));
        J_D(2*l,2*l-1)= d21*(V/(5.615*delta_t));
        J_D(2*l,2*l)= d22*(V/(5.615*delta_t));
    end
end


% contribution from T matrix
J_T=zeros(2*Nx*Ny,2*Nx*Ny);

for i=1:Ny
    for j=1:Nx
        % dgamma_o/dP : i-1/2,j
        if i==1
            gamma_oil_top=0;
            dgamma_oil_dp_top=0;
            dgamma_oil_ds_top=0;
            gamma_gas_top=0;
            dgamma_gas_dp_top=0;
            dgamma_gas_ds_top=0;
            Rs_top=0;
            Rs_prime_top=0;
        else
            ky=2/(1/ky_res_orig(floor(i-1/2),j)+1/ky_res_orig(ceil(i-1/2),j));
            Ty=ky*delta_x*Lz/delta_y;
            P_upwinded=Py_upwinded(floor(i-1/2),j);
            gamma_oil_top=alpha*Ty*Hy_oil(floor(i-1/2),j);
            gamma_gas_top=alpha*Ty*Hy_gas(floor(i-1/2),j);
            Rs_top=Rs(P_upwinded);
            %             Rs_prime_top=Rs_prime(P_upwinded);
            if P_res_orig(i,j)<=P_res_orig(i-1,j)
                dgamma_oil_dp_top=0;
                dgamma_oil_ds_top=0;
                dgamma_oil_dp_top_top=alpha*Ty*kro(S_res_n1(i-1,j))*(bo_small_prime(P_res_n1(i-1,j))/mu_oil);
                dgamma_oil_ds_top_top=alpha*Ty*kro_prime(S_res_n1(i-1,j))/(Bo(P_res_n1(i-1,j))*mu_oil);
                dgamma_gas_dp_top=0;
                dgamma_gas_ds_top=0;
                dgamma_gas_dp_top_top=alpha*Ty*krg(S_res_n1(i-1,j))*(bg_small_prime(P_res_n1(i-1,j))/(mu_g(P_res_n1(i-1,j)))-(mu_g_prime(P_res_n1(i-1,j))*bg_small(P_res_n1(i-1,j))/(mu_g(P_res_n1(i-1,j)))^2));
                dgamma_gas_ds_top_top=alpha*Ty*krg_prime(S_res_n1(i-1,j))/(Bg(P_res_n1(i-1,j))*mu_g(P_res_n1(i-1,j)));
                Rs_prime_top=0;
                Rs_prime_top_top=Rs_prime(P_res_n1(i-1,j));
                mu_g_prime_top=0;
                mu_g_prime_top_top=mu_g_prime(P_res_n1(i-1,j));
                
            else
                dgamma_oil_dp_top=alpha*Ty*kro(S_res_n1(i,j))*(bo_small_prime(P_res_n1(i,j))/mu_oil);
                dgamma_oil_ds_top=alpha*Ty*kro_prime(S_res_n1(i,j))/(Bo(P_res_n1(i,j))*mu_oil);
                dgamma_oil_dp_top_top=0;
                dgamma_oil_ds_top_top=0;
                dgamma_gas_dp_top=alpha*Ty*krg(S_res_n1(i,j))*(bg_small_prime(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))-(mu_g_prime(P_res_n1(i,j))*bg_small(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))^2));
                dgamma_gas_ds_top=alpha*Ty*krg_prime(S_res_n1(i,j))/(Bg(P_res_n1(i,j))*mu_g(P_res_n1(i,j)));
                dgamma_gas_dp_top_top=0;
                dgamma_gas_ds_top_top=0;
                Rs_prime_top=Rs_prime(P_res_n1(i,j));
                Rs_prime_top_top=0;
                mu_g_prime_top=mu_g_prime(P_res_n1(i,j));
                mu_g_prime_top_top=0;
            end
        end
        % dgamma_o/dP : i+1/2,j
        if i==Ny
            gamma_oil_bottom=0;
            dgamma_oil_dp_bottom=0;
            dgamma_oil_ds_bottom=0;
            gamma_gas_bottom=0;
            dgamma_gas_dp_bottom=0;
            dgamma_gas_ds_bottom=0;
            Rs_bottom=0;
            Rs_prime_bottom=0;
        else
            ky=2/(1/ky_res_orig(floor(i+1/2),j)+1/ky_res_orig(ceil(i+1/2),j));
            Ty=ky*delta_x*Lz/delta_y;
            P_upwinded=Py_upwinded(floor(i+1/2),j);
            gamma_oil_bottom=alpha*Ty*Hy_oil(floor(i+1/2),j);
            gamma_gas_bottom=alpha*Ty*Hy_gas(floor(i+1/2),j);
            Rs_bottom=Rs(P_upwinded);
            %             Rs_prime_bottom=Rs_prime(P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i+1,j)
                dgamma_oil_dp_bottom=0;
                dgamma_oil_ds_bottom=0;
                dgamma_oil_dp_bottom_bottom=alpha*Ty*kro(S_res_n1(i+1,j))*(bo_small_prime(P_res_n1(i+1,j))/mu_oil);
                dgamma_oil_ds_bottom_bottom=alpha*Ty*kro_prime(S_res_n1(i+1,j))/(Bo(P_res_n1(i+1,j))*mu_oil);
                dgamma_gas_dp_bottom=0;
                dgamma_gas_ds_bottom=0;
                dgamma_gas_dp_bottom_bottom=alpha*Ty*krg(S_res_n1(i+1,j))*(bg_small_prime(P_res_n1(i+1,j))/(mu_g(P_res_n1(i+1,j)))-(mu_g_prime(P_res_n1(i+1,j))*bg_small(P_res_n1(i+1,j))/(mu_g(P_res_n1(i+1,j)))^2));
                dgamma_gas_ds_bottom_bottom=alpha*Ty*kro_prime(S_res_n1(i+1,j))/(Bg(P_res_n1(i+1,j))*mu_g(P_res_n1(i+1,j)));
                Rs_prime_bottom=0;
                Rs_prime_bottom_bottom=Rs_prime(P_res_n1(i+1,j));
                mu_g_prime_bottom=0;
                mu_g_prime_bottom_bottom=mu_g_prime(P_res_n1(i+1,j));
                
            else
                dgamma_oil_dp_bottom=alpha*Ty*kro(S_res_n1(i,j))*(bo_small_prime(P_res_n1(i,j))/mu_oil);
                dgamma_oil_ds_bottom=alpha*Ty*kro_prime(S_res_n1(i,j))/(Bo(P_res_n1(i,j))*mu_oil);
                dgamma_oil_dp_bottom_bottom=0;
                dgamma_oil_ds_bottom_bottom=0;
                dgamma_gas_dp_bottom=alpha*Ty*krg(S_res_n1(i,j))*(bg_small_prime(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))-(mu_g_prime(P_res_n1(i,j))*bg_small(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))^2));
                dgamma_gas_ds_bottom=alpha*Ty*krg_prime(S_res_n1(i,j))/(Bg(P_res_n1(i,j))*mu_g(P_res_n1(i,j)));
                dgamma_gas_dp_bottom_bottom=0;
                dgamma_gas_ds_bottom_bottom=0;
                Rs_prime_bottom=Rs_prime(P_res_n1(i,j));
                Rs_prime_bottom_bottom=0;
                mu_g_prime_bottom=mu_g_prime(P_res_n1(i,j));
                mu_g_prime_bottom_bottom=0;
            end
        end
        % dgamma_o/dP : i,j-1/2
        if j==1
            gamma_oil_left=0;
            dgamma_oil_dp_left=0;
            dgamma_oil_ds_left=0;
            gamma_gas_left=0;
            dgamma_gas_dp_left=0;
            dgamma_gas_ds_left=0;
            Rs_left=0;
            Rs_prime_left=0;
        else
            kx=(kx_res_orig(i,floor(j-1/2))+kx_res_orig(i,ceil(j-1/2)))/2;
            Tx=kx*delta_y*Lz/delta_x;
            P_upwinded=Px_upwinded(i,floor(j-1/2));
            gamma_oil_left=alpha*Tx*Hx_oil(i,floor(j-1/2));
            gamma_gas_left=alpha*Tx*Hx_gas(i,floor(j-1/2));
            Rs_left=Rs(P_upwinded);
            %             Rs_prime_left=Rs_prime(P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i,j-1)
                dgamma_oil_dp_left=0;
                dgamma_oil_ds_left=0;
                dgamma_oil_dp_left_left=alpha*Tx*kro(S_res_n1(i,j-1))*(bo_small_prime(P_res_n1(i,j-1))/mu_oil);
                dgamma_oil_ds_left_left=alpha*Tx*kro_prime(S_res_n1(i,j-1))/(Bo(P_res_n1(i,j-1))*mu_oil);
                dgamma_gas_dp_left=0;
                dgamma_gas_ds_left=0;
                dgamma_gas_dp_left_left=alpha*Tx*krg(S_res_n1(i,j-1))*(bg_small_prime(P_res_n1(i,j-1))/(mu_g(P_res_n1(i,j-1)))-(mu_g_prime(P_res_n1(i,j-1))*bg_small(P_res_n1(i,j-1))/(mu_g(P_res_n1(i,j-1)))^2));
                dgamma_gas_ds_left_left=alpha*Tx*krg_prime(S_res_n1(i,j-1))/(Bg(P_res_n1(i,j-1))*mu_g(P_res_n1(i,j-1)));
                Rs_prime_left=0;
                Rs_prime_left_left=Rs_prime(P_res_n1(i,j-1));
                mu_g_prime_left=0;
                mu_g_prime_left_left=mu_g_prime(P_res_n1(i,j-1));
            else
                dgamma_oil_dp_left=alpha*Tx*kro(S_res_n1(i,j))*(bo_small_prime(P_res_n1(i,j))/mu_oil);
                dgamma_oil_ds_left=alpha*Tx*kro_prime(S_res_n1(i,j))/(Bo(P_res_n1(i,j))*mu_oil);
                dgamma_oil_dp_left_left=0;
                dgamma_oil_ds_left_left=0;
                dgamma_gas_dp_left=alpha*Tx*krg(S_res_n1(i,j))*(bg_small_prime(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))-(mu_g_prime(P_res_n1(i,j))*bg_small(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))^2));
                dgamma_gas_ds_left=alpha*Tx*krg_prime(S_res_n1(i,j))/(Bg(P_res_n1(i,j))*mu_g(P_res_n1(i,j)));
                dgamma_gas_dp_left_left=0;
                dgamma_gas_ds_left_left=0;
                Rs_prime_left=Rs_prime(P_res_n1(i,j));
                Rs_prime_left_left=0;
                mu_g_prime_left=mu_g_prime(P_res_n1(i,j));
                mu_g_prime_left_left=0;
            end
        end
        % dgamma_o/dP : i,j+1/2
        if j==Nx
            gamma_oil_right=0;
            dgamma_oil_dp_right=0;
            dgamma_oil_ds_right=0;
            gamma_gas_right=0;
            dgamma_gas_dp_right=0;
            dgamma_gas_ds_right=0;
            Rs_right=0;
            Rs_prime_right=0;
        else
            kx=(kx_res_orig(i,floor(j+1/2))+kx_res_orig(i,ceil(j+1/2)))/2;
            Tx=kx*delta_y*Lz/delta_x;
            P_upwinded=Px_upwinded(i,floor(j+1/2));
            gamma_oil_right=alpha*Tx*Hx_oil(i,floor(j+1/2));
            gamma_gas_right=alpha*Tx*Hx_gas(i,floor(j+1/2));
            Rs_right=Rs(P_upwinded);
            %             Rs_prime_right=Rs_prime(P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i,j+1)
                dgamma_oil_dp_right=0;
                dgamma_oil_ds_right=0;
                dgamma_oil_dp_right_right=alpha*Tx*kro(S_res_n1(i,j+1))*(bo_small_prime(P_res_n1(i,j+1))/mu_oil);
                dgamma_oil_ds_right_right=alpha*Tx*kro_prime(S_res_n1(i,j+1))/(Bo(P_res_n1(i,j+1))*mu_oil);
                dgamma_gas_dp_right=0;
                dgamma_gas_ds_right=0;
                dgamma_gas_dp_right_right=alpha*Tx*krg(S_res_n1(i,j+1))*(bg_small_prime(P_res_n1(i,j+1))/(mu_g(P_res_n1(i,j+1)))-(mu_g_prime(P_res_n1(i,j+1))*bg_small(P_res_n1(i,j+1))/(mu_g(P_res_n1(i,j+1)))^2));
                dgamma_gas_ds_right_right=alpha*Tx*krg_prime(S_res_n1(i,j+1))/(Bg(P_res_n1(i,j+1))*mu_g(P_res_n1(i,j+1)));
                Rs_prime_right=0;
                Rs_prime_right_right=Rs_prime(P_res_n1(i,j+1));
                mu_g_prime_right=0;
                mu_g_prime_right_right=mu_g_prime(P_res_n1(i,j));
            else
                dgamma_oil_dp_right=alpha*Tx*kro(S_res_n1(i,j))*(bo_small_prime(P_res_n1(i,j))/mu_oil);
                dgamma_oil_ds_right=alpha*Tx*kro_prime(S_res_n1(i,j))/(Bo(P_res_n1(i,j))*mu_oil);
                dgamma_oil_dp_right_right=0;
                dgamma_oil_ds_right_right=0;
                dgamma_gas_dp_right=alpha*Tx*krg(S_res_n1(i,j))*(bg_small_prime(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))-(mu_g_prime(P_res_n1(i,j))*bg_small(P_res_n1(i,j))/(mu_g(P_res_n1(i,j)))^2));
                dgamma_gas_ds_right=alpha*Tx*krg_prime(S_res_n1(i,j))/(Bg(P_res_n1(i,j))*mu_g(P_res_n1(i,j)));
                dgamma_gas_dp_right_right=0;
                dgamma_gas_ds_right_right=0;
                Rs_prime_right=Rs_prime(P_res_n1(i,j));
                Rs_prime_right_right=0;
                mu_g_prime_right=mu_g_prime(P_res_n1(i,j));
                mu_g_prime_right_right=0;
            end
        end
        
        P_mid=P_res_n1(i,j);
        if i-1>=1
            P_top=P_res_n1(i-1,j);
        else
            P_top=0;
        end
        if i+1<=Ny
            P_bottom=P_res_n1(i+1,j);
        else
            P_bottom=0;
        end
        if j-1>=1
            P_left=P_res_n1(i,j-1);
        else
            P_left=0;
        end
        if j+1<=Nx
            P_right=P_res_n1(i,j+1);
        else
            P_right=0;
        end
        sum_gamma_oil=gamma_oil_left+gamma_oil_right+gamma_oil_top+gamma_oil_bottom;
        sum_dgamma_oil_dp=dgamma_oil_dp_left+dgamma_oil_dp_right+dgamma_oil_dp_top+dgamma_oil_dp_bottom;
        sum_dgamma_oil_ds=dgamma_oil_ds_left+dgamma_oil_ds_right+dgamma_oil_ds_top+dgamma_oil_ds_bottom;
        dR_dP_oil=(dgamma_oil_dp_left*(P_left-P_mid)+dgamma_oil_dp_right*(P_right-P_mid)+dgamma_oil_dp_top*(P_top-P_mid)+dgamma_oil_dp_bottom*(P_bottom-P_mid)-sum_gamma_oil);
        dR_dS_oil=(dgamma_oil_ds_left*(P_left-P_mid)+dgamma_oil_ds_right*(P_right-P_mid)+dgamma_oil_ds_top*(P_top-P_mid)+dgamma_oil_ds_bottom*(P_bottom-P_mid));
        
        sum_gamma_gas=gamma_gas_left+gamma_gas_right+gamma_gas_top+gamma_gas_bottom;
        sum_dgamma_gas_dp=dgamma_gas_dp_left+dgamma_gas_dp_right+dgamma_gas_dp_top+dgamma_gas_dp_bottom;
        sum_dgamma_gas_ds=dgamma_gas_ds_left+dgamma_gas_ds_right+dgamma_gas_ds_top+dgamma_gas_ds_bottom;
        
        dR_dP_gas=((Rs_left*dgamma_oil_dp_left+gamma_oil_left*Rs_prime_left+dgamma_gas_dp_left)*(P_left-P_mid)+(Rs_right*dgamma_oil_dp_right+gamma_oil_right*Rs_prime_right+dgamma_gas_dp_right)*(P_right-P_mid)+(Rs_top*dgamma_oil_dp_top+gamma_oil_top*Rs_prime_top+dgamma_gas_dp_top)*(P_top-P_mid)+(Rs_bottom*dgamma_oil_dp_bottom+gamma_oil_bottom*Rs_prime_bottom+dgamma_gas_dp_bottom)*(P_bottom-P_mid)-((Rs_left*gamma_oil_left+gamma_gas_left)+(Rs_right*gamma_oil_right+gamma_gas_right)+(Rs_top*gamma_oil_top+gamma_gas_top)+(Rs_bottom*gamma_oil_bottom+gamma_gas_bottom)));
        dR_dS_gas=((Rs_left*dgamma_oil_ds_left+dgamma_gas_ds_left)*(P_left-P_mid)+(Rs_right*dgamma_oil_ds_right+dgamma_gas_ds_right)*(P_right-P_mid)+(Rs_top*dgamma_oil_ds_top+dgamma_gas_ds_top)*(P_top-P_mid)+(Rs_bottom*dgamma_oil_ds_bottom+dgamma_gas_ds_bottom)*(P_bottom-P_mid));
        
        l=(i-1)*Nx+j;
        J_T(2*l,2*l-1)= dR_dP_oil;
        J_T(2*l,2*l)= dR_dS_oil;
        J_T(2*l-1,2*l-1)= dR_dP_gas;
        J_T(2*l-1,2*l)= dR_dS_gas;
        
%         if l-1>=1
       if j>1
            %             dR_dP_oil_left=(1/V)*(dgamma_oil_dp_left*(P_left-P_mid)+gamma_oil_left);
            %             dR_dS_oil_left=(1/V)*dgamma_oil_ds_left*(P_left-P_mid)+gamma_oil_left;
            dR_dP_oil_left=(dgamma_oil_dp_left_left*(P_left-P_mid)+gamma_oil_left);
            dR_dS_oil_left=dgamma_oil_ds_left_left*(P_left-P_mid);
            J_T(2*l,2*l-3)= dR_dP_oil_left;
            J_T(2*l,2*l-2)= dR_dS_oil_left;
            %             dR_dP_gas_left=((Rs_left*gamma_oil_left+gamma_gas_left)+(Rs_left*dgamma_oil_dp_left_left+gamma_oil_left*Rs_prime_left+dgamma_gas_dp_left_left)*(P_left-P_mid));
            %             dR_dS_gas_left=((Rs_prime_left*gamma_oil_left+Rs_left*dgamma_oil_ds_left_left+dgamma_gas_ds_left_left)*(P_left-P_mid));
            %             dR_dP_gas_left=((Rs_left*gamma_oil_left+gamma_gas_left)+(Rs_left*dgamma_oil_dp_left_left+gamma_oil_left*Rs_prime_left_left+dgamma_gas_dp_left_left)*(P_left-P_mid));
            %             dR_dS_gas_left=((Rs_prime_left_left*gamma_oil_left+Rs_left*dgamma_oil_ds_left_left+dgamma_gas_ds_left_left)*(P_left-P_mid));
            dR_dP_gas_left=((Rs_left*gamma_oil_left+gamma_gas_left)+(Rs_left*dgamma_oil_dp_left_left+gamma_oil_left*Rs_prime_left_left+dgamma_gas_dp_left_left)*(P_left-P_mid));
            dR_dS_gas_left=((Rs_left*dgamma_oil_ds_left_left+dgamma_gas_ds_left_left)*(P_left-P_mid));
            J_T(2*l-1,2*l-3)= dR_dP_gas_left;
            J_T(2*l-1,2*l-2)= dR_dS_gas_left;
        end
%         if l+1<=Nx*Ny
        if j<Nx
            %             dR_dP_oil_right=(1/V)*(dgamma_oil_dp_right*(P_right-P_mid)+gamma_oil_right);
            %             dR_dS_oil_right=(1/V)*dgamma_oil_ds_right*(P_right-P_mid)+gamma_oil_right;
            dR_dP_oil_right=(dgamma_oil_dp_right_right*(P_right-P_mid)+gamma_oil_right);
            dR_dS_oil_right=dgamma_oil_ds_right_right*(P_right-P_mid);
            J_T(2*l,2*l+1)= dR_dP_oil_right;
            J_T(2*l,2*l+2)= dR_dS_oil_right;
            %             dR_dP_gas_right=((Rs_right*gamma_oil_right+gamma_gas_right)+(Rs_right*dgamma_oil_dp_right_right+gamma_oil_right*Rs_prime_right+dgamma_gas_dp_right_right)*(P_right-P_mid));
            %             dR_dS_gas_right=((Rs_prime_right*gamma_oil_right+Rs_right*dgamma_oil_ds_right_right+dgamma_gas_ds_right_right)*(P_right-P_mid));
            %             dR_dP_gas_right=((Rs_right*gamma_oil_right+gamma_gas_right)+(Rs_right*dgamma_oil_dp_right_right+gamma_oil_right*Rs_prime_right_right+dgamma_gas_dp_right_right)*(P_right-P_mid));
            %             dR_dS_gas_right=((Rs_prime_right_right*gamma_oil_right+Rs_right*dgamma_oil_ds_right_right+dgamma_gas_ds_right_right)*(P_right-P_mid));
            dR_dP_gas_right=((Rs_right*gamma_oil_right+gamma_gas_right)+(Rs_right*dgamma_oil_dp_right_right+gamma_oil_right*Rs_prime_right_right+dgamma_gas_dp_right_right)*(P_right-P_mid));
            dR_dS_gas_right=((Rs_right*dgamma_oil_ds_right_right+dgamma_gas_ds_right_right)*(P_right-P_mid));
            J_T(2*l-1,2*l+1)= dR_dP_gas_right;
            J_T(2*l-1,2*l+2)= dR_dS_gas_right;
        end
%         if l-Nx>=1
        if i>1
            %             dR_dP_oil_top=(1/V)*(dgamma_oil_dp_top*(P_top-P_mid)+gamma_oil_top);
            %             dR_dS_oil_top=(1/V)*dgamma_oil_ds_top*(P_top-P_mid)+gamma_oil_top;
            dR_dP_oil_top=(dgamma_oil_dp_top_top*(P_top-P_mid)+gamma_oil_top);
            dR_dS_oil_top=dgamma_oil_ds_top_top*(P_top-P_mid);
            J_T(2*l,2*l-1-2*Nx)= dR_dP_oil_top;
            J_T(2*l,2*l-2*Nx)= dR_dS_oil_top;
            %             dR_dP_gas_top=((Rs_top*gamma_oil_top+gamma_gas_top)+(Rs_top*dgamma_oil_dp_top_top+gamma_oil_top*Rs_prime_top+dgamma_gas_dp_top_top)*(P_top-P_mid));
            %             dR_dS_gas_top=((Rs_prime_top*gamma_oil_top+Rs_top*dgamma_oil_ds_top_top+dgamma_gas_ds_top_top)*(P_top-P_mid));
            %             dR_dP_gas_top=((Rs_top*gamma_oil_top+gamma_gas_top)+(Rs_top*dgamma_oil_dp_top_top+gamma_oil_top*Rs_prime_top_top+dgamma_gas_dp_top_top)*(P_top-P_mid));
            %             dR_dS_gas_top=((Rs_prime_top_top*gamma_oil_top+Rs_top*dgamma_oil_ds_top_top+dgamma_gas_ds_top_top)*(P_top-P_mid));
            dR_dP_gas_top=((Rs_top*gamma_oil_top+gamma_gas_top)+(Rs_top*dgamma_oil_dp_top_top+gamma_oil_top*Rs_prime_top_top+dgamma_gas_dp_top_top)*(P_top-P_mid));
            dR_dS_gas_top=((Rs_top*dgamma_oil_ds_top_top+dgamma_gas_ds_top_top)*(P_top-P_mid));
            J_T(2*l-1,2*l-1-2*Nx)= dR_dP_gas_top;
            J_T(2*l-1,2*l-2*Nx)= dR_dS_gas_top;
        end
%         if l+Nx<=Nx*Ny
        if i<Ny
            %             dR_dP_oil_bottom=(1/V)*(dgamma_oil_dp_bottom*(P_bottom-P_mid)+gamma_oil_bottom);
            %             dR_dS_oil_bottom=(1/V)*dgamma_oil_ds_bottom*(P_bottom-P_mid)+gamma_oil_bottom;
            dR_dP_oil_bottom=(dgamma_oil_dp_bottom_bottom*(P_bottom-P_mid)+gamma_oil_bottom);
            dR_dS_oil_bottom=dgamma_oil_ds_bottom_bottom*(P_bottom-P_mid);
            J_T(2*l,2*l-1+2*Nx)= dR_dP_oil_bottom;
            J_T(2*l,2*l+2*Nx)= dR_dS_oil_bottom;
            %             dR_dP_gas_bottom=((Rs_bottom*gamma_oil_bottom+gamma_gas_bottom)+(Rs_bottom*dgamma_oil_dp_bottom_bottom+gamma_oil_bottom*Rs_prime_bottom+dgamma_gas_dp_bottom_bottom)*(P_bottom-P_mid));
            %             dR_dS_gas_bottom=((Rs_prime_bottom*gamma_oil_bottom+Rs_bottom*dgamma_oil_ds_bottom_bottom+dgamma_gas_ds_bottom_bottom)*(P_bottom-P_mid));
            %             dR_dP_gas_bottom=((Rs_bottom*gamma_oil_bottom+gamma_gas_bottom)+(Rs_bottom*dgamma_oil_dp_bottom_bottom+gamma_oil_bottom*Rs_prime_bottom_bottom+dgamma_gas_dp_bottom_bottom)*(P_bottom-P_mid));
            %             dR_dS_gas_bottom=((Rs_prime_bottom_bottom*gamma_oil_bottom+Rs_bottom*dgamma_oil_ds_bottom_bottom+dgamma_gas_ds_bottom_bottom)*(P_bottom-P_mid));
            dR_dP_gas_bottom=((Rs_bottom*gamma_oil_bottom+gamma_gas_bottom)+(Rs_bottom*dgamma_oil_dp_bottom_bottom+gamma_oil_bottom*Rs_prime_bottom_bottom+dgamma_gas_dp_bottom_bottom)*(P_bottom-P_mid));
            dR_dS_gas_bottom=((Rs_bottom*dgamma_oil_ds_bottom_bottom+dgamma_gas_ds_bottom_bottom)*(P_bottom-P_mid));
            J_T(2*l-1,2*l-1+2*Nx)= dR_dP_gas_bottom;
            J_T(2*l-1,2*l+2*Nx)= dR_dS_gas_bottom;
        end
        
        
    end
end
% spy(J_T)
J_matrix=J_T-J_D;

% compare the result with the reference case
dR_matrix=abs((R_matrix-R_excel)./R_excel);
[row,col]=size(J_matrix);
dJ_matrix=zeros(row,col);
for i=1:row
    for j=1:col
        if J_matrix(i,j)~=0
            dJ_matrix(i,j)=abs((J_matrix(i,j)-J_excel(i,j))/J_excel(i,j));
        end
    end
end


% % figure
% % imagesc(Pn_field/6894.76)
% % title('Pressure map for the entire field at the end of the simulation run, t = 400 days');
% % colorbar



% %
% % for i=normal_block
% % byu=b((Pn(i-1)+Pn(i))/2); % b y-axis up
% % byd=b((Pn(i)+Pn(i+1))/2); % b y-axis down
% % bxl=b((Pn(i-Ny)+Pn(i))/2); % b x-axis left
% % bxr=b((Pn(i)+Pn(i+Ny))/2); % b x-axis right
% % T(i,i)=-(gamma_y*byu+gamma_y*byd+gamma_x*bxl+gamma_x*bxr);
% % T(i,i-1)=gamma_y*byu;
% % T(i,i+1)=gamma_y*byd;
% % T(i,i-Ny)=gamma_x*bxl;
% % T(i,i+Ny)=gamma_x*bxr;
% % end
% %
% %
% % % % boundary condition, modify all edged block
% % % modify each row one by one
% % % left edge (not including corners)
% % for i=setdiff(left_edge,corner)
% % byu=b((Pn(i-1)+Pn(i))/2); % b y-axis up
% % byd=b((Pn(i)+Pn(i+1))/2); % b y-axis down
% % bxr=b((Pn(i)+Pn(i+Ny))/2); % b x-axis right
% % T(i,i)=-(gamma_y*byu+gamma_y*byd+gamma_x*bxr);
% % T(i,i-1)=gamma_y*byu;
% % T(i,i+1)=gamma_y*byd;
% % T(i,i+Ny)=gamma_x*bxr;
% % end
% %
% % % right edge (not including corners)
% % for i=setdiff(right_edge,corner)
% % byu=b((Pn(i-1)+Pn(i))/2); % b y-axis up
% % byd=b((Pn(i)+Pn(i+1))/2); % b y-axis down
% % bxl=b((Pn(i-Ny)+Pn(i))/2); % b x-axis left
% % T(i,i)=-(gamma_y*byu+gamma_y*byd+gamma_x*bxl);
% % T(i,i-1)=gamma_y*byu;
% % T(i,i+1)=gamma_y*byd;
% % T(i,i-Ny)=gamma_x*bxl;
% % end
% %
% % % top edge (not including corners)
% % for i=setdiff(top_edge,corner)
% % byd=b((Pn(i)+Pn(i+1))/2); % b y-axis down
% % bxl=b((Pn(i-Ny)+Pn(i))/2); % b x-axis left
% % bxr=b((Pn(i)+Pn(i+Ny))/2); % b x-axis right
% % T(i,i)=-(gamma_y*byd+gamma_x*bxl+gamma_x*bxr);
% % T(i,i+1)=gamma_y*byd;
% % T(i,i-Ny)=gamma_x*bxl;
% % T(i,i+Ny)=gamma_x*bxr;
% % end
% %
% % % bottom edge (not including corners)
% % for i=setdiff(bottom_edge,corner)
% % byu=b((Pn(i-1)+Pn(i))/2); % b y-axis up
% % bxl=b((Pn(i-Ny)+Pn(i))/2); % b x-axis left
% % bxr=b((Pn(i)+Pn(i+Ny))/2); % b x-axis right
% % T(i,i)=-(gamma_y*byu+gamma_x*bxl+gamma_x*bxr);
% % T(i,i-1)=gamma_y*byu;
% % T(i,i-Ny)=gamma_x*bxl;
% % T(i,i+Ny)=gamma_x*bxr;
% % end
% %
% % % top left corner
% % for i=1
% % byd=b((Pn(i)+Pn(i+1))/2); % b y-axis down
% % bxr=b((Pn(i)+Pn(i+Ny))/2); % b x-axis right
% % T(i,i)=-(gamma_y*byd+gamma_x*bxr);
% % T(i,i+1)=gamma_y*byd;
% % T(i,i+Ny)=gamma_x*bxr;
% % end
% %
% % % bottom left corner
% % for i=Ny
% % byu=b((Pn(i-1)+Pn(i))/2); % b y-axis up
% % bxr=b((Pn(i)+Pn(i+Ny))/2); % b x-axis right
% % T(i,i)=-(gamma_y*byu+gamma_x*bxr);
% % T(i,i-1)=gamma_y*byu;
% % T(i,i+Ny)=gamma_x*bxr;
% % end
% %
% % % top right corner
% % for i=Nx*Ny-Ny+1
% % byd=b((Pn(i)+Pn(i+1))/2); % b y-axis down
% % bxl=b((Pn(i-Ny)+Pn(i))/2); % b x-axis left
% % T(i,i)=-(gamma_y*byd+gamma_x*bxl);
% % T(i,i+1)=gamma_y*byd;
% % T(i,i-Ny)=gamma_x*bxl;
% % end
% %
% % % bottom right corner
% % for i=Nx*Ny
% % byu=b((Pn(i-1)+Pn(i))/2); % b y-axis up
% % bxl=b((Pn(i-Ny)+Pn(i))/2); % b x-axis left
% % T(i,i)=-(gamma_y*byu+gamma_x*bxl);
% % T(i,i-1)=gamma_y*byu;
% % T(i,i-Ny)=gamma_x*bxl;
% % end
% %
% % % buid Q matrix
% % Q_res_orig=zeros(Ny,Nx);
% % % input a well, comment tis section out when do some test case 2
% % Q_res_orig((Ny+1)/2,(Nx+1)/2)=q;
% %
% % % % ************************************************************************
% % % % Test case 1
% % % Q_res_orig((Ny+1)/2,(Nx+1)/2+5)=q;
% % % Q_res_orig((Ny+1)/2,(Nx+1)/2-5)=q;
% % % % ************************************************************************
% %
% % % ************************************************************************
% % % Test case 2
% % % % refine x
% % % Q_res_orig((Ny+1)/2,Nx/2)=q/2;
% % % Q_res_orig((Ny+1)/2,Nx/2+1)=q/2;
% % % % refine y
% % % Q_res_orig(Ny/2,(Nx+1)/2)=q/2;
% % % Q_res_orig(Ny/2+1,(Nx+1)/2)=q/2;
% % % ************************************************************************
% %
% % Q_res=reshape(Q_res_orig,[],1);
% %
% % % buid B matrix
% % B0=1;
% % Ci=(cf*phi/B0);
% % Vi=delta_x*delta_y*Lz;
% % B_res=diag((Ci*Vi)*ones(1,Nx*Ny));
% %
% %
% % A=T-(1/delta_t)*B_res;
% %
% %
% % b_matrix=Q_res-(1/delta_t)*(B_res*Pn);
% % Pn1=A\b_matrix;
% % Pn=Pn1;
% % Pn_field=reshape(Pn,Ny,Nx);
% % P_field_all(:,:,index)=Pn_field;
% % end
% %
% % P_well_all=squeeze(P_field_all((Ny+1)/2,(Nx+1)/2,:));
% % % *************************************************************
% % % test case 2
% % % % refine x
% % % P_well_all=squeeze(((P_field_all((Ny+1)/2,Nx/2,:))+(P_field_all((Ny+1)/2,Nx/2+1,:)))/2);
% % % % refine y
% % % P_well_all=squeeze(((P_field_all(Ny/2,(Nx+1)/2,:))+(P_field_all(Ny/2+1,(Nx+1)/2+1,:)))/2);
% % % *************************************************************
% %
% %
% % figure
% % hold on
% % plot(t_vec/(24*60*60),P_well_all/6894.76);% convert s to day and Pa to psi
% % plot(Eclipse_time,Eclipse_pressure);
% % legend('MATLAB result','Eclipse result');
% % title('Pressure drop over time');
% % xlabel('time (day)');
% % ylabel('pressure (psi)');
% % hold off
% %
% % % % *************************************************************************
% % % % save matrix for plot later, change file name when save
% % % t_save=t_vec/(24*60*60);
% % % p_save=P_well_all/6894.76;
% % % save('pVSt_refine_t','t_save','p_save')
% % % % *************************************************************************
% %
% % % % *************************************************************************
% % % % % load matrix to plot: test case 1
% % % load = matfile('pVSt_one_well.mat');
% % % t_save=load.t_save;
% % % p_save=load.p_save;
% % %
% % % % load_3 = matfile('pVSt_three_well.mat');
% % % % t_save_3=load_3.t_save;
% % % % p_save_3=load_3.p_save;
% % %
% % % figure
% % % hold on
% % % plot(t_save,p_save);
% % % plot(t_save_3,p_save_3);
% % % legend('Base case (one well)','Test case (three wells)');
% % % title('Pressure drop over time');
% % % xlabel('time (day)');
% % % ylabel('pressure (psi)');
% % % hold off
% % % % % *************************************************************************
% %
% %
% % % % *************************************************************************
% % % % load matrix to plot: test case 2
% % %
% % % load = matfile('pVSt_one_well.mat');
% % % t_save=load.t_save;
% % % p_save=load.p_save;
% % %
% % % % load = matfile('pVSt_refine_x.mat');
% % % % t_save_x=load.t_save;
% % % % p_save_x=load.p_save;
% % % %
% % % % load = matfile('pVSt_refine_y.mat');
% % % % t_save_y=load.t_save;
% % % % p_save_y=load.p_save;
% % % %
% % % load = matfile('pVSt_refine_t.mat');
% % % t_save_t=load.t_save;
% % % p_save_t=load.p_save;
% % %
% % % figure
% % % hold on
% % % plot(Eclipse_time,Eclipse_pressure);
% % % plot(t_save,p_save);
% % % % plot(t_save_x,p_save_x);
% % % % plot(t_save_y,p_save_y);
% % % plot(t_save_t,p_save_t);
% % % % legend('Base case','Refine in x','Refine in y','Refine in t');
% % % legend('Eclipse','Base case','Refine in t');
% % % title('Refinement Study');
% % % xlabel('time (day)');
% % % ylabel('pressure (psi)');
% % % hold off
% % % % *************************************************************************
% %
% % plot_row_num=3;
% % plot_col_num=3;
% % tol_plot=plot_row_num*plot_col_num;
% % [p,q,r]=size(P_field_all);
% % interval=floor(r/(tol_plot-1));
% % time_plot_index=[1,1+interval:interval:1+interval*(tol_plot-2),r];
% % time_plot=zeros(1,length(time_plot_index));
% %
% % for i=1:1:length(time_plot_index)
% % time_plot(i)=t_vec(time_plot_index(i));
% % end
% %
% % % figure
% % % for i=1:1:plot_row_num*plot_col_num
% % % subplot(plot_row_num,plot_col_num,i)
% % % imagesc(P_field_all(:,:,time_plot_index(i))/6894.76);
% % % title(sprintf('t = %d days',time_plot(i)/(24*60*60)));
% % % colorbar
% % % end
% %
% % figure
% % imagesc(Pn_field/6894.76)
% % title('Pressure map for the entire field at the end of the simulation run, t = 400 days');
% % colorbar
