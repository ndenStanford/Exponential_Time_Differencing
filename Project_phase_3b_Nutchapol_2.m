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


% test case
P_res_orig=[2400,2425,2450;2350,2375,2400];
S_res_orig=[0,0.15,0.25;0.35,0.45,0.55];
P_res_n1=[2405,2430,2455;2355,2380,2405];
S_res_n1=[0.15,0.25,0.35;0.45,0.55,0.65];



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
R_matrix=(T_matrix*P_matrix_n1-D_matrix*(P_matrix_n1-P_matrix_n));

% make J matrix
J_matrix=J_matrix(P_res_n1,S_res_n1);


% compare the result with the reference case
dR_matrix=abs(R_matrix-R_excel);
dR_rel_matrix=abs((R_matrix-R_excel)./R_excel);

[row,col]=size(J_matrix);
dJ_matrix=zeros(row,col);
dJ_rel_matrix=zeros(row,col);

for i=1:row
    for j=1:col
        if J_matrix(i,j)~=0
            dJ_matrix(i,j)=abs((J_matrix(i,j)-J_excel(i,j)));
        end
    end
end

for i=1:row
    for j=1:col
        if J_matrix(i,j)~=0
            dJ_rel_matrix(i,j)=abs((J_matrix(i,j)-J_excel(i,j))/J_excel(i,j));
        end
    end
end

% figure
% % imagesc(dR_matrix)
% plot(dR_matrix)
% ylabel('absolute error of the residual')
% xlabel('block number')
% title('absolute error of the residual: g');
% % colorbar
% 
% figure
% % imagesc(dR_rel_matrix)
% plot(dR_rel_matrix)
% ylabel('relative error of the residual')
% xlabel('block number')
% title('relative error of the residual: g');
% % colorbar
% 
% figure
% imagesc(sparse(dJ_matrix))
% title('absolute error of the jacobian: J');
% colorbar
% 
% figure
% imagesc(sparse(dJ_rel_matrix))
% title('relative error of the jacobian: J');
% colorbar



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
