function info=initializeParameters()
% unit conversion factor
info.alpha=0.001127;

% Grid and simulation parameters
info.Nx=23;
info.Ny=23;
% info.t=1;
info.delta_t=1;%1;
info.t_runtime=600;

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
% info.Pinit=2000;
info.P_res_orig=info.Pinit*ones(info.Ny,info.Nx);
info.Sinit=0;
info.S_res_orig=info.Sinit*ones(info.Ny,info.Nx);
info.P_res_n1=info.P_res_orig;
info.S_res_n1=info.S_res_orig;

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

% implicit only
info.delta_t_max=5;%50;
end