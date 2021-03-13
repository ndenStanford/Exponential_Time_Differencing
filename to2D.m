function [P_res_2D,S_res_2D]=to2D(info,P_matrix)
% extract p and s from P_matrix
% info.Nx=23;
% info.Ny=23;
% Nx=7;
% Ny=7;

P_res_1D=zeros((info.Nx)*(info.Ny),1);
S_res_1D=zeros(info.Nx*info.Ny,1);
for i=1:info.Nx*info.Ny
P_res_1D(i,1)=P_matrix(2*i-1);
S_res_1D(i,1)=P_matrix(2*i);
end
P_res_2D=transpose(reshape(P_res_1D,[info.Nx,info.Ny]));
S_res_2D=transpose(reshape(S_res_1D,[info.Nx,info.Ny]));
end