function P_matrix =P_matrix(P_res_orig,S_res_orig)
P_res=reshape(P_res_orig.',[],1);
S_res=reshape(S_res_orig.',[],1);
% l=(i-1)*nx+j;
[Ny,Nx]=size(P_res_orig);
P_matrix=zeros(2*Nx*Ny,1);
for i=1:1:Nx*Ny
P_matrix(2*i-1,1)=P_res(i,1);
P_matrix(2*i,1)=S_res(i,1);
end
end
