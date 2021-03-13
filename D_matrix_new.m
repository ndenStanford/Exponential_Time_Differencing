function D_matrix_new=D_matrix_new(info,P_matrix_n,P_matrix_n1)
[P_res_n,S_res_n]=to2D(info,P_matrix_n);
[P_res_n1,S_res_n1]=to2D(info,P_matrix_n1);

[info.Ny,info.Nx]=size(P_res_n);
D_matrix_new=zeros(2*info.Nx*info.Ny,1);
for i=1:info.Ny
    for j=1:info.Nx
        Pn=P_res_n(i,j);
        Sn=S_res_n(i,j);
        Pn1=P_res_n1(i,j);
        Sn1=S_res_n1(i,j);
        
        Dg=(-((1-Sn1)*info.phi*Rs(info,Pn1)*bo_small(info,Pn1)+info.phi*Sn1*bg_small(info,Pn1))+((1-Sn)*info.phi*Rs(info,Pn)*bo_small(info,Pn)+info.phi*Sn*bg_small(info,Pn)))*(info.V/(5.615*info.delta_t));        
        Do=-(info.phi*(1-Sn1)*bo_small(info,Pn1)-info.phi*(1-Sn)*bo_small(info,Pn))*(info.V/(5.615*info.delta_t));
        
        
        l=(i-1)*info.Nx+j;
        D_matrix_new(2*l-1,1)=Dg;
        D_matrix_new(2*l,1)=Do;      
        
    end
end