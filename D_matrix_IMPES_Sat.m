function [D_matrix,d11_matrix,d12_matrix,d21_matrix,d22_matrix]=D_matrix_IMPES_Sat(info,P_matrix_n,Pressure_matrix_n1,delta_t_n) % can get rid of everything
[Pressure_n,Saturation_n]=to2D(info,P_matrix_n);
Pressure_n1=transpose(reshape(Pressure_matrix_n1,[info.Nx,info.Ny]));

[info.Ny,info.Nx]=size(Pressure_n);
% D_matrix=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
D_matrix=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
d11_matrix=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
d12_matrix=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
d21_matrix=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
d22_matrix=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
for i=1:info.Ny
    for j=1:info.Nx
        Pn=Pressure_n(i,j);
        Sn=Saturation_n(i,j);
        Pn1=Pressure_n1(i,j);
%         Sn1=Saturation_n1(i,j);
        if (Pn==Pn1)
            bo_small_prime=0;
            bg_small_prime=0;
            Rs_prime=0; 
        else
            bo_small_prime=(bo_small(info,Pn1)-bo_small(info,Pn))/(Pn1-Pn);
            bg_small_prime=(bg_small(info,Pn1)-bg_small(info,Pn))/(Pn1-Pn);
            Rs_prime=(Rs(info,Pn1)-Rs(info,Pn))/(Pn1-Pn);            
        end
        d21=(1-Sn)*(info.phi*bo_small_prime)*(info.V/(5.615*delta_t_n));
        d22=-(info.phi*bo_small(info,Pn1))*(info.V/(5.615*delta_t_n));
        d11=(info.phi*bo_small(info,Pn1)*(1-Sn)*Rs_prime+Rs(info,Pn)*(1-Sn)*info.phi*bo_small_prime+Sn*info.phi*bg_small_prime)*(info.V/(5.615*delta_t_n));
        d12=(info.phi*bg_small(info,Pn1)-info.phi*bo_small(info,Pn1)*Rs(info,Pn))*(info.V/(5.615*delta_t_n));
        
        % store d
        d11_matrix(i,j)=d11;
        d12_matrix(i,j)=d12;
        d21_matrix(i,j)=d21;
        d22_matrix(i,j)=d22;
        
        l=(i-1)*info.Nx+j;
        %         D_matrix(2*l-1,2*l-1)= d11*(info.V/(5.615*info.delta_t));
        %         D_matrix(2*l-1,2*l)= d12*(info.V/(5.615*info.delta_t));
        %         D_matrix(2*l,2*l-1)= d21*(info.V/(5.615*info.delta_t));
        %         D_matrix(2*l,2*l)= d22*(info.V/(5.615*info.delta_t));
        D_matrix(l,l)= (d21/d22);
    end
end
% spy(D_matrix)