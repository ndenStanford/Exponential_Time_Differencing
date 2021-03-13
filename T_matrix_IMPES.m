function T_matrix=T_matrix_IMPES(info,P_matrix,d11_matrix,d12_matrix,d21_matrix,d22_matrix)
% [P_res_n1,S_res_n1]=to2D(info,P_matrix);

[Px_upwinded,Py_upwinded,Sx_upwinded,Sy_upwinded,Hx_oil,Hy_oil,Hx_gas,Hy_gas]=connection_list(info,P_matrix);
% calculate the necessary elements
% T_matrix=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
T_matrix=zeros(info.Nx*info.Ny,info.Nx*info.Ny);

for i=1:info.Ny
    for j=1:info.Nx
        %         Tx=kx_res_orig(i,j)*delta_y*Lz/delta_x;
        %         Ty=ky_res_orig(i,j)*delta_x*Lz/delta_y;
        % beta : i-1/2,j
        if i==1
            gamma_top=0;
            beta_top=0;
        else
            ky=2/(1/info.ky_res_orig(floor(i-1/2),j)+1/info.ky_res_orig(ceil(i-1/2),j));
            Ty=ky*info.delta_x*info.Lz/info.delta_y;
            gamma_top=info.alpha*Ty*Hy_oil(floor(i-1/2),j);
            beta_top=info.alpha*Ty*(Rs(info,Py_upwinded(floor(i-1/2),j))*Hy_oil(floor(i-1/2),j)+Hy_gas(floor(i-1/2),j));
            %             gamma_top=1;
            %             beta_top=1;
        end
        % beta : i+1/2,j
        if i==info.Ny
            gamma_bottom=0;
            beta_bottom=0;
        else
            ky=2/(1/info.ky_res_orig(floor(i+1/2),j)+1/info.ky_res_orig(ceil(i+1/2),j));
            Ty=ky*info.delta_x*info.Lz/info.delta_y;
            gamma_bottom=info.alpha*Ty*Hy_oil(floor(i+1/2),j);
            beta_bottom=info.alpha*Ty*(Rs(info,Py_upwinded(floor(i+1/2),j))*Hy_oil(floor(i+1/2),j)+Hy_gas(floor(i+1/2),j));
            %             gamma_bottom=1;
            %             beta_bottom=1;
        end
        % beta : i,j-1/2
        if j==1
            gamma_left=0;
            beta_left=0;
        else
            kx=(info.kx_res_orig(i,floor(j-1/2))+info.kx_res_orig(i,ceil(j-1/2)))/2;
            Tx=kx*info.delta_y*info.Lz/info.delta_x;
            gamma_left=info.alpha*Tx*Hx_oil(i,floor(j-1/2));
            beta_left=info.alpha*Tx*(Rs(info,Px_upwinded(i,floor(j-1/2)))*Hx_oil(i,floor(j-1/2))+Hx_gas(i,floor(j-1/2)));
            %             gamma_left=1;
            %             beta_left=1;
        end
        % beta : i,j+1/2
        if j==info.Nx
            gamma_right=0;
            beta_right=0;
        else
            kx=(info.kx_res_orig(i,floor(j+1/2))+info.kx_res_orig(i,ceil(j+1/2)))/2;
            Tx=kx*info.delta_y*info.Lz/info.delta_x;
            gamma_right=info.alpha*Tx*Hx_oil(i,floor(j+1/2));
            beta_right=info.alpha*Tx*(Rs(info,Px_upwinded(i,floor(j+1/2)))*Hx_oil(i,floor(j+1/2))+Hx_gas(i,floor(j+1/2)));
            %             gamma_right=1;
            %             beta_right=1;
        end
        l=(i-1)*info.Nx+j;
        %         T_matrix(2*l-1,2*l-1)= -(beta_left+beta_right+beta_bottom+beta_top);
        %         T_matrix(2*l,2*l-1)= -(gamma_left+gamma_right+gamma_bottom+gamma_top);
        T_matrix(l,l)=-(-(d22_matrix(i,j)/d12_matrix(i,j))*(beta_left+beta_right+beta_bottom+beta_top)+(gamma_left+gamma_right+gamma_bottom+gamma_top));
        %         T_matrix(2*l-1,2*l-1)= 1;
        %         T_matrix(2*l,2*l-1)= 1;
        if l-1>=1
            %             T_matrix(2*l-1,2*l-3)= beta_left;
            %             T_matrix(2*l,2*l-3)= gamma_left;
            T_matrix(l,l-1)= -(d22_matrix(i,j)/d12_matrix(i,j))*beta_left+gamma_left;
        end
        if l+1<=info.Nx*info.Ny
            %             T_matrix(2*l-1,2*l+1)= beta_right;
            %             T_matrix(2*l,2*l+1)= gamma_right;
            T_matrix(l,l+1)= -(d22_matrix(i,j)/d12_matrix(i,j))*beta_right+gamma_right;
        end
        if l-info.Nx>=1
            %             T_matrix(2*l-1,2*l-1-2*info.Nx)= beta_top;
            %             T_matrix(2*l,2*l-1-2*info.Nx)= gamma_top;
            T_matrix(l,l-info.Nx)= -(d22_matrix(i,j)/d12_matrix(i,j))*beta_top+gamma_top;
        end
        if l+info.Nx<=info.Nx*info.Ny
            %             T_matrix(2*l-1,2*l-1+2*info.Nx)= beta_bottom;
            %             T_matrix(2*l,2*l-1+2*info.Nx)= gamma_bottom;
            T_matrix(l,l+info.Nx)= -(d22_matrix(i,j)/d12_matrix(i,j))*beta_bottom+gamma_bottom;
        end
        %         a=T_matrix(2*l-1,2*l-1);
        %         b=T_matrix(2*l,2*l-1);
        
        
    end
end
% spy(T_matrix)

end