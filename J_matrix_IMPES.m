% get rid of d matrix
function J_matrix=J_matrix_IMPES(info,P_matrix,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix,delta_t_n)
[P_res_n1,S_res_n1]=to2D(info,P_matrix);
[Px_upwinded,Py_upwinded,Sx_upwinded,Sy_upwinded,Hx_oil,Hy_oil,Hx_gas,Hy_gas]=connection_list(info,P_matrix);
% contribution from D matrix
[info.Ny,info.Nx]=size(P_res_n1);
% J_D=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
J_D=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
for i=1:info.Ny
    for j=1:info.Nx
        Pn1=P_res_n1(i,j);
        Sn1=S_res_n1(i,j);
        d21=(1-Sn1)*(info.phi*bo_small_prime(info,Pn1));
        d22=-(info.phi*bo_small(info,Pn1));
        d11=(1-Sn1)*(info.phi*bo_small(info,Pn1)*Rs_prime(info,Pn1)+info.phi*bo_small_prime(info,Pn1)*Rs(info,Pn1))+Sn1*(info.phi*bg_small_prime(info,Pn1));
        d12=-info.phi*bo_small(info,Pn1)*Rs(info,Pn1)+info.phi*bg_small(info,Pn1);
        l=(i-1)*info.Nx+j;
        %         J_D(2*l-1,2*l-1)= d11*(info.V/(5.615*info.delta_t));
        %         J_D(2*l-1,2*l)= d12*(info.V/(5.615*info.delta_t));
        %         J_D(2*l,2*l-1)= d21*(info.V/(5.615*info.delta_t));
        %         J_D(2*l,2*l)= d22*(info.V/(5.615*info.delta_t));
        J_D(l,l)= (-(d22_matrix(i,j)/d12_matrix(i,j))*d11+d21)*(info.V/(5.615*delta_t_n));
        
    end
end


% contribution from T matrix
% J_T=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
J_T=zeros(info.Nx*info.Ny,info.Nx*info.Ny);

for i=1:info.Ny
    for j=1:info.Nx
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
            ky=2/(1/info.ky_res_orig(floor(i-1/2),j)+1/info.ky_res_orig(ceil(i-1/2),j));
            Ty=ky*info.delta_x*info.Lz/info.delta_y;
            P_upwinded=Py_upwinded(floor(i-1/2),j);
            gamma_oil_top=info.alpha*Ty*Hy_oil(floor(i-1/2),j);
            gamma_gas_top=info.alpha*Ty*Hy_gas(floor(i-1/2),j);
            Rs_top=Rs(info,P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i-1,j)
                dgamma_oil_dp_top=0;
                dgamma_oil_ds_top=0;
                dgamma_oil_dp_top_top=info.alpha*Ty*kro(S_res_n1(i-1,j))*(bo_small_prime(info,P_res_n1(i-1,j))/info.mu_oil);
                dgamma_oil_ds_top_top=info.alpha*Ty*kro_prime(S_res_n1(i-1,j))/(Bo(info,P_res_n1(i-1,j))*info.mu_oil);
                dgamma_gas_dp_top=0;
                dgamma_gas_ds_top=0;
                dgamma_gas_dp_top_top=info.alpha*Ty*krg_new(S_res_n1(i-1,j))*(bg_small_prime(info,P_res_n1(i-1,j))/(mu_g(info,P_res_n1(i-1,j)))-(mu_g_prime(info,P_res_n1(i-1,j))*bg_small(info,P_res_n1(i-1,j))/(mu_g(info,P_res_n1(i-1,j)))^2));
                dgamma_gas_ds_top_top=info.alpha*Ty*krg_prime(S_res_n1(i-1,j))/(Bg(info,P_res_n1(i-1,j))*mu_g(info,P_res_n1(i-1,j)));
                Rs_prime_top=0;
                Rs_prime_top_top=Rs_prime(info,P_res_n1(i-1,j));
                mu_g_prime_top=0;
                mu_g_prime_top_top=mu_g_prime(info,P_res_n1(i-1,j));
                
            else
                dgamma_oil_dp_top=info.alpha*Ty*kro(S_res_n1(i,j))*(bo_small_prime(info,P_res_n1(i,j))/info.mu_oil);
                dgamma_oil_ds_top=info.alpha*Ty*kro_prime(S_res_n1(i,j))/(Bo(info,P_res_n1(i,j))*info.mu_oil);
                dgamma_oil_dp_top_top=0;
                dgamma_oil_ds_top_top=0;
                dgamma_gas_dp_top=info.alpha*Ty*krg_new(S_res_n1(i,j))*(bg_small_prime(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))-(mu_g_prime(info,P_res_n1(i,j))*bg_small(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))^2));
                dgamma_gas_ds_top=info.alpha*Ty*krg_prime(S_res_n1(i,j))/(Bg(info,P_res_n1(i,j))*mu_g(info,P_res_n1(i,j)));
                dgamma_gas_dp_top_top=0;
                dgamma_gas_ds_top_top=0;
                Rs_prime_top=Rs_prime(info,P_res_n1(i,j));
                Rs_prime_top_top=0;
                mu_g_prime_top=mu_g_prime(info,P_res_n1(i,j));
                mu_g_prime_top_top=0;
            end
        end
        % dgamma_o/dP : i+1/2,j
        if i==info.Ny
            gamma_oil_bottom=0;
            dgamma_oil_dp_bottom=0;
            dgamma_oil_ds_bottom=0;
            gamma_gas_bottom=0;
            dgamma_gas_dp_bottom=0;
            dgamma_gas_ds_bottom=0;
            Rs_bottom=0;
            Rs_prime_bottom=0;
        else
            ky=2/(1/info.ky_res_orig(floor(i+1/2),j)+1/info.ky_res_orig(ceil(i+1/2),j));
            Ty=ky*info.delta_x*info.Lz/info.delta_y;
            P_upwinded=Py_upwinded(floor(i+1/2),j);
            gamma_oil_bottom=info.alpha*Ty*Hy_oil(floor(i+1/2),j);
            gamma_gas_bottom=info.alpha*Ty*Hy_gas(floor(i+1/2),j);
            Rs_bottom=Rs(info,P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i+1,j)
                dgamma_oil_dp_bottom=0;
                dgamma_oil_ds_bottom=0;
                dgamma_oil_dp_bottom_bottom=info.alpha*Ty*kro(S_res_n1(i+1,j))*(bo_small_prime(info,P_res_n1(i+1,j))/info.mu_oil);
                dgamma_oil_ds_bottom_bottom=info.alpha*Ty*kro_prime(S_res_n1(i+1,j))/(Bo(info,P_res_n1(i+1,j))*info.mu_oil);
                dgamma_gas_dp_bottom=0;
                dgamma_gas_ds_bottom=0;
                dgamma_gas_dp_bottom_bottom=info.alpha*Ty*krg_new(S_res_n1(i+1,j))*(bg_small_prime(info,P_res_n1(i+1,j))/(mu_g(info,P_res_n1(i+1,j)))-(mu_g_prime(info,P_res_n1(i+1,j))*bg_small(info,P_res_n1(i+1,j))/(mu_g(info,P_res_n1(i+1,j)))^2));
                dgamma_gas_ds_bottom_bottom=info.alpha*Ty*kro_prime(S_res_n1(i+1,j))/(Bg(info,P_res_n1(i+1,j))*mu_g(info,P_res_n1(i+1,j)));
                Rs_prime_bottom=0;
                Rs_prime_bottom_bottom=Rs_prime(info,P_res_n1(i+1,j));
                mu_g_prime_bottom=0;
                mu_g_prime_bottom_bottom=mu_g_prime(info,P_res_n1(i+1,j));
                
            else
                dgamma_oil_dp_bottom=info.alpha*Ty*kro(S_res_n1(i,j))*(bo_small_prime(info,P_res_n1(i,j))/info.mu_oil);
                dgamma_oil_ds_bottom=info.alpha*Ty*kro_prime(S_res_n1(i,j))/(Bo(info,P_res_n1(i,j))*info.mu_oil);
                dgamma_oil_dp_bottom_bottom=0;
                dgamma_oil_ds_bottom_bottom=0;
                dgamma_gas_dp_bottom=info.alpha*Ty*krg_new(S_res_n1(i,j))*(bg_small_prime(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))-(mu_g_prime(info,P_res_n1(i,j))*bg_small(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))^2));
                dgamma_gas_ds_bottom=info.alpha*Ty*krg_prime(S_res_n1(i,j))/(Bg(info,P_res_n1(i,j))*mu_g(info,P_res_n1(i,j)));
                dgamma_gas_dp_bottom_bottom=0;
                dgamma_gas_ds_bottom_bottom=0;
                Rs_prime_bottom=Rs_prime(info,P_res_n1(i,j));
                Rs_prime_bottom_bottom=0;
                mu_g_prime_bottom=mu_g_prime(info,P_res_n1(i,j));
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
            kx=(info.kx_res_orig(i,floor(j-1/2))+info.kx_res_orig(i,ceil(j-1/2)))/2;
            Tx=kx*info.delta_y*info.Lz/info.delta_x;
            P_upwinded=Px_upwinded(i,floor(j-1/2));
            gamma_oil_left=info.alpha*Tx*Hx_oil(i,floor(j-1/2));
            gamma_gas_left=info.alpha*Tx*Hx_gas(i,floor(j-1/2));
            Rs_left=Rs(info,P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i,j-1)
                dgamma_oil_dp_left=0;
                dgamma_oil_ds_left=0;
                dgamma_oil_dp_left_left=info.alpha*Tx*kro(S_res_n1(i,j-1))*(bo_small_prime(info,P_res_n1(i,j-1))/info.mu_oil);
                dgamma_oil_ds_left_left=info.alpha*Tx*kro_prime(S_res_n1(i,j-1))/(Bo(info,P_res_n1(i,j-1))*info.mu_oil);
                dgamma_gas_dp_left=0;
                dgamma_gas_ds_left=0;
                dgamma_gas_dp_left_left=info.alpha*Tx*krg_new(S_res_n1(i,j-1))*(bg_small_prime(info,P_res_n1(i,j-1))/(mu_g(info,P_res_n1(i,j-1)))-(mu_g_prime(info,P_res_n1(i,j-1))*bg_small(info,P_res_n1(i,j-1))/(mu_g(info,P_res_n1(i,j-1)))^2));
                dgamma_gas_ds_left_left=info.alpha*Tx*krg_prime(S_res_n1(i,j-1))/(Bg(info,P_res_n1(i,j-1))*mu_g(info,P_res_n1(i,j-1)));
                Rs_prime_left=0;
                Rs_prime_left_left=Rs_prime(info,P_res_n1(i,j-1));
                mu_g_prime_left=0;
                mu_g_prime_left_left=mu_g_prime(info,P_res_n1(i,j-1));
            else
                dgamma_oil_dp_left=info.alpha*Tx*kro(S_res_n1(i,j))*(bo_small_prime(info,P_res_n1(i,j))/info.mu_oil);
                dgamma_oil_ds_left=info.alpha*Tx*kro_prime(S_res_n1(i,j))/(Bo(info,P_res_n1(i,j))*info.mu_oil);
                dgamma_oil_dp_left_left=0;
                dgamma_oil_ds_left_left=0;
                dgamma_gas_dp_left=info.alpha*Tx*krg_new(S_res_n1(i,j))*(bg_small_prime(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))-(mu_g_prime(info,P_res_n1(i,j))*bg_small(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))^2));
                dgamma_gas_ds_left=info.alpha*Tx*krg_prime(S_res_n1(i,j))/(Bg(info,P_res_n1(i,j))*mu_g(info,P_res_n1(i,j)));
                dgamma_gas_dp_left_left=0;
                dgamma_gas_ds_left_left=0;
                Rs_prime_left=Rs_prime(info,P_res_n1(i,j));
                Rs_prime_left_left=0;
                mu_g_prime_left=mu_g_prime(info,P_res_n1(i,j));
                mu_g_prime_left_left=0;
            end
        end
        % dgamma_o/dP : i,j+1/2
        if j==info.Nx
            gamma_oil_right=0;
            dgamma_oil_dp_right=0;
            dgamma_oil_ds_right=0;
            gamma_gas_right=0;
            dgamma_gas_dp_right=0;
            dgamma_gas_ds_right=0;
            Rs_right=0;
            Rs_prime_right=0;
        else
            kx=(info.kx_res_orig(i,floor(j+1/2))+info.kx_res_orig(i,ceil(j+1/2)))/2;
            Tx=kx*info.delta_y*info.Lz/info.delta_x;
            P_upwinded=Px_upwinded(i,floor(j+1/2));
            gamma_oil_right=info.alpha*Tx*Hx_oil(i,floor(j+1/2));
            gamma_gas_right=info.alpha*Tx*Hx_gas(i,floor(j+1/2));
            Rs_right=Rs(info,P_upwinded);
            if P_res_n1(i,j)<=P_res_n1(i,j+1)
                dgamma_oil_dp_right=0;
                dgamma_oil_ds_right=0;
                dgamma_oil_dp_right_right=info.alpha*Tx*kro(S_res_n1(i,j+1))*(bo_small_prime(info,P_res_n1(i,j+1))/info.mu_oil);
                dgamma_oil_ds_right_right=info.alpha*Tx*kro_prime(S_res_n1(i,j+1))/(Bo(info,P_res_n1(i,j+1))*info.mu_oil);
                dgamma_gas_dp_right=0;
                dgamma_gas_ds_right=0;
                dgamma_gas_dp_right_right=info.alpha*Tx*krg_new(S_res_n1(i,j+1))*(bg_small_prime(info,P_res_n1(i,j+1))/(mu_g(info,P_res_n1(i,j+1)))-(mu_g_prime(info,P_res_n1(i,j+1))*bg_small(info,P_res_n1(i,j+1))/(mu_g(info,P_res_n1(i,j+1)))^2));
                dgamma_gas_ds_right_right=info.alpha*Tx*krg_prime(S_res_n1(i,j+1))/(Bg(info,P_res_n1(i,j+1))*mu_g(info,P_res_n1(i,j+1)));
                Rs_prime_right=0;
                Rs_prime_right_right=Rs_prime(info,P_res_n1(i,j+1));
                mu_g_prime_right=0;
                mu_g_prime_right_right=mu_g_prime(info,P_res_n1(i,j));
            else
                dgamma_oil_dp_right=info.alpha*Tx*kro(S_res_n1(i,j))*(bo_small_prime(info,P_res_n1(i,j))/info.mu_oil);
                dgamma_oil_ds_right=info.alpha*Tx*kro_prime(S_res_n1(i,j))/(Bo(info,P_res_n1(i,j))*info.mu_oil);
                dgamma_oil_dp_right_right=0;
                dgamma_oil_ds_right_right=0;
                dgamma_gas_dp_right=info.alpha*Tx*krg_new(S_res_n1(i,j))*(bg_small_prime(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))-(mu_g_prime(info,P_res_n1(i,j))*bg_small(info,P_res_n1(i,j))/(mu_g(info,P_res_n1(i,j)))^2));
                dgamma_gas_ds_right=info.alpha*Tx*krg_prime(S_res_n1(i,j))/(Bg(info,P_res_n1(i,j))*mu_g(info,P_res_n1(i,j)));
                dgamma_gas_dp_right_right=0;
                dgamma_gas_ds_right_right=0;
                Rs_prime_right=Rs_prime(info,P_res_n1(i,j));
                Rs_prime_right_right=0;
                mu_g_prime_right=mu_g_prime(info,P_res_n1(i,j));
                mu_g_prime_right_right=0;
            end
        end
        
        P_mid=P_res_n1(i,j);
        if i-1>=1
            P_top=P_res_n1(i-1,j);
        else
            P_top=0;
        end
        if i+1<=info.Ny
            P_bottom=P_res_n1(i+1,j);
        else
            P_bottom=0;
        end
        if j-1>=1
            P_left=P_res_n1(i,j-1);
        else
            P_left=0;
        end
        if j+1<=info.Nx
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
        
        l=(i-1)*info.Nx+j;
        %         J_T(2*l,2*l-1)= dR_dP_oil;
        %         J_T(2*l,2*l)= dR_dS_oil;
        %         J_T(2*l-1,2*l-1)= dR_dP_gas;
        %         J_T(2*l-1,2*l)= dR_dS_gas;
        J_T(l,l)=-(d22_matrix(i,j)/d12_matrix(i,j))*dR_dP_gas+dR_dP_oil;
        
        
        if j>1
            dR_dP_oil_left=(dgamma_oil_dp_left_left*(P_left-P_mid)+gamma_oil_left);
            dR_dS_oil_left=dgamma_oil_ds_left_left*(P_left-P_mid);
            %             J_T(2*l,2*l-3)= dR_dP_oil_left;
            %             J_T(2*l,2*l-2)= dR_dS_oil_left;
            dR_dP_gas_left=((Rs_left*gamma_oil_left+gamma_gas_left)+(Rs_left*dgamma_oil_dp_left_left+gamma_oil_left*Rs_prime_left_left+dgamma_gas_dp_left_left)*(P_left-P_mid));
            dR_dS_gas_left=((Rs_left*dgamma_oil_ds_left_left+dgamma_gas_ds_left_left)*(P_left-P_mid));
            %             J_T(2*l-1,2*l-3)= dR_dP_gas_left;
            %             J_T(2*l-1,2*l-2)= dR_dS_gas_left;
            J_T(l,l-1)=-(d22_matrix(i,j)/d12_matrix(i,j))*dR_dP_gas_left+dR_dP_oil_left;
        end
        
        if j<info.Nx
            dR_dP_oil_right=(dgamma_oil_dp_right_right*(P_right-P_mid)+gamma_oil_right);
            dR_dS_oil_right=dgamma_oil_ds_right_right*(P_right-P_mid);
            %             J_T(2*l,2*l+1)= dR_dP_oil_right;
            %             J_T(2*l,2*l+2)= dR_dS_oil_right;
            dR_dP_gas_right=((Rs_right*gamma_oil_right+gamma_gas_right)+(Rs_right*dgamma_oil_dp_right_right+gamma_oil_right*Rs_prime_right_right+dgamma_gas_dp_right_right)*(P_right-P_mid));
            dR_dS_gas_right=((Rs_right*dgamma_oil_ds_right_right+dgamma_gas_ds_right_right)*(P_right-P_mid));
            %             J_T(2*l-1,2*l+1)= dR_dP_gas_right;
            %             J_T(2*l-1,2*l+2)= dR_dS_gas_right;
            J_T(l,l+1)=-(d22_matrix(i,j)/d12_matrix(i,j))*dR_dP_gas_right+dR_dP_oil_right;
        end
        
        if i>1
            dR_dP_oil_top=(dgamma_oil_dp_top_top*(P_top-P_mid)+gamma_oil_top);
            dR_dS_oil_top=dgamma_oil_ds_top_top*(P_top-P_mid);
            %             J_T(2*l,2*l-1-2*info.Nx)= dR_dP_oil_top;
            %             J_T(2*l,2*l-2*info.Nx)= dR_dS_oil_top;
            dR_dP_gas_top=((Rs_top*gamma_oil_top+gamma_gas_top)+(Rs_top*dgamma_oil_dp_top_top+gamma_oil_top*Rs_prime_top_top+dgamma_gas_dp_top_top)*(P_top-P_mid));
            dR_dS_gas_top=((Rs_top*dgamma_oil_ds_top_top+dgamma_gas_ds_top_top)*(P_top-P_mid));
            %             J_T(2*l-1,2*l-1-2*info.Nx)= dR_dP_gas_top;
            %             J_T(2*l-1,2*l-2*info.Nx)= dR_dS_gas_top;
            J_T(l,l-info.Nx)=-(d22_matrix(i,j)/d12_matrix(i,j))*dR_dP_gas_top+dR_dP_oil_top;
        end
        
        if i<info.Ny
            dR_dP_oil_bottom=(dgamma_oil_dp_bottom_bottom*(P_bottom-P_mid)+gamma_oil_bottom);
            dR_dS_oil_bottom=dgamma_oil_ds_bottom_bottom*(P_bottom-P_mid);
            %             J_T(2*l,2*l-1+2*info.Nx)= dR_dP_oil_bottom;
            %             J_T(2*l,2*l+2*info.Nx)= dR_dS_oil_bottom;
            dR_dP_gas_bottom=((Rs_bottom*gamma_oil_bottom+gamma_gas_bottom)+(Rs_bottom*dgamma_oil_dp_bottom_bottom+gamma_oil_bottom*Rs_prime_bottom_bottom+dgamma_gas_dp_bottom_bottom)*(P_bottom-P_mid));
            dR_dS_gas_bottom=((Rs_bottom*dgamma_oil_ds_bottom_bottom+dgamma_gas_ds_bottom_bottom)*(P_bottom-P_mid));
            %             J_T(2*l-1,2*l-1+2*info.Nx)= dR_dP_gas_bottom;
            %             J_T(2*l-1,2*l+2*info.Nx)= dR_dS_gas_bottom;
            J_T(l,l+info.Nx)=-(d22_matrix(i,j)/d12_matrix(i,j))*dR_dP_gas_bottom+dR_dP_oil_bottom;
        end
        
        
    end
    
    % producer well
    % contribution from Q matrix
%     J_Q1=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
%     J_Q2=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
    J_Q1=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
    J_Q2=zeros(info.Nx*info.Ny,info.Nx*info.Ny);
    p_well_block=P_res_n1(info.yw,info.xw);
    s_well_block=S_res_n1(info.yw,info.xw);
    % ro=0.28*(((ky/kx)^(1/2)*delta_x^2+(kx/ky)^(1/2)*delta_y^2)^(1/2))/((ky/kx)^(1/4)+(kx/ky)^(1/4));
    % WI=(alpha*2*pi*(kx*ky)^(1/2)*Lz/(log(ro/rw)+s));
    To=info.WI*(kro(s_well_block)/(info.mu_oil*Bo(info,p_well_block)));
    Tg=info.WI*(krg_new(s_well_block)/(mu_g(info,p_well_block)*Bg(info,p_well_block)));
    
    % case 2 rate control
    % oil flow
    dRo_dp_2=0;
    dRo_ds_2=0;
    
    % case 1 pressure control
    % oil flow
    % dTo_dp=WI*kro(s_well_block)*(-1/(mu_oil*Bo(p_well_block))^2)*(mu_oil*Bo_prime(p_well_block));
    % dRo_dp=-(To+(p_well_block-p_bhp)*dTo_dp);
    % dTo_ds=WI*kro_prime(s_well_block)/(mu_oil*Bo(p_well_block));
    % dRo_ds=-((p_well_block-p_bhp)*dTo_ds);
    dTo_dp=info.WI*kro(s_well_block)*(1/info.mu_oil*bo_small_prime(info,p_well_block));
    dRo_dp_1=-(To+(p_well_block-info.p_bhp)*dTo_dp);
    dTo_ds=info.WI*kro_prime(s_well_block)*(bo_small(info,p_well_block)/info.mu_oil);
    dRo_ds_1=-((p_well_block-info.p_bhp)*dTo_ds);
    
    % case 2 rate control
    % gas flow
    dTg_dp=info.WI*krg_new(s_well_block)*(1/mu_g(info,p_well_block)*bg_small_prime(info,p_well_block)-bg_small(info,p_well_block)*mu_g_prime(info,p_well_block)/(mu_g(info,p_well_block))^2);
    dRg_dp_2=-(1/To*dTg_dp-(Tg/To^2)*dTo_dp+Rs_prime(info,p_well_block))*info.qwo_control;
    dTg_ds=info.WI*bg_small(info,p_well_block)/mu_g(info,p_well_block)*krg_prime(s_well_block);
    dRg_ds_2=-(1/To*dTg_ds-(Tg/To^2)*dTo_ds)*info.qwo_control;
    
    % case 1 pressure control
    % gas flow
    % dTg_dp=WI*krg(s_well_block)*(-1/(mu_g(p_well_block)*Bg(p_well_block))^2)*(mu_g_prime(p_well_block)*Bg(p_well_block)+mu_g(p_well_block)*Bg_prime(p_well_block));
    % dRg_dp=-((Rs_prime(p_well_block)*To+Rs(p_well_block)*dTo_dp)*(p_well_block-p_bhp)+Rs(p_well_block)*To+Tg+(p_well_block-p_bhp)*dTg_dp);
    % dTg_ds=WI*krg_prime(s_well_block)/(mu_g(p_well_block)*Bg(p_well_block));
    % dRg_ds=-(Rs(p_well_block)*(p_well_block-p_bhp)*dTo_ds+(p_well_block-p_bhp)*dTg_ds);
    dRg_dp_1=-((Rs_prime(info,p_well_block)*To+Rs(info,p_well_block)*dTo_dp)*(p_well_block-info.p_bhp)+Rs(info,p_well_block)*To+Tg+(p_well_block-info.p_bhp)*dTg_dp);
    dRg_ds_1=-(Rs(info,p_well_block)*(p_well_block-info.p_bhp)*dTo_ds+(p_well_block-info.p_bhp)*dTg_ds);
    
    
    % case 1 pressure control
    % l=(info.ywa-1)*info.Nx+info.xwa;
    l=(info.yw-1)*info.Nx+info.xw;
%     J_Q1(2*l-1,2*l-1)=dRg_dp_1;
%     J_Q1(2*l-1,2*l)=dRg_ds_1;
%     J_Q1(2*l,2*l-1)=dRo_dp_1;
%     J_Q1(2*l,2*l)=dRo_ds_1;
    J_Q1(l,l)=-(d22_matrix(i,j)/d12_matrix(i,j))*dRg_dp_1+dRo_dp_1;
    
    % case 2 rate control
    % l=(info.ywa-1)*info.Nx+info.xwa;
    l=(info.yw-1)*info.Nx+info.xw;
%     J_Q2(2*l-1,2*l-1)=dRg_dp_2;
%     J_Q2(2*l-1,2*l)=dRg_ds_2;
%     J_Q2(2*l,2*l-1)=dRo_dp_2;
%     J_Q2(2*l,2*l)=dRo_ds_2;
    J_Q2(l,l)=-(d22_matrix(i,j)/d12_matrix(i,j))*dRg_dp_2+dRo_dp_2;
    
    
    % % injector well
    % % contribution from Q matrix
    % % J_Q1=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
    % % J_Q2=zeros(2*info.Nx*info.Ny,2*info.Nx*info.Ny);
    % p_well_block=P_res_n1(info.ywa,info.xwa);
    % s_well_block=S_res_n1(info.ywa,info.xwa);
    % % ro=0.28*(((ky/kx)^(1/2)*delta_x^2+(kx/ky)^(1/2)*delta_y^2)^(1/2))/((ky/kx)^(1/4)+(kx/ky)^(1/4));
    % % WI=(alpha*2*pi*(kx*ky)^(1/2)*Lz/(log(ro/rw)+s));
    % To=info.WI*(kro(s_well_block)/(info.mu_oil*Bo(info,p_well_block)));
    % Tg=info.WI*(krg_new(s_well_block)/(mu_g(info,p_well_block)*Bg(info,p_well_block)));
    %
    % % case 2 rate control
    % % oil flow
    % dRo_dp_2=0;
    % dRo_ds_2=0;
    %
    % % case 1 pressure control
    % % oil flow
    % % dTo_dp=WI*kro(s_well_block)*(-1/(mu_oil*Bo(p_well_block))^2)*(mu_oil*Bo_prime(p_well_block));
    % % dRo_dp=-(To+(p_well_block-p_bhp)*dTo_dp);
    % % dTo_ds=WI*kro_prime(s_well_block)/(mu_oil*Bo(p_well_block));
    % % dRo_ds=-((p_well_block-p_bhp)*dTo_ds);
    % dTo_dp=info.WI*kro(s_well_block)*(1/info.mu_oil*bo_small_prime(info,p_well_block));
    % dRo_dp_1=-(To+(p_well_block-info.p_bhp)*dTo_dp);
    % dTo_ds=info.WI*kro_prime(s_well_block)*(bo_small(info,p_well_block)/info.mu_oil);
    % dRo_ds_1=-((p_well_block-info.p_bhp)*dTo_ds);
    %
    % % case 2 rate control
    % % gas flow
    % dTg_dp=info.WI*krg_new(s_well_block)*(1/mu_g(info,p_well_block)*bg_small_prime(info,p_well_block)-bg_small(info,p_well_block)*mu_g_prime(info,p_well_block)/(mu_g(info,p_well_block))^2);
    % dRg_dp_2=-(1/To*dTg_dp-(Tg/To^2)*dTo_dp+Rs_prime(info,p_well_block))*info.qwo_control;
    % dTg_ds=info.WI*bg_small(info,p_well_block)/mu_g(info,p_well_block)*krg_prime(s_well_block);
    % dRg_ds_2=-(1/To*dTg_ds-(Tg/To^2)*dTo_ds)*info.qwo_control;
    %
    % % case 1 pressure control
    % % gas flow
    % % dTg_dp=WI*krg(s_well_block)*(-1/(mu_g(p_well_block)*Bg(p_well_block))^2)*(mu_g_prime(p_well_block)*Bg(p_well_block)+mu_g(p_well_block)*Bg_prime(p_well_block));
    % % dRg_dp=-((Rs_prime(p_well_block)*To+Rs(p_well_block)*dTo_dp)*(p_well_block-p_bhp)+Rs(p_well_block)*To+Tg+(p_well_block-p_bhp)*dTg_dp);
    % % dTg_ds=WI*krg_prime(s_well_block)/(mu_g(p_well_block)*Bg(p_well_block));
    % % dRg_ds=-(Rs(p_well_block)*(p_well_block-p_bhp)*dTo_ds+(p_well_block-p_bhp)*dTg_ds);
    % dRg_dp_1=-((Rs_prime(info,p_well_block)*To+Rs(info,p_well_block)*dTo_dp)*(p_well_block-info.p_bhp)+Rs(info,p_well_block)*To+Tg+(p_well_block-info.p_bhp)*dTg_dp);
    % dRg_ds_1=-(Rs(info,p_well_block)*(p_well_block-info.p_bhp)*dTo_ds+(p_well_block-info.p_bhp)*dTg_ds);
    %
    %
    % % case 1 pressure control
    % l=(info.ywa-1)*info.Nx+info.xwa;
    % J_Q1(2*l-1,2*l-1)=dRg_dp_1;
    % J_Q1(2*l-1,2*l)=dRg_ds_1;
    % J_Q1(2*l,2*l-1)=dRo_dp_1;
    % J_Q1(2*l,2*l)=dRo_ds_1;
    %
    % % case 2 rate control
    % l=(info.ywa-1)*info.Nx+info.xwa;
    % J_Q2(2*l-1,2*l-1)=dRg_dp_2;
    % J_Q2(2*l-1,2*l)=dRg_ds_2;
    % J_Q2(2*l,2*l-1)=dRo_dp_2;
    % J_Q2(2*l,2*l)=dRo_ds_2;
    
end
% spy(J_T)
if control=='p'
    % case 1 pressure control
    J_matrix=J_T-J_D+J_Q1;
else
    % case 2 rate control
    J_matrix=J_T-J_D+J_Q2;
end
% J_matrix=J_T-J_D;
end