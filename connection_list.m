% make connection list
function[Px_upwinded,Py_upwinded,Sx_upwinded,Sy_upwinded,Hx_oil,Hy_oil,Hx_gas,Hy_gas]=connection_list(info,P_matrix)
[P_res_orig,S_res_orig]=to2D(info,P_matrix);

[row,col]=size(P_res_orig);
% Hx
Px_upwinded=zeros(row,col-1);
Sx_upwinded=zeros(row,col-1);
for i=1:row
    for j=1:col-1
        if (P_res_orig(i,j)>P_res_orig(i,j+1))
            Px_upwinded(i,j)=P_res_orig(i,j);
            Sx_upwinded(i,j)=S_res_orig(i,j);
        else
            Px_upwinded(i,j)=P_res_orig(i,j+1);
            Sx_upwinded(i,j)=S_res_orig(i,j+1);
        end        
    end
end
% Hy
Py_upwinded=zeros(row-1,col);
Sy_upwinded=zeros(row-1,col);
for i=1:row-1
    for j=1:col
        if (P_res_orig(i,j)>P_res_orig(i+1,j))
            Py_upwinded(i,j)=P_res_orig(i,j);
            Sy_upwinded(i,j)=S_res_orig(i,j);
        else
            Py_upwinded(i,j)=P_res_orig(i+1,j);
            Sy_upwinded(i,j)=S_res_orig(i+1,j);
        end        
    end
end


Hx_oil=kro(Sx_upwinded)./(info.mu_oil*Bo(info,Px_upwinded));
Hy_oil=kro(Sy_upwinded)./(info.mu_oil*Bo(info,Py_upwinded));
Hx_gas=krg_new(Sx_upwinded)./(mu_g(info,Px_upwinded).*Bg(info,Px_upwinded));
Hy_gas=krg_new(Sy_upwinded)./(mu_g(info,Py_upwinded).*Bg(info,Py_upwinded));
end