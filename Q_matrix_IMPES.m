function Q_matrix =Q_matrix_IMPES_Sat(info,P_matrix,control,d11_matrix,d12_matrix,d21_matrix,d22_matrix)
% producer well
[P_res_orig,S_res_orig]=to2D(info,P_matrix);    
p_well_block=P_res_orig(info.yw,info.xw);
s_well_block=S_res_orig(info.yw,info.xw);
ro=0.28*(((info.ky/info.kx)^(1/2)*info.delta_x^2+(info.kx/info.ky)^(1/2)*info.delta_y^2)^(1/2))/((info.ky/info.kx)^(1/4)+(info.kx/info.ky)^(1/4));
WI=(info.alpha*2*pi*(info.kx*info.ky)^(1/2)*info.Lz/(log(ro/info.rw)+info.s));
To=WI*(kro(s_well_block)/(info.mu_oil*Bo(info,p_well_block)));
qwo=-To*(p_well_block-info.p_bhp);
Tg=WI*(krg_new(s_well_block)/(mu_g(info,p_well_block)*Bg(info,p_well_block)));
qwg=-(Rs(info,p_well_block)*To*(p_well_block-info.p_bhp)+Tg*(p_well_block-info.p_bhp));
% Q_matrix=zeros(2*info.Nx*info.Ny,1);
Q_matrix=zeros(info.Nx*info.Ny,1);
l=(info.yw-1)*info.Nx+info.xw;
if control=='p'
% specify pressure
% Q_matrix(2*l-1,1)=qwg;
% Q_matrix(2*l,1)=qwo;
qg=qwg;
qo=qwo;
else
% specify flowrate
% Q_matrix(2*l-1,1)=-info.qwo_control*(Tg/To+Rs(info,p_well_block));
% Q_matrix(2*l,1)=-info.qwo_control;
qg=-info.qwo_control*(Tg/To+Rs(info,p_well_block));
qo=-info.qwo_control;
end
Q_matrix(l,1)=-(d22_matrix(info.yw,info.xw)/d12_matrix(info.yw,info.xw))*qg+qo;

% % injector well
% [P_res_orig,S_res_orig]=to2D(info,P_matrix);    
% p_well_block=P_res_orig(info.ywa,info.xwa);
% s_well_block=S_res_orig(info.ywa,info.xwa);
% ro=0.28*(((info.ky/info.kx)^(1/2)*info.delta_x^2+(info.kx/info.ky)^(1/2)*info.delta_y^2)^(1/2))/((info.ky/info.kx)^(1/4)+(info.kx/info.ky)^(1/4));
% WI=(info.alpha*2*pi*(info.kx*info.ky)^(1/2)*info.Lz/(log(ro/info.rw)+info.s));
% To=WI*(kro(s_well_block)/(info.mu_oil*Bo(info,p_well_block)));
% qwo=-To*(p_well_block-info.p_bhp);
% Tg=WI*(krg_new(s_well_block)/(mu_g(info,p_well_block)*Bg(info,p_well_block)));
% qwg=-(Rs(info,p_well_block)*To*(p_well_block-info.p_bhp)+Tg*(p_well_block-info.p_bhp));
% % Q_matrix=zeros(2*info.Nx*info.Ny,1);
% l=(info.ywa-1)*info.Nx+info.xwa;
% if control=='p'
% % specify pressure
% Q_matrix(2*l-1,1)=qwg;
% Q_matrix(2*l,1)=qwo;
% else
% % specify flowrate
% Q_matrix(2*l-1,1)=-info.qwo_control*(Tg/To+Rs(info,p_well_block));
% Q_matrix(2*l,1)=-info.qwo_control;
% end

end