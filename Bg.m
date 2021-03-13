function Bg=Bg(info,P)
% % P=P/6894.76; % Pa to psi 
% P_bub=3400; % psi 
% P_atm=14.7; % psi 

if P<info.P_bub
    Dp=(info.P_atm-P); % psi to Pa
else
    Dp=(info.P_atm-info.P_bub); % psi to Pa
end

Bg = exp(1.7e-3*Dp); % cp tp Pa*s
end