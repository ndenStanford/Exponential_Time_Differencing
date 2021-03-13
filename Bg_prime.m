function Bg_prime=Bg_prime(info,P)
% % P=P/6894.76; % Pa to psi 
% P_bub=3400; % psi 
% P_atm=14.7; % psi 

if P<info.P_bub
    Dp=(info.P_atm-P); % psi to Pa
    dDp_dP=-1; 
else
    Dp=(info.P_atm-info.P_bub); % psi to Pa
    dDp_dP=0; 
end

Bg_prime = (1.7e-3)*exp(1.7e-3*Dp)*dDp_dP;% psi-1 
end