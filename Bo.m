function Bo=Bo(info,P)
% % P=P/6894.76; % Pa to psi 
% P_bub=3400; % psi 
% P_atm=14.7; % psi 
% co=0.8e-5; % psi-1 to Pa-1

if P<info.P_bub
    Bo=exp(-8e-5*(info.P_atm-P)); % psi to Pa
else
    Bo=exp(-8e-5*(info.P_atm-info.P_bub))*exp(-info.co*(P-info.P_bub)); % psi to Pa
end

end