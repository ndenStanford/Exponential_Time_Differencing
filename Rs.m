function Rs=Rs(info,P)
% P_bub=3400;%*6894.76; % psi to Pa
if P>info.P_bub
Rs=(178.11^2/5.615);    
else    
Rs=(178.11^2/5.615)*(P/info.P_bub)^1.3;
end
end