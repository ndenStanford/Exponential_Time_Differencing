function mu_g_prime=mu_g_prime(info,P)
% P=P/6894.76; % Pa tp psi
mu_g_prime = ((3e-10)*2*P + 1e-6);%*10^-3; % cp tp Pa*s
end