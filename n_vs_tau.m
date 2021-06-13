clear 

n = [21.3 76.9 201.2 564.6 1944.3];
tau_c = [58.1 17.6 8.0 4.2 1.4];

nlog10 = log10(n);
tau_clog10 = log10(tau_c);

% loglog(n, tau_c, 'LineWidth', 2)

plot(nlog10,tau_clog10,'LineWidth',2)

%%
p = polyfit(nlog10, tau_clog10, 1);

tau_c_hat = 10^p(2)*n.^(p(1));

plot(n, tau_c, n, tau_c_hat,'LineWidth',2)