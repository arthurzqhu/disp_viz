% nu = 1:1:60;
% D_eOverD = (nu+2)./nu;
% reldisp = sqrt(1./nu);
% plot(D_eOverD, reldisp)
x = 1:0.01:3;
y = sqrt((x-1)/2);
plot(x,y)
xlim([1 3])