% function C = adiab_ql_aigel(z_CB, z_CT, thet_CB, qv_CB, p_CB)

z_CB = vocalspdi_flight_basics(iday).z_CB;
z_CT = vocalspdi_flight_basics(iday).z_CT;
thet_CB = nanmean(gen(iday).s_thet(gen(iday).normAC<0.05 & gen(iday).normAC>-0.05));
qv_CB = nanmean(gen(iday).s_mr(gen(iday).normAC<0.05 & gen(iday).normAC>-0.05));
p_CB = nanmean(gen(iday).s_ps(gen(iday).normAC<0.05 & gen(iday).normAC>-0.05));

C(1) = 0;
thetap(1)= thet_CB; % cloud base potential temperature
qv(1) = qv_CB; % cloud base water vapor mixing ratio
kappa = 287/1004;

dz = 0.1;

z = z_CB:dz:z_CT;
H = 8100;

p0 = 1000;

p = p_CB*exp(-(z-z_CB)/H);
T(1) = thetap(1)*(p0/p(1))^(-kappa);

for k=2:length(z)
    
	thetap(k)=thetap(k-1);
    T(k)=T(k-1);
	qv(k)=qv(k-1);
	qs(k)=wsat_tp(T(k)-273.15,p(k))*1000;
    
	phi=qs(k)*(17.27*237*2.536/(1004*(T(k)-36)^2));
	C(k)=(qv(k)-qs(k))/(1+phi); %condensed liquid
	thetap(k)=thetap(k)+2.5e6*C(k)/(1004*(p0/p(k))^(-kappa)); %potential temperature update
	qv(k)=qv(k)-C(k); %water vapor update

end

cumul_cond_liq = cumsum(C);