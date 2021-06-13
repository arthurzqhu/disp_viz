function ql_adb = adiab_ql(z_CB,z_input,T_CB,r_CB,p_CB)

% z_CB = 710;
% z_input = 1025;
% T_CB = 8.72+273.15;
% r_CB = 7.36/1000;
% p_CB = 934;

% T in Kelvin
% r in kg/kg
% p in mb

H = 8100;
g = 9.8;
cp = 1004;
Lv = 2.5e6;
Rd = 287;
Rv = 462;

dz = 0.1;

z = z_CB:dz:z_input;

p = p_CB*exp(-(z-z_CB)/H);

T(1) = T_CB;
r(1) = r_CB;

for ialt = 1:length(z)-1
    
    dTdz(ialt) = -g*(1+(Lv*r(ialt))/(Rd*T(ialt)))/(cp + Lv^2*r(ialt)/(Rv*T(ialt)^2));
    drdT(ialt) = dwsatdT_tp(T(ialt)-273.15,p(ialt));
    
    T(ialt+1) = T(ialt) + dTdz(ialt)*dz;
    r(ialt+1) = r(ialt) + drdT(ialt)*dTdz(ialt)*dz;
    
end

q_CB = r_CB/(r_CB+1);
q = r./(r+1);

ql_adb = q_CB-q;

% wsat = wsat_tp(0:0.01:25,934);
% dT = 0.01;
% dwsatdT = (wsat(2:end)-wsat(1:end-1))/dT;