c=3;

camps={'vocalspdi','masepdi','postpdi'};
campaign=camps{c};

fb = load([campaign,'_flight_basics.mat']);
fbvar = [campaign,'_flight_basics'];

iday = 16;


jleg = 1;
                
t = clouds.(campaign)(iday).s_t;
z = clouds.(campaign)(iday).s_ap;
lwc = clouds.(campaign)(iday).s_lwc_pdi;
rhoa = clouds.(campaign)(iday).s_rhoa;
ql_obs = lwc./rhoa;
ql_adb = clouds.(campaign)(iday).ql_adb_prof;

ti = fb.(fbvar)(iday).ti(jleg);
tf = fb.(fbvar)(iday).tf(jleg);

% ti = 4.9e4;
% tf = 4.92e4;

z_CB = fb.(fbvar)(iday).z_CB(jleg);
z_CT = fb.(fbvar)(iday).z_CT(jleg);

idx_leg = t>=ti & t<=tf;
% idx_leg = 1:length(z);

% plot(t,z)
plot(ql_obs(idx_leg),z(idx_leg),'.'); hold on
plot(ql_adb(idx_leg),z(idx_leg)); hold off


% plot(ql_adb_lin,z_lin-40); hold off
% plot(ql_adb_prof,z)
% ylim([z_CB z_CT])