clear
load clouds.mat

%%
close all

for iday = 1:13
    s_t = clouds.masecas(iday).s_t;
    s_ap = clouds.masecas(iday).s_ap;
    
    s_thet = clouds.masecas(iday).s_thet;
    s_lwc_cas = clouds.masecas(iday).s_lwc_cas;
    s_lwc_gerb = clouds.masecas(iday).s_lwc_gerb;
    s_lwc_hotw = clouds.masecas(iday).s_lwc_hotw;
    a_ntot = clouds.masecas(iday).a_ntot;
    ccn_a = clouds.masecas(iday).ccn_a;
    ccn_b = clouds.masecas(iday).ccn_b;
    
    figure
%     hold on
%     plot(s_lwc_cas,s_lwc_gerb,'.')
%     plot(s_lwc_cas, s_lwc_hotw,'.')
%     legend('Gerber', 'Hot wire')
%     xlabel('CAS')
%     refline(1, 0)
%     ylim([0 1])
%     xlim([0 1])
% %     plot(s_t, s_lwc_hotw)
%     hold off
    
    subplot(2,1,1)
    plot(s_t, s_ap)
    yyaxis right
    plot(s_t, a_ntot, 'LineWidth', 2)
    ylabel('aerosol conc')
    
    subplot(2,1,2)
    plot(s_t, s_ap)
    yyaxis right
    plot(s_t, s_lwc_cas, 'LineWidth', 2)
    ylabel('LWC')
end


%%
ti = {[]; []; []; [62800 68970]; [62900 71510]; [61310 67430]; 65400; [61350 65530];...
    [61560 64360]; [61990 65710]; []; []; 61200}; 
tf = {[]; []; []; [67760 71980]; [69270 72590]; [66150 71030]; 71400; [63960 68060];...
    [63450 67800]; [62980 67690]; []; []; 67320};

for iday = 1:13
    masecas_flight_basics(iday).ti = ti{iday};
    masecas_flight_basics(iday).tf = tf{iday};
end

save('masecas_flight_basics.mat', 'masecas_flight_basics')