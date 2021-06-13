clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load masepdi_flight_basics.mat
load clouds.mat

%%

close all

thet_BL = [286;286;287;286];
thet_FB = [300;302;303;302];

qt_BL = [12.0;12.3;12.0;12.3];
qt_FT = [5.0 ;3.0 ;8.0 ;5.0 ];

%%
close all
iday=4;
figure('Position',[1732 209 794 630])
yyaxis right
plot(clouds.masepdi(iday).s_ps,'.')
yyaxis left
plot(clouds.masepdi(iday).s_qt)
grid('on')
%%

for iday = 1:length(clouds.masepdi)
    ti{iday} = masepdi_flight_basics(iday).ti;
    tf{iday} = masepdi_flight_basics(iday).tf;
    z_CB{iday} = masepdi_flight_basics(iday).z_CB;
    z_CT{iday} = masepdi_flight_basics(iday).z_CT;
    
    t = floor(clouds.masepdi(iday).s_t);
    lwc = clouds.masepdi(iday).s_lwc_pdi;
    z = clouds.masepdi(iday).s_ap;
    s_ntot = clouds.masepdi(iday).s_ntot_pdi;
    ta = clouds.masepdi(iday).s_ta;
    thet = clouds.masepdi(iday).s_thet;
    thete = clouds.masepdi(iday).s_thete;
    s_qt = clouds.masepdi(iday).s_qt;
    
    thet(thet<0) = nan;
    ta(ta<0) = nan;
    thete(thete<0) = nan;
    
%     thet_BL(iday,1) = nanmean(thet(z<nanmean(z_CB{iday})));
%     thet_FB(iday,1) = nanmean(thet(z<nanmean(z_CT{iday})*1.3 & z>nanmean(z_CT{iday})*1.2));
    
    T_BL(iday,1) = nanmean(ta(z<nanmean(z_CB{iday})));
    T_FB(iday,1) = nanmean(ta(z<nanmean(z_CT{iday})*1.3 & z>nanmean(z_CT{iday})*1.2));
    
%     qt_BL(iday,1) = nanmean(s_qt(z<nanmean(z_CB{iday})));
%     qt_FT(iday,1) = nanmean(s_qt(z<nanmean(z_CT{iday})*1.3 & z>nanmean(z_CT{iday})*1.2));
    
    ent_ratio_T = (thet-thet_BL(iday))./(thet_FB(iday)-thet_BL(iday));
    ent_ratio_qt = (s_qt-qt_BL(iday))./(qt_FT(iday)-qt_BL(iday));
    
%     figure
%     plot(t,lwc,'.')
%     yyaxis right
%     plot(t,z)
    
    ti_idx{iday} = arrayfun(@(x) findInSorted(t, ti{iday}(x)), 1:length(ti{iday}));
    tf_idx{iday} = arrayfun(@(x) findInSorted(t, tf{iday}(x)), 1:length(tf{iday}));
    
%     figure
%     plot(t, ent_ratio_T)
%     yyaxis right
%     plot(t, z)
%     plot(thet,z)
%     xline(thet_BL(iday));
%     xline(thet_FB(iday));
%     hold on
%     scatter(ti{iday}, z(ti_idx{iday}), 'g')
%     scatter(tf{iday}, z(tf_idx{iday}), 'b')
%     for ileg = 1:length(ti{iday})
%         z_leg = z(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         lwc_leg = lwc(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         t_leg = t(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         thet_leg = thet(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         thete_leg = thete(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         wz_leg = wz(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         s_ntot_leg = s_ntot(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%         
%         figure
%         plot(thete_leg,z_leg,'.')
%         title(['day ' num2str(iday) ', leg' num2str(ileg)])
%         CBline = refline(0, z_CB{iday}(ileg));
%         CTline = refline(0, z_CT{iday}(ileg));
%         
%         CBline.Color = 'r'; CBline.LineWidth = 1;
%         CTline.Color = 'g'; CTline.LineWidth = 1;
%         
%     end
    
    masepdi_flight_basics(iday).ti = ti{iday};
    masepdi_flight_basics(iday).tf = tf{iday};
    masepdi_flight_basics(iday).ti_idx = ti_idx{iday};
    masepdi_flight_basics(iday).tf_idx = tf_idx{iday};
    masepdi_flight_basics(iday).z_CB = z_CB{iday};
    masepdi_flight_basics(iday).z_CT = z_CT{iday};
    masepdi_flight_basics(iday).thet_BL = thet_BL(iday);
    masepdi_flight_basics(iday).thet_FB = thet_FB(iday);
    masepdi_flight_basics(iday).T_BL = T_BL(iday);
    masepdi_flight_basics(iday).T_FB = T_FB(iday);
    masepdi_flight_basics(iday).qt_BL = qt_BL(iday);
    masepdi_flight_basics(iday).qt_FB = qt_FT(iday);
%     
    clouds.masepdi(iday).ent_ratio_T = ent_ratio_T;
    clouds.masepdi(iday).ent_ratio_qt = ent_ratio_qt;
    clouds.masepdi(iday).s_ta = ta;
    clouds.masepdi(iday).s_thet = thet;
    clouds.masepdi(iday).s_thete = thete;
    
    %% pcasp validity testing
    a_t = clouds.masepdi(iday).a_t;
    a_ntot = clouds.masepdi(iday).a_ntot;
    ai_ntot = clouds.masepdi(iday).a_ntot_ex;
    
%     [cmt, cmt_ipdi, cmt_ipcasp] = intersect(t, a_t);
%     figure
%     plot(lwc(cmt_ipdi),ai_ntot(cmt_ipcasp),'.')
end
%% 
% save('masepdi_flight_basics.mat', 'masepdi_flight_basics')
save('clouds.mat','clouds', '-v7.3')