clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load vocalspdi_flight_basics.mat
load clouds.mat

%%

thet_BL = [287.3;287.5;288;288.5;288.8;289;289.0;288;288.8;289.0;290.0;...
    290.2;290.5;290.0;290.0];
thet_FT = [309.0;308.3;311;311.0;311.5;310;308.5;309;309.3;308.5;310.5;...
    310.0;312.5;311.5;310.5];

qt_BL = [7.5;7.5;7.4;7.8;8.0;8.0;8.0;8.5;8.5;8.8;8.5;8.0;8.1 ;8.2;8.1];
qt_FT = [0.9;0.7;1.0;1.1;1.2;1.2;1.1;3.2;4.5;3.5;0.7;0.7;0.66;0.7;2.0];

%%
close all
iday=15;
figure('Position',[1732 209 794 630])
yyaxis right
plot(clouds.vocalspdi(iday).s_ap,'.')
yyaxis left
plot(clouds.vocalspdi(iday).s_qt)
grid('on')

%%


for iday = 1:length(clouds.vocalspdi)
    ti{iday} = vocalspdi_flight_basics(iday).ti;
    tf{iday} = vocalspdi_flight_basics(iday).tf;
    z_CB{iday} = vocalspdi_flight_basics(iday).z_CB;
    z_CT{iday} = vocalspdi_flight_basics(iday).z_CT;
    
    t = floor(clouds.vocalspdi(iday).s_t);
    lwc = clouds.vocalspdi(iday).s_lwc_pdi;
    z = clouds.vocalspdi(iday).s_ap;
    s_ntot = clouds.vocalspdi(iday).s_ntot_pdi;
    T = clouds.vocalspdi(iday).s_ta;
    thet = clouds.vocalspdi(iday).s_thet;
    thete = clouds.vocalspdi(iday).s_thete;
    wz = clouds.vocalspdi(iday).s_wz;
    s_qt = clouds.vocalspdi(iday).s_qt;
    Lv = 2.5e6; 
    cp = 1005;
    s_rhoa = clouds.vocalspdi(iday).s_rhoa;
    
    dthet_evap = Lv/cp*lwc/1000./s_rhoa;
%     
%     s_mr = clouds.vocalspdi(iday).s_mr/1000; % convert to unitless
    
%     
%     g = 9.8; 
%     
%     s_q = s_mr./(1+s_mr);
%     S = cp*ta + g*z + Lv*s_q; % calc moist static energy
    
%     thet_BL(iday,1) = mean(thet(z<mean(z_CB{iday})));
%     thet_FB(iday,1) = mean(thet(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
    
%     T_BL(iday,1) = mean(T(z<mean(z_CB{iday})));
%     T_FB(iday,1) = mean(T(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
    
%     qt_BL(iday,1) = mean(s_qt(z<mean(z_CB{iday})));
%     qt_FT(iday,1) = mean(s_qt(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
    
%     ent_ratio_thet = (thet-thet_BL(iday))./(thet_FT(iday)-thet_BL(iday));
%     ent_ratio_thet_corr = (thet+dthet_evap-thet_BL(iday))./(thet_FT(iday)-thet_BL(iday));
    
    ent_ratio_T = (thet-thet_BL(iday))./(thet_FT(iday)-thet_BL(iday));
    
    ent_ratio_qt = (s_qt-qt_BL(iday))./(qt_FT(iday)-qt_BL(iday));
    
%     figure
%     
%     plot(ent_ratio_qt, ent_ratio_thet,'.')
%     refline(1,0)
%     xlabel('z')
%     plot(z, T,'.')
%     ylabel('\theta')
%     yyaxis right
%     plot(z, s_qt,'.')
%     ylabel('q_t')

%     plot(t,lwc,'.')
%     yyaxis right
%     plot(t,z)
    
    ti_idx{iday} = arrayfun(@(x) findInSorted(t, ti{iday}(x)), 1:length(ti{iday}));
    tf_idx{iday} = arrayfun(@(x) findInSorted(t, tf{iday}(x)), 1:length(tf{iday}));
    
%     figure
%     plot(t, ent_ratio_T)
%     yyaxis right
%     plot(t, z)
%     plot(S,z)
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
    
%     vocalspdi_flight_basics(iday).ti = ti{iday};
%     vocalspdi_flight_basics(iday).tf = tf{iday};
%     vocalspdi_flight_basics(iday).ti_idx = ti_idx{iday};
%     vocalspdi_flight_basics(iday).tf_idx = tf_idx{iday};
%     vocalspdi_flight_basics(iday).z_CB = z_CB{iday};
%     vocalspdi_flight_basics(iday).z_CT = z_CT{iday};
%     vocalspdi_flight_basics(iday).thet_BL = thet_BL(iday);
%     vocalspdi_flight_basics(iday).thet_FB = thet_FB(iday);
%     vocalspdi_flight_basics(iday).T_BL = T_BL(iday);
%     vocalspdi_flight_basics(iday).T_FB = T_FB(iday);
%     vocalspdi_flight_basics(iday).qt_BL = qt_BL(iday);
%     vocalspdi_flight_basics(iday).qt_FB = qt_FB(iday);
%     
    clouds.vocalspdi(iday).ent_ratio_T = ent_ratio_T;
    clouds.vocalspdi(iday).ent_ratio_qt = ent_ratio_qt;
    %% pcasp validity testing
%     a_t = clouds.vocalspdi(iday).a_t;
%     a_ntot = clouds.vocalspdi(iday).a_ntot;
%     ai_ntot = clouds.vocalspdi(iday).a_ntot_ex;
%     
%     [cmt, cmt_ipdi, cmt_ipcasp] = intersect(t, a_t);
%     figure
%     plot(lwc(cmt_ipdi),ai_ntot(cmt_ipcasp),'.')
end

%%
% save('vocalspdi_flight_basics.mat', 'vocalspdi_flight_basics')
save('clouds.mat','clouds', '-v7.3')